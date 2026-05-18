# model.py -- GWSWEXmodel: end-user model interface.
from __future__ import annotations

import types
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from .config import (
    InitialConditions,
    LateralFluxes,
    Material,
    ModelParams,
    SolverConfig,
    SpatialDomain,
    TemporalDomain,
    Vegetation,
    _broadcast_to_ne,
    _broadcast_to_nl_ne,
)
from .io import GwswexNCWriter

# f2py compiled module (import deferred to init to allow pre-config)
_F: Optional[types.ModuleType] = None  # set to the f2py module at kernel init time

_UNIT_SCALES_T = {"s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}
_UNIT_SCALES_L = {"m": 1.0, "cm": 0.01, "mm": 0.001}


class GWSWEXmodel:
    """
    GWSWEX model interface.

    Parameters
    ----------
    name : str
        Model run name (used to derive the default NetCDF output filename).
    T : str
        Time unit for user-facing I/O ('s', 'min', 'h', 'd').
    L : str
        Length unit for user-facing I/O ('m', 'cm', 'mm').
    write_output : bool, default True
        If True, a NetCDF output file is opened during :meth:`init` and each
        call to :meth:`run_step` / :meth:`run` writes a timestep to it.  The
        file is closed on :meth:`deinit`.
    output_fpath : str | pathlib.Path | None, default None
        Path to the NetCDF output file.  When ``None`` (default) and
        ``write_output=True``, the file is written to the current working
        directory using a filename derived from ``name`` via
        :meth:`_default_output_filename`.
    """

    # Characters stripped from the model name when deriving a default output
    # filename.  Covers filesystem-unsafe characters on Windows/macOS/Linux
    # plus single/double quotes.
    _FILENAME_UNSAFE = set('/\\:*?"<>|\'')

    @classmethod
    def _default_output_filename(cls, name: str) -> str:
        """Sanitise *name* into a safe lowercase ``.nc`` filename (<=64 chars)."""
        s = name.strip().lower()
        s = "".join("_" if c.isspace() else c for c in s)
        s = "".join(c for c in s if c not in cls._FILENAME_UNSAFE)
        # collapse runs of underscores and trim leading/trailing underscores
        while "__" in s:
            s = s.replace("__", "_")
        s = s.strip("_.")
        if not s:
            s = "gwswex"
        return s[:64] + ".nc"

    def __init__(
        self,
        name: str = "gwswex",
        T: str = "s",
        L: str = "m",
        write_output: bool = True,
        output_fpath: Optional[str | Path] = None,
        flush_nc: bool = False,
    ):
        self.name = name
        self._T_scale = _UNIT_SCALES_T[T.lower()]
        self._L_scale = _UNIT_SCALES_L[L.lower()]
        self._T_unit = T.lower()
        self._L_unit = L.lower()

        # NetCDF output
        self._write_output: bool = bool(write_output)
        self._flush_nc: bool = bool(flush_nc)
        if output_fpath is not None:
            self._output_fpath: Optional[Path] = Path(output_fpath)
        elif self._write_output:
            self._output_fpath = Path.cwd() / self._default_output_filename(name)
        else:
            self._output_fpath = None

        # configuration slots (set via register_* methods)
        self.space: Optional[SpatialDomain] = None
        self.time: Optional[TemporalDomain] = None
        self.solver: SolverConfig = SolverConfig()
        self.model_params: ModelParams = ModelParams()
        self.ic: Optional[InitialConditions] = None
        self._lateral: Optional[LateralFluxes] = None

        self._is_initialised = False
        self._writer: Optional[GwswexNCWriter] = None
        self._vegetation: dict[int, Vegetation] = {}

        # root growth state (set by _register_space at init time)
        self._has_dynamic_growth: bool = False
        self._root_depth_initial_per_elem: Optional[np.ndarray] = None  # [ne] user L
        self._root_depth_final_per_elem: Optional[np.ndarray] = None  # [ne] user L

        # stored forcing (set by set_forcing, consumed by run / run_step)
        self._forcing: Optional[dict] = None

        # per-step mass balance history (populated during run / run_step)
        self._mass_balance_history: list[dict] = []

    # ------------------------------------------------------------------
    # Configuration registration
    # ------------------------------------------------------------------
    def add_material(self, **kwargs) -> None:
        """Register a soil material. Must be called before init()."""
        mat = Material(**kwargs)
        if self.space is None:
            raise RuntimeError("Call init_space before adding materials")
        self.space.add_material(mat)

    def add_vegetation(self, **kwargs) -> None:
        """Register a vegetation type. id and root are required; et_stress, root_growth are optional."""
        veg = Vegetation(**kwargs)
        self._vegetation[veg.id] = veg

    def init_space(self, **kwargs) -> None:
        """Create (unfrozen) spatial domain."""
        self.space = SpatialDomain(**kwargs)

    def _register_space(self) -> None:
        """Validate the spatial configuration and derive the rooting mask.

        Idempotent: safe to call multiple times before :meth:`init`.
        Does *not* freeze the spatial domain; freezing is performed by
        :meth:`init` after all wrapper calls succeed.
        """
        if self.space is None:
            raise RuntimeError("Call init_space before init")

        # check all sID reference registered materials
        for sid in np.unique(self.space.sID):
            if sid not in self.space._materials:
                raise ValueError(f"Material ID {sid} referenced in sID but not registered")

        # validate vID references
        for vid in np.unique(self.space.vID):
            if vid not in self._vegetation:
                raise ValueError(f"Vegetation ID {vid} referenced in vID but not registered")

        # Derive is_root from vegetation parameters. Transpiration demand is
        # apportioned uniformly across rooted layers at solve time, so no
        # per-layer weight array is built here.
        nl, ne = self.space.nl, self.space.ne
        is_root = np.zeros((nl, ne), dtype=np.int32)
        has_dynamic_growth: bool = False
        rd_init_arr = np.zeros(ne, dtype=float)
        rd_final_arr = np.zeros(ne, dtype=float)

        for ex in range(ne):
            vid = int(self.space.vID[ex])
            veg = self._vegetation[vid]
            assert self.space.bnds is not None, "space.bnds is None; spatial domain may be incomplete"
            surface = float(self.space.bnds[ex, 0])

            if veg.root_growth_model != "static" and veg.root_depth_initial is not None:
                # Depth-based dynamic growth.
                # At t=0 roots occupy layers whose midpoint depth <= root_depth_initial;
                # at t=end they occupy layers whose midpoint depth <= root_depth_final.
                has_dynamic_growth = True
                d_init = float(veg.root_depth_initial)
                d_final = float(veg.root_depth_final) if veg.root_depth_final is not None else d_init
                d_init_si = d_init / self._L_scale
                rd_init_arr[ex] = d_init
                rd_final_arr[ex] = d_final
                for l in range(nl):
                    z_mid = 0.5 * (self.space.bnds[ex, l] + self.space.bnds[ex, l + 1])
                    depth_m = (surface - z_mid) / self._L_scale
                    # Initial mask reflects the t=0 rooting state; subsequent
                    # macro-steps push the current mask through update_is_root.
                    if 0.0 <= depth_m <= d_init_si:
                        is_root[l, ex] = 1
            else:
                # Static specification via RootParams(depth=...)
                if veg.root is None:
                    raise ValueError(
                        f"Vegetation ID {veg.id}: 'root' must be provided for static "
                        "vegetation, or root_depth_initial/final with a non-static "
                        "root_growth_model must be supplied."
                    )
                root_depth = float(veg.root.depth)
                rd_init_arr[ex] = root_depth
                rd_final_arr[ex] = root_depth
                for l in range(nl):
                    z_mid = 0.5 * (self.space.bnds[ex, l] + self.space.bnds[ex, l + 1])
                    depth = surface - z_mid
                    if 0.0 <= depth <= root_depth:
                        is_root[l, ex] = 1

        self.space.is_root = is_root
        self._has_dynamic_growth = has_dynamic_growth
        self._root_depth_initial_per_elem = rd_init_arr
        self._root_depth_final_per_elem = rd_final_arr

    def register_space(self) -> None:
        """Deprecated: registration now happens automatically inside :meth:`init`.

        Retained as a no-op shim that runs the validation eagerly so that
        existing user scripts and notebooks continue to work. The spatial
        domain remains unfrozen until :meth:`init` succeeds.
        """
        import warnings

        warnings.warn(
            "register_space() is deprecated and is no longer required; "
            "init() now performs registration automatically.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._register_space()

    def init_time(self, start=None, stop=None, dt=None, dt_min=None, adaptive=True, n_steps=None) -> None:
        """Create the temporal domain.

        dt and dt_min may be numeric (in user T units) or datetime.timedelta.
        start and stop may be datetime.datetime or numeric; n_steps is derived
        from (stop - start) / dt when not supplied explicitly.
        """
        from datetime import datetime as _dt
        from datetime import timedelta as _td

        # Convert dt
        if isinstance(dt, _td):
            dt_user = dt.total_seconds() / self._T_scale
        elif dt is not None:
            dt_user = float(dt)
        else:
            raise ValueError("dt must be provided to init_time()")

        # Convert dt_min
        if isinstance(dt_min, _td):
            dt_min_user = dt_min.total_seconds() / self._T_scale
        elif dt_min is not None:
            dt_min_user = float(dt_min)
        else:
            dt_min_user = 1e-6

        # Compute n_steps
        if n_steps is None:
            if start is not None and stop is not None:
                if isinstance(start, _dt) and isinstance(stop, _dt):
                    duration_s = (stop - start).total_seconds()
                else:
                    duration_s = (float(stop) - float(start)) * self._T_scale
                n_steps = int(round(duration_s / (dt_user * self._T_scale)))
            else:
                raise ValueError("Either n_steps or both start and stop must be provided to init_time()")

        self.time = TemporalDomain(dt=dt_user, dt_min=dt_min_user, n_steps=n_steps, adaptive=adaptive)

    def _register_time(self) -> None:
        """Validate the temporal domain. Idempotent; does not freeze."""
        assert self.time is not None, "Call init_time() before init()"

    def register_time(self) -> None:
        """Deprecated: registration now happens automatically inside :meth:`init`."""
        import warnings

        warnings.warn(
            "register_time() is deprecated and is no longer required; "
            "init() now performs registration automatically.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._register_time()

    def set_initial_conditions(self, gw, sw, uz) -> None:
        """Set initial conditions with automatic shape broadcasting.

        Parameters
        ----------
        gw : scalar or array-like
            Initial GW table elevations [user L], broadcastable to ``(ne,)``.
        sw : scalar or array-like
            Initial SW depths [user L], broadcastable to ``(ne,)``.
        uz : scalar or array-like
            Initial UZ pressure heads [user L], broadcastable to ``(nl, ne)``.
            Provide in ``(ne, nl)`` orientation (row = element) or as a scalar
            sentinel ``-999`` to initialise at the equilibrium profile.
        """
        if self.space is None:
            raise RuntimeError("Call init_space() before set_initial_conditions()")
        ne, nl = self.space.ne, self.space.nl
        gw_arr = _broadcast_to_ne(gw, ne, "gw").astype(np.float64)
        sw_arr = _broadcast_to_ne(sw, ne, "sw").astype(np.float64)
        uz_arr = _broadcast_to_nl_ne(uz, nl, ne, "uz").astype(np.float64)
        self.ic = InitialConditions(gw=gw_arr, sw=sw_arr, uz=uz_arr)

    def _register_initial_conditions(self) -> None:
        """Validate the initial conditions. Idempotent; does not freeze."""
        assert self.ic is not None, "Call set_initial_conditions() before init()"

    def register_initial_conditions(self) -> None:
        """Deprecated: registration now happens automatically inside :meth:`init`."""
        import warnings

        warnings.warn(
            "register_initial_conditions() is deprecated and is no longer required; "
            "init() now performs registration automatically.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._register_initial_conditions()

    def set_solver(self, **kwargs) -> None:
        if self._is_initialised:
            raise RuntimeError(
                "set_solver() must be called before init(). " "Deinit the model and re-initialise to change the solver."
            )
        self.solver = SolverConfig(**kwargs)

    def switch_solver(
        self,
        *,
        warm_start: str = "proxy",
        icratio_init: float | None = None,
        f_ga_init: float | None = None,
        **kwargs,
    ) -> None:
        """Switch the active solver (and its parameters) mid-simulation.

        Translates the in-kernel state from the current solver's
        representation into the new one and updates the solver-specific
        parameters (Picard tolerance, Courant number, etc.).  Safe to call
        between :meth:`run_step` invocations and across a checkpoint
        restart (see :meth:`get_checkpoint` / :meth:`load_checkpoint`).

        Parameters
        ----------
        warm_start : {"proxy", "cold", "manual"}, default "proxy"
            Selects how the explicit solver's persistent Green-Ampt /
            connectivity state is seeded after a switch *to* explicit;
            ignored when the destination is implicit.

            - ``"proxy"`` (default): per-layer ``IC`` and ``ICratio`` are
              derived from the converged effective saturation
              :math:`S_e = (\\theta - \\theta_r)/(\\theta_s - \\theta_r)`
              of the implicit profile; per-element ``F_GA`` is the
              column-integral of :math:`(\\theta - \\theta_r)\\Delta z`,
              capped at :math:`5\\,\\psi_f` so the next-step
              Green-Ampt infiltration capacity does not collapse.
            - ``"cold"``: leave the GA state at the cold defaults
              applied by ``kernel_switch_solver`` (``IC = 0``,
              ``ICratio = ICratio_min``, ``F_GA = F_min``).
            - ``"manual"``: apply the user-supplied scalars
              ``icratio_init`` (uniform, ``IC = icratio_init * d_a``)
              and ``f_ga_init`` (uniform per element).  At least one of
              the two must be provided; the unspecified field falls
              back to the cold default.
        icratio_init, f_ga_init : float, optional
            Only meaningful with ``warm_start="manual"``.  See above.
        **kwargs : dict
            Same fields accepted by :class:`SolverConfig`.  ``solver``
            (``"explicit"`` or ``"implicit"``) is required and selects the
            target solver; all other fields default to the same values
            applied at the original :meth:`init`.

        Notes
        -----
        - When switching to ``"implicit"``, the per-layer matric head
          warm start is derived from the current ``theta_curr`` profile
          via the analytical van-Genuchten inverse so that no spurious
          storage transient is introduced on the first Picard solve.
          The ``warm_start`` argument is ignored in this direction.
        - When switching to ``"explicit"``, the persistent Green-Ampt
          and connectivity state is reset to the cold defaults inside
          the kernel translator, then optionally overwritten according
          to ``warm_start``.
        - The OMP-thread count and model parameters set via
          :meth:`set_model_params` are preserved across the switch.
        """
        if warm_start not in ("proxy", "cold", "manual"):
            raise ValueError(f"warm_start must be 'proxy', 'cold', or 'manual'; got {warm_start!r}")
        if not self._is_initialised:
            raise RuntimeError(
                "switch_solver() requires the model to be initialised. "
                "Use set_solver() before init() for the initial choice."
            )
        global _F
        new_cfg = SolverConfig(**kwargs)

        ierr = _F.gwswex_wrapper.switch_solver(new_cfg.solver_type_id)
        if ierr != 0:
            raise RuntimeError(f"Kernel solver switch failed (ierr={ierr})")

        # Push the new solver-specific parameters into the kernel.
        _F.gwswex_wrapper.set_solver_params(
            new_cfg.courant_number,
            self.time.dt_min * self._T_scale,
            new_cfg.beta_hyst,
            new_cfg.n_trapz,
            new_cfg.h_min * self._L_scale,
        )
        if new_cfg.solver == "implicit":
            _F.gwswex_wrapper.set_picard_params(
                new_cfg.picard_tol * self._L_scale,
                new_cfg.picard_max_iter,
            )
        # Warm-start of the explicit Green-Ampt / connectivity state.
        # Only meaningful when the destination is the explicit solver.
        if new_cfg.solver == "explicit":
            if warm_start == "proxy":
                if icratio_init is not None or f_ga_init is not None:
                    raise ValueError("icratio_init / f_ga_init are only honoured with " "warm_start='manual'")
                ierr_ws = _F.gwswex_wrapper.warm_start_explicit_proxy()
                if ierr_ws != 0:
                    raise RuntimeError(f"Explicit proxy warm-start failed (ierr={ierr_ws})")
            elif warm_start == "manual":
                if icratio_init is None and f_ga_init is None:
                    raise ValueError("warm_start='manual' requires at least one of " "icratio_init or f_ga_init")
                r_in = float(icratio_init) if icratio_init is not None else -1.0
                f_in = float(f_ga_init) * self._L_scale if f_ga_init is not None else -1.0
                ierr_ws = _F.gwswex_wrapper.warm_start_explicit(r_in, f_in)
                if ierr_ws != 0:
                    raise RuntimeError(f"Explicit manual warm-start failed (ierr={ierr_ws})")
            # warm_start == "cold": kernel translator already applied the
            # cold defaults; nothing further to do.
        # OMP threads are preserved across the switch (unchanged in kernel).
        new_cfg.omp_threads = self.solver.omp_threads

        # Mirror into the Python-side config object.
        self.solver.unfreeze()
        self.solver = new_cfg
        self.solver.freeze()

    def set_omp_threads(self, n: int) -> None:
        """Set the OpenMP thread count for the element-parallel loop.

        Unlike set_solver(), this may be called at any time — before or after
        init() — because changing the OMP thread count requires no
        reallocation and does not alter the model state.  The new value is
        applied immediately to the OpenMP runtime and mirrored into the
        in-kernel model singleton.  When called before init(), the value is
        also stored in self.solver so that a subsequent init() call sees the
        same setting.

        Has no effect on the OMP runtime when the binary was built without
        -fopenmp (thread count remains 1); the value is stored regardless.
        """
        if n < 1:
            raise ValueError(f"omp_threads must be >= 1, got {n}")
        if not self._is_initialised:
            # Pre-init path: solver is not yet frozen.
            self.solver.omp_threads = n
        else:
            # Post-init path: update the kernel runtime directly.
            global _F
            _F.gwswex_wrapper.set_omp_threads(n)
            # Mirror into the Python-side config for introspection / serialisation.
            self.solver.unfreeze()
            self.solver.omp_threads = n
            self.solver.freeze()

    def set_model_params(self, **kwargs) -> None:
        self.model_params = ModelParams(**kwargs)

    # ------------------------------------------------------------------
    # Kernel lifecycle
    # ------------------------------------------------------------------
    def init(self) -> None:
        """Initialise the Fortran kernel with all registered configuration."""
        global _F
        from . import f_gwswex as _F  # type: ignore[attr-defined]  # f2py .so has no stub

        assert _F is not None, "f2py kernel module failed to load"

        assert self.space is not None, "Call init_space() before init()"
        assert self.time is not None, "Call init_time() before init()"
        assert self.ic is not None, "Call set_initial_conditions() before init()"

        # Run validation / derivation steps automatically. These are
        # idempotent and safe to invoke even if the user (legacy code) has
        # already called the deprecated public register_* shims.
        self._register_space()
        self._register_time()
        self._register_initial_conditions()

        ne = self.space.ne
        sp = self.space._to_fortran()

        # Convert dimensional quantities to SI (metres, seconds).
        # _L_scale converts user L → metres; _T_scale converts user T → seconds.
        rate_scale = self._L_scale / self._T_scale  # [user LT⁻¹] → [m s⁻¹]
        L = self._L_scale  # [user L]    → [m]
        L_inv = 1.0 / self._L_scale if self._L_scale != 0 else 1.0  # [user L⁻¹] → [m⁻¹]

        ierr = _F.gwswex_wrapper.init(
            self.solver.solver_type_id,
            sp["bnds"] * L,  # layer boundaries [L] → [m]
            sp["sID"],
            sp["K_sat"] * rate_scale,  # K_sat [LT⁻¹] → [m s⁻¹]
            sp["theta_s"],  # [-] dimensionless
            sp["theta_r"],  # [-] dimensionless
            sp["alpha"] * L_inv,  # [L⁻¹] → [m⁻¹]
            sp["vg_n"],  # [-] dimensionless
            sp["lam"],  # [-] dimensionless
            sp["is_root"],
        )
        if ierr != 0:
            raise RuntimeError(f"Kernel init failed (ierr={ierr})")

        # set solver and ET params
        _F.gwswex_wrapper.set_solver_params(
            self.solver.courant_number,
            self.time.dt_min * self._T_scale,
            self.solver.beta_hyst,
            self.solver.n_trapz,
            self.solver.h_min * self._L_scale,  # [L] → [m]
        )
        _F.gwswex_wrapper.set_omp_threads(self.solver.omp_threads)
        if self.solver.solver == "implicit":
            _F.gwswex_wrapper.set_picard_params(
                self.solver.picard_tol * self._L_scale,  # [L] → [m]
                self.solver.picard_max_iter,
            )
        _F.gwswex_wrapper.set_model_params(
            self.model_params.psi_f * self._L_scale,  # [L] → [m]
            self.model_params.F_min * self._L_scale,  # [L] → [m]
            self.model_params.ICratio_min,  # [-] dimensionless
        )

        # set initial conditions (convert units to SI; preserve -999 sentinel)
        gw = self.ic.gw * self._L_scale
        sw = self.ic.sw * self._L_scale
        uz_raw = self.ic.uz.copy()
        uz_sentinel = uz_raw == -999.0
        uz = uz_raw * self._L_scale
        uz[uz_sentinel] = -999.0  # restore sentinel after unit conversion
        ierr = _F.gwswex_wrapper.set_ic(gw, sw, uz)
        if ierr != 0:
            raise RuntimeError(f"IC setup failed (ierr={ierr})")

        # store vegetation library in Fortran kernel so it can own root distribution
        self._call_set_vegetation(_F, ne)

        # All wrapper calls succeeded: freeze every configuration object now
        # so that subsequent mutations are caught. Freezing is deferred to
        # this point so that user-side helpers and the deprecated
        # register_* shims can run repeatedly during setup without
        # tripping the frozen guard.
        self.space.freeze()
        self.time.freeze()
        self.ic.freeze()
        self.solver.freeze()
        self.model_params.freeze()

        self._is_initialised = True

        # Open the built-in NetCDF writer if requested.  Opened here
        # so that it is ready for per-step writes by run_step(); closed in
        # deinit().  The caller's explicit run(output_file=...) path opens
        # its own writer and closes it at the end of run(); the two paths are
        # mutually exclusive per instance.
        if self._write_output and self._writer is None and self._output_fpath is not None:
            self._output_fpath.parent.mkdir(parents=True, exist_ok=True)
            self._writer = GwswexNCWriter(str(self._output_fpath), self, flush_each_step=self._flush_nc)
            self._writer.write_config()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _call_set_vegetation(self, _F_mod, ne: int) -> None:
        """Package Python vegetation objects and call wrapper.set_vegetation().

        For static vegetation types (``veg.root`` set), the rooting depth is
        forwarded to Fortran, which recomputes the per-element ``is_root`` mask
        from the column geometry. For dynamic-growth types, ``root_depth`` is
        passed as 0.0 so that the Python-supplied mask (carrying every layer
        that is ever rooted during the run) is preserved on the kernel side;
        the per-step mask is updated through :meth:`update_is_root` driven by
        :meth:`_recompute_is_root_at_frac`.
        """
        assert self.space is not None
        nveg_ids = sorted(self._vegetation.keys())
        nveg = len(nveg_ids)

        vID_arr = np.array([int(self.space.vID[ex]) for ex in range(ne)], dtype=np.int32)
        root_depth = np.zeros(nveg, dtype=np.float64)
        s_star_arr = np.zeros(nveg, dtype=np.float64)
        s_w_arr = np.zeros(nveg, dtype=np.float64)
        s_h_arr = np.zeros(nveg, dtype=np.float64)
        s_e_arr = np.zeros(nveg, dtype=np.float64)

        for i, vid in enumerate(nveg_ids):
            veg = self._vegetation[vid]
            if veg.root_growth_model == "static" and veg.root is not None and veg.root.depth > 0.0:
                root_depth[i] = veg.root.depth * self._L_scale
            # else: root_depth stays 0.0 → Fortran keeps Python-supplied mask
            ets = veg.et_stress
            s_star_arr[i] = ets.s_star
            s_w_arr[i] = ets.s_w
            s_h_arr[i] = ets.s_h
            s_e_arr[i] = ets.s_e

        ierr = _F_mod.gwswex_wrapper.set_vegetation(
            vID_arr,
            root_depth,
            s_star_arr,
            s_w_arr,
            s_h_arr,
            s_e_arr,
        )
        if ierr != 0:
            raise RuntimeError(f"Vegetation setup failed (ierr={ierr})")

    def _recompute_is_root_at_frac(self, frac: float) -> np.ndarray:
        """Recompute the rooting mask at an interpolated root depth.

        For dynamic-growth vegetation, the current rooting depth at
        ``frac \u2208 [0, 1]`` is obtained by linear interpolation between the
        per-element initial and final depths; a layer is rooted whenever its
        midpoint depth lies within the current rooting depth.
        """
        assert self.space is not None
        assert self.space.bnds is not None
        assert self._root_depth_initial_per_elem is not None
        assert self._root_depth_final_per_elem is not None

        ne, nl = self.space.ne, self.space.nl
        is_root = np.zeros((nl, ne), dtype=np.int32)

        for ex in range(ne):
            rd_i = self._root_depth_initial_per_elem[ex]
            rd_f = self._root_depth_final_per_elem[ex]
            rd_curr = rd_i + (rd_f - rd_i) * frac  # current root depth [user L]
            surface = float(self.space.bnds[ex, 0])

            for l in range(nl):
                z_mid = 0.5 * (self.space.bnds[ex, l] + self.space.bnds[ex, l + 1])
                depth_l = surface - z_mid
                if 0.0 <= depth_l <= rd_curr:
                    is_root[l, ex] = 1

        return is_root

    def deinit(self) -> None:
        """Clean up Fortran kernel memory."""
        if self._is_initialised:
            assert _F is not None, "Kernel module not loaded"
            _F.gwswex_wrapper.deinit()
            self._is_initialised = False
        if self._writer is not None:
            self._writer.close()
            self._writer = None

    # ------------------------------------------------------------------
    # Stepping
    # ------------------------------------------------------------------
    def set_lateral(self, gw: np.ndarray, sw: np.ndarray) -> None:
        """Set lateral flux rates [user L/T units] for the next step."""
        self._lateral = LateralFluxes(
            gw=gw * (self._L_scale / self._T_scale),
            sw=sw * (self._L_scale / self._T_scale),
        )

    def step(self, dt: float, precip: np.ndarray, pet: np.ndarray, ptt: np.ndarray) -> None:
        """Execute one macro time-step."""
        if not self._is_initialised:
            raise RuntimeError("Call model.init() before stepping")
        assert _F is not None
        assert self.space is not None

        dt_si = dt * self._T_scale
        rate_scale = self._L_scale / self._T_scale

        lat_gw = self._lateral.gw if self._lateral else np.zeros(self.space.ne)
        lat_sw = self._lateral.sw if self._lateral else np.zeros(self.space.ne)

        ierr = _F.gwswex_wrapper.step(
            dt_si,
            precip * rate_scale,
            pet * rate_scale,
            ptt * rate_scale,
            lat_gw,
            lat_sw,
        )
        if ierr != 0:
            raise RuntimeError(f"Kernel step failed (ierr={ierr})")

        self._lateral = None  # consume lateral fluxes

    def set_forcing(
        self,
        precip,
        pet,
        ptt,
        lat_gw=None,
        lat_sw=None,
    ) -> None:
        """Store atmospheric and lateral forcing for use by :meth:`run` / :meth:`run_step`.

        Each forcing array is coerced to shape ``(n_steps, ne)`` [user L/T units].
        Scalars and 1-D arrays (per-step or per-element) are broadcast accordingly.

        Call this after :meth:`init_time` and :meth:`init_space`.
        """
        if self.space is None or self.time is None:
            raise RuntimeError("Call init_space() and init_time() before set_forcing()")

        ne = self.space.ne
        n_steps = self.time.n_steps

        def _coerce(arr, name: str) -> np.ndarray:
            a = np.asarray(arr, dtype=float)
            if a.ndim == 0:
                return np.full((n_steps, ne), float(a))
            if a.ndim == 1:
                if a.size == ne:
                    return np.broadcast_to(a[np.newaxis, :], (n_steps, ne)).copy()
                if a.size == n_steps:
                    return np.broadcast_to(a[:, np.newaxis], (n_steps, ne)).copy()
                if a.size == 1:
                    return np.full((n_steps, ne), a[0])
            if a.ndim == 2 and a.shape == (n_steps, ne):
                return a.copy()
            raise ValueError(f"'{name}': cannot coerce shape {a.shape} to ({n_steps}, {ne})")

        self._forcing = {
            "precip": _coerce(precip, "precip"),
            "pet": _coerce(pet, "pet"),
            "ptt": _coerce(ptt, "ptt"),
            "lat_gw": _coerce(lat_gw, "lat_gw") if lat_gw is not None else np.zeros((n_steps, ne)),
            "lat_sw": _coerce(lat_sw, "lat_sw") if lat_sw is not None else np.zeros((n_steps, ne)),
        }

    def run(
        self,
        n_steps: Optional[int] = None,
        output_file: Optional[str] = None,
        callback: Optional[Callable] = None,
    ) -> None:
        """
        Run the full simulation using the forcing stored by :meth:`set_forcing`.

        Parameters
        ----------
        n_steps : int, optional
            Run only this many steps.  Defaults to the full temporal domain.
            A value exceeding ``self.time.n_steps`` runs to completion with a warning.
        output_file : str, optional
            Path to NetCDF output file.
        callback : callable, optional
            Called after each step with ``(step_index, state_dict)``.
        """
        import warnings

        if self._forcing is None:
            raise RuntimeError("No forcing available. Call set_forcing() before run().")

        assert self.time is not None, "Call init() before run()"
        assert self.space is not None, "Call init() before run()"

        forcing = self._forcing

        total_steps = self.time.n_steps
        if n_steps is None:
            actual_steps = total_steps
        elif n_steps > total_steps:
            warnings.warn(f"n_steps={n_steps} exceeds total steps={total_steps}; running to completion.")
            actual_steps = total_steps
        else:
            actual_steps = n_steps

        if output_file:
            self._writer = GwswexNCWriter(output_file, self, flush_each_step=self._flush_nc)
            self._writer.write_config()

        ne = self.space.ne
        precip = forcing["precip"]
        pet = forcing["pet"]
        ptt = forcing["ptt"]
        lat_gw = forcing.get("lat_gw", np.zeros((actual_steps, ne)))
        lat_sw = forcing.get("lat_sw", np.zeros((actual_steps, ne)))

        # reset per-run history
        self._mass_balance_history = []

        for t in range(actual_steps):
            # Update root mask if using a dynamic growth model
            if self._has_dynamic_growth and actual_steps > 1:
                frac = t / (actual_steps - 1)
                ir = self._recompute_is_root_at_frac(frac)
                self.update_is_root(ir)

            self.set_lateral(gw=lat_gw[t], sw=lat_sw[t])
            self.step(self.time.dt, precip[t], pet[t], ptt[t])

            state = self.get_state()
            mb = self.get_mass_balance()
            self._mass_balance_history.append(mb)
            if self._writer:
                self._writer.write_timestep(
                    t * self.time.dt,
                    state,
                    {"precip": precip[t], "pet": pet[t], "ptt": ptt[t]},
                )
            if callback:
                callback(t, state)

        if self._writer:
            self._writer.close()
            self._writer = None

    def run_step(self, t: int, *, track: bool = True) -> None:
        """Execute a single macro step at index *t* using the stored forcing.

        Intended for manual step-by-step loops::

            for t in model.Time.steps:
                model.run_step(t)

        Must be preceded by :meth:`set_forcing` and :meth:`init`.

        Parameters
        ----------
        t : int
            Zero-based step index.
        track : bool, optional
            If ``True`` (default), append a mass-balance snapshot to
            ``_mass_balance_history`` and write the current state to the
            built-in NetCDF writer (if open) after each step.  Set to
            ``False`` in performance-sensitive loops where per-step
            diagnostics are not needed; the NetCDF file (including its
            configuration metadata) is still created and closed correctly
            by :meth:`init` / :meth:`deinit`.
        """
        if not self._is_initialised:
            raise RuntimeError("Call model.init() before run_step()")
        if self._forcing is None:
            raise RuntimeError("Call set_forcing() before run_step()")
        assert self.time is not None, "Call init() before run_step()"
        if t < 0 or t >= self.time.n_steps:
            raise IndexError(f"Step t={t} is out of range [0, {self.time.n_steps})")

        if self._has_dynamic_growth and self.time.n_steps > 1:
            frac = t / (self.time.n_steps - 1)
            ir = self._recompute_is_root_at_frac(frac)
            self.update_is_root(ir)

        self.set_lateral(
            gw=self._forcing["lat_gw"][t],
            sw=self._forcing["lat_sw"][t],
        )
        self.step(
            self.time.dt,
            self._forcing["precip"][t],
            self._forcing["pet"][t],
            self._forcing["ptt"][t],
        )

        if track:
            # capture mass balance for this step
            self._mass_balance_history.append(self.get_mass_balance())

            # Write timestep to the built-in NetCDF writer, if opened.
            if self._writer is not None:
                self._writer.write_timestep(
                    t * self.time.dt,
                    self.get_state(),
                    {
                        "precip": self._forcing["precip"][t],
                        "pet": self._forcing["pet"][t],
                        "ptt": self._forcing["ptt"][t],
                    },
                )

    # ------------------------------------------------------------------
    # State access
    # ------------------------------------------------------------------
    def get_state(self) -> dict:
        """Retrieve current model state (in user units)."""
        assert _F is not None, "Kernel not loaded; call init() first"
        assert self.space is not None
        ne, nl = self.space.ne, self.space.nl
        inv_L = 1.0 / self._L_scale
        gw = _F.gwswex_wrapper.get_gw(ne) * inv_L
        gwv = _F.gwswex_wrapper.get_gwv(ne) * inv_L
        sw = _F.gwswex_wrapper.get_sw(ne) * inv_L
        uz = _F.gwswex_wrapper.get_uz(nl, ne) * inv_L
        theta = _F.gwswex_wrapper.get_theta(nl, ne)
        return {"GWH": gw, "GWV": gwv, "SW": sw, "UZ": uz, "theta": theta}

    def get_mass_balance(self) -> dict:
        """Retrieve accumulated flux diagnostics for the last completed step.

        All flux and storage terms are returned in SI (m, m/s where applicable).
        Both solvers populate all fields; for the explicit solver, recharge and
        runoff are derived from storage changes; for the implicit solver they are
        accumulated directly from the Richards flux terms.
        """
        assert _F is not None, "Kernel not loaded; call init() first"
        if self.space is None:
            raise RuntimeError("No spatial domain registered")
        ne = self.space.ne
        acc_p, acc_inf, acc_e, acc_t, acc_r, acc_ro, acc_lg, acc_ls, acc_dgw, acc_dsw, acc_duz, nsub = (
            _F.gwswex_wrapper.get_accumulators(ne)
        )
        return {
            "precip": acc_p,
            "infiltration": acc_inf,
            "evap": acc_e,
            "transp": acc_t,
            "recharge": acc_r,
            "runoff": acc_ro,
            "lat_gw": acc_lg,  # actual applied GW lateral volume change [L]
            "lat_sw": acc_ls,  # actual applied SW lateral depth change [L]
            "delta_gw": acc_dgw,  # total GW storage change this step [L]
            "delta_sw": acc_dsw,  # total SW storage change this step [L]
            "delta_uz": acc_duz,  # total UZ storage change this step [L]
            "n_substeps": nsub,
        }

    # ------------------------------------------------------------------
    # Checkpointing and restart
    # ------------------------------------------------------------------
    def save_checkpoint(self, filepath: str, t: Optional[int] = None) -> None:
        """Save full model state to a NetCDF checkpoint file for restart.

        The checkpoint captures every kernel state field needed to resume the
        simulation from this point under either solver:

        * State pair (prev/curr) for ``GWH``, ``GWV``, ``SW``, ``UZ``, ``theta``
          — written as ``*_prev`` and ``*_curr`` so the first post-restart step
          produces the correct storage deltas.
        * Explicit-solver persistent fields ``IC``, ``ICratio``, ``F_GA``
          (Green–Ampt cumulative infiltration tracker; meaningless for the
          implicit solver but harmless to round-trip).
        * Implicit-solver matric-head field ``h`` (nl, ne) — preserves the
          Picard warm-start and avoids the spurious hydrostatic fallback that
          ``kernel_set_h`` triggers when ``h_prev`` is all-zero.
        * Solver identity (``solver`` global attribute) so :meth:`load_checkpoint`
          can refuse cross-solver loads, and the timestep index ``t``
          (``timestep`` global attribute) so :meth:`list_checkpoints` can report
          where in the run a checkpoint was taken.

        Parameters
        ----------
        filepath : str
            Destination NetCDF file (will be overwritten).
        t : int, optional
            Macro-step index this checkpoint corresponds to.  Stored as a
            global attribute ``timestep``.  Pass the index of the step that
            *just completed* so a restart re-enters the run loop at ``t + 1``.
        """
        import netCDF4 as nc

        assert self.space is not None, "No spatial domain; call init_space() first"
        assert _F is not None, "Kernel not loaded; call init() first"

        ne, nl = self.space.ne, self.space.nl
        # State on the kernel side (in SI). get_state already returns user units;
        # we want to store SI internally so we round-trip cleanly across user
        # unit choices: divide-out the user scale.
        state_user = self.get_state()
        state_si = {k: v * (self._L_scale if k != "theta" else 1.0) for k, v in state_user.items()}

        ic_arr, icratio_arr, f_ga_arr = _F.gwswex_wrapper.get_ic_state(nl, ne, ne)

        with nc.Dataset(filepath, "w") as ds:
            ds.createDimension("ne", ne)
            ds.createDimension("nl", nl)

            # --- Global metadata for list_checkpoints / restart safety ---
            ds.gwswex_checkpoint_version = 1
            ds.solver = self.solver.solver
            ds.timestep = -1 if t is None else int(t)
            if self.time is not None:
                ds.dt_seconds = float(self.time.dt * self._T_scale)
                ds.n_steps_total = int(self.time.n_steps)
            ds.T_unit = self._T_unit
            ds.L_unit = self._L_unit

            # --- Current state (post-step) ---
            for k, v in state_si.items():
                dims = ("ne",) if v.ndim == 1 else ("nl", "ne")
                ds.createVariable(k, "f8", dims)[:] = v

            # --- Explicit-solver persistent fields (always stored; ignored on
            #     load by the implicit solver). ---
            ds.createVariable("IC", "f8", ("nl", "ne"))[:] = ic_arr
            ds.createVariable("ICratio", "f8", ("nl", "ne"))[:] = icratio_arr
            ds.createVariable("F_GA", "f8", ("ne",))[:] = f_ga_arr

            # --- Implicit-solver head profile (only when this run is implicit). ---
            if self.solver.solver == "implicit":
                h_arr = _F.gwswex_wrapper.get_h(nl, ne)
                ds.createVariable("h", "f8", ("nl", "ne"))[:] = h_arr

    def load_checkpoint(self, filepath: str) -> int:
        """Restore model state from a checkpoint file written by
        :meth:`save_checkpoint`.

        After this call the kernel state matches the checkpoint exactly and
        the next :meth:`step` / :meth:`run_step` will advance the simulation
        from the saved point.  Forcing arrays held in ``self._forcing`` are
        unchanged; pass new ``precip``/``pet``/``ptt``/lateral arrays via
        :meth:`set_forcing` (whole-run replacement) or :meth:`update_forcing`
        (per-step replacement) to alter the forcings of the resumed run.

        The returned integer is the timestep index recorded in the checkpoint
        (``-1`` if not recorded).  Resume the run loop at ``returned_t + 1``.

        Cross-solver restart (e.g. saving from explicit and loading into an
        implicit run) is rejected: the implicit solver requires ``h``, which
        the explicit solver does not store, and the explicit solver's IC /
        ICratio / F_GA tracker would be undefined coming from an implicit
        run.

        Parameters
        ----------
        filepath : str
            Path to a NetCDF checkpoint produced by :meth:`save_checkpoint`.

        Returns
        -------
        int
            Timestep index of the checkpoint (``-1`` if absent).
        """
        import netCDF4 as nc

        assert self.space is not None, "No spatial domain; call init_space() first"
        assert _F is not None, "Kernel not loaded; call init() first"

        ne, nl = self.space.ne, self.space.nl
        with nc.Dataset(filepath, "r") as ds:
            ckpt_solver = getattr(ds, "solver", None)
            if ckpt_solver is not None and ckpt_solver != self.solver.solver:
                raise RuntimeError(
                    f"Checkpoint was saved with solver='{ckpt_solver}' but the "
                    f"active model uses solver='{self.solver.solver}'. "
                    "Reconfigure the model with the matching solver before loading."
                )

            gw = np.array(ds["GWH"][:])
            sw = np.array(ds["SW"][:])
            uz = np.array(ds["UZ"][:])
            _F.gwswex_wrapper.set_ic(gw, sw, uz)

            ic_arr = np.array(ds["IC"][:])
            icratio_arr = np.array(ds["ICratio"][:])
            f_ga_arr = np.array(ds["F_GA"][:])
            _F.gwswex_wrapper.set_ic_state(ic_arr, icratio_arr, f_ga_arr)

            if self.solver.solver == "implicit":
                if "h" not in ds.variables:
                    raise RuntimeError(
                        "Implicit-solver restart requires the 'h' field but the "
                        "checkpoint does not contain it (likely written by an "
                        "explicit-solver run)."
                    )
                h_arr = np.array(ds["h"][:])
                _F.gwswex_wrapper.set_h(h_arr)

            return int(getattr(ds, "timestep", -1))

    @staticmethod
    def list_checkpoints(directory: str, pattern: str = "*.nc") -> list[dict]:
        """List GWSWEX checkpoint files in ``directory`` with summary metadata.

        Scans ``directory`` for files matching ``pattern`` (default ``*.nc``)
        and returns one dict per file containing the metadata embedded by
        :meth:`save_checkpoint`.  Files that are not GWSWEX checkpoints
        (missing the ``gwswex_checkpoint_version`` global attribute) are
        skipped silently.

        Returned dicts contain at minimum::

            {
                "path":     <absolute path>,
                "filename": <basename>,
                "timestep": <int, -1 if absent>,
                "solver":   <"explicit" | "implicit" | None>,
                "ne":       <int>,
                "nl":       <int>,
                "dt_seconds":   <float, optional>,
                "n_steps_total":<int,   optional>,
                "T_unit":   <str, optional>,
                "L_unit":   <str, optional>,
            }

        The list is sorted by timestep ascending, then by filename, so that
        printing it directly gives a chronological run timeline.

        Parameters
        ----------
        directory : str
            Directory to scan.
        pattern : str
            Glob pattern within ``directory`` (default ``"*.nc"``).
        """
        import glob
        import os

        import netCDF4 as nc

        results: list[dict] = []
        for path in glob.glob(os.path.join(directory, pattern)):
            try:
                with nc.Dataset(path, "r") as ds:
                    if not hasattr(ds, "gwswex_checkpoint_version"):
                        continue
                    info: dict = {
                        "path": os.path.abspath(path),
                        "filename": os.path.basename(path),
                        "timestep": int(getattr(ds, "timestep", -1)),
                        "solver": getattr(ds, "solver", None),
                        "ne": ds.dimensions["ne"].size if "ne" in ds.dimensions else None,
                        "nl": ds.dimensions["nl"].size if "nl" in ds.dimensions else None,
                    }
                    for opt in ("dt_seconds", "n_steps_total", "T_unit", "L_unit"):
                        if hasattr(ds, opt):
                            info[opt] = getattr(ds, opt)
                    results.append(info)
            except (OSError, KeyError):
                # Not a readable NetCDF, or missing the dimensions we need.
                continue

        results.sort(key=lambda r: (r["timestep"], r["filename"]))
        return results

    # ------------------------------------------------------------------
    # Root mask update (for time-varying vegetation)
    # ------------------------------------------------------------------
    def update_is_root(self, is_root: np.ndarray) -> None:
        """Update the per-layer rooting mask (calls Fortran only when values change)."""
        assert self.space is not None, "No spatial domain"
        assert _F is not None, "Kernel not loaded"
        new_mask = np.asarray(is_root, dtype=np.int32)
        if self.space.is_root is None or not np.array_equal(new_mask, self.space.is_root):
            _F.gwswex_wrapper.set_is_root(np.asfortranarray(new_mask))
            # Mirror the new mask on the Python-side spatial domain so that
            # subsequent equality checks (and any external readers) see the
            # current state.
            object.__setattr__(self.space, "is_root", new_mask)

    # ------------------------------------------------------------------
    # Step-loop helper methods (for use with run_step)
    # ------------------------------------------------------------------
    def update_root_mask(self, is_root: np.ndarray) -> None:
        """Alias for :meth:`update_is_root`, for use inside step-by-step loops."""
        self.update_is_root(is_root)

    def update_lateral_fluxes(self, gw: np.ndarray, sw: np.ndarray) -> None:
        """Pre-set lateral fluxes for the next :meth:`run_step` call.

        This is an alternative to calling :meth:`set_lateral` directly;
        provided for API symmetry with the other ``update_*`` methods.
        """
        self.set_lateral(gw=gw, sw=sw)

    def update_forcing(self, t: int, **kwargs) -> None:
        """Update stored forcing values at step index *t*.

        Only keys present in the stored forcing dict are updated.  Useful for
        injecting time-varying observations into an otherwise pre-loaded forcing
        array during a manual step loop.
        """
        if self._forcing is None:
            raise RuntimeError("Call set_forcing() before update_forcing()")
        assert self.space is not None, "No spatial domain"
        ne = self.space.ne
        for key, val in kwargs.items():
            if key in self._forcing:
                self._forcing[key][t] = (
                    np.broadcast_to(np.asarray(val, dtype=float).ravel()[-1:], (ne,)).copy()
                    if np.asarray(val).size == 1
                    else np.asarray(val, dtype=float).ravel()[:ne]
                )

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------
    @property
    def Time(self) -> Optional[TemporalDomain]:
        """Alias for :attr:`time`, for use in ``for t in model.Time.steps:`` loops."""
        return self.time

    @property
    def mass_balance_history(self) -> list[dict]:
        """Per-step mass balance diagnostics accumulated during :meth:`run` or :meth:`run_step`.

        Each entry is a dict with keys: precip, infiltration, evap, transp,
        recharge, runoff, lat_gw, lat_sw, delta_gw, delta_sw, delta_uz, n_substeps.
        All values are numpy arrays of shape ``(ne,)`` in SI units.

        The list is reset at the start of each :meth:`run` call.
        Returns an empty list if no steps have been executed.
        """
        return self._mass_balance_history
