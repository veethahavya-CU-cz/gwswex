# core.py -- Pydantic configuration and validation models.
from __future__ import annotations

import warnings
from typing import Any, Optional

import numpy as np
import psutil
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


def _cpu_counts() -> tuple[int, int]:
    """Return (physical_cores, logical_cores) for the current machine via psutil."""
    physical = psutil.cpu_count(logical=False) or 1
    logical = psutil.cpu_count(logical=True) or physical
    return (physical, logical)


# ---------------------------------------------------------------------------
# Broadcast utilities (used internally and by model.py)
# ---------------------------------------------------------------------------


def _broadcast_to_ne(val: Any, ne: int, name: str) -> np.ndarray:
    """Coerce *val* to a 1-D float64 array of length *ne*.

    Accepts scalars, length-1 sequences (broadcast to all elements), and
    length-*ne* sequences/arrays.  Raises ``ValueError`` for anything else.
    """
    arr = np.asarray(val, dtype=float).squeeze()
    if arr.ndim == 0:  # scalar
        return np.full(ne, float(arr))
    arr = arr.ravel()
    if arr.size == 1:
        return np.full(ne, arr[0])
    if arr.size == ne:
        return arr.copy()
    raise ValueError(f"'{name}': cannot broadcast size-{arr.size} array to ({ne},)")


def _broadcast_to_nl_ne(val: Any, nl: int, ne: int, name: str) -> np.ndarray:
    """Coerce *val* to a float64 array of shape *(nl, ne)*.

    Accepts:
    * scalar / length-1 array  → broadcast to all layers and elements
    * 1-D length-*nl* array    → same layer values for all elements
    * 2-D shape ``(ne, nl)``   → transposed (user-natural row=element layout)
    * 2-D shape ``(nl, ne)``   → used as-is
    """
    arr = np.asarray(val, dtype=float)
    arr = arr.squeeze()  # collapse size-1 leading/trailing dims
    if arr.ndim == 0:
        return np.full((nl, ne), float(arr))
    if arr.ndim == 1:
        if arr.size == 1:
            return np.full((nl, ne), arr[0])
        if arr.size == nl:
            return np.broadcast_to(arr[:, np.newaxis], (nl, ne)).copy()
        if arr.size == ne:
            return np.broadcast_to(arr[np.newaxis, :], (nl, ne)).copy()
        raise ValueError(f"'{name}': 1-D size {arr.size} cannot broadcast to ({nl}, {ne})")
    if arr.ndim == 2:
        if arr.shape == (nl, ne):
            return arr.copy()
        if arr.shape == (ne, nl):
            return arr.T.copy()
    raise ValueError(f"'{name}': shape {arr.shape} cannot broadcast to ({nl}, {ne})")


class Freezable(BaseModel):
    """Base model that can be frozen to prevent further mutation."""

    model_config = {"arbitrary_types_allowed": True}
    _is_frozen: bool = False

    def freeze(self) -> None:
        object.__setattr__(self, "_is_frozen", True)

    def unfreeze(self) -> None:
        object.__setattr__(self, "_is_frozen", False)

    def __setattr__(self, name, value):
        if self._is_frozen and name != "_is_frozen":
            raise AttributeError(f"Cannot modify frozen {type(self).__name__}.{name}")
        super().__setattr__(name, value)


class VanGenuchtenParams(Freezable):
    """Van Genuchten retention curve parameters for a single material."""

    alpha: float = Field(gt=0, description="VG alpha [1/m]")
    n: float = Field(gt=1, description="VG n [-]")
    theta_r: float = Field(ge=0, description="Residual VWC [-]")
    theta_s: float = Field(gt=0, description="Saturated VWC (porosity) [-]")

    @property
    def m(self) -> float:
        return 1.0 - 1.0 / self.n

    @property
    def Sy(self) -> float:
        return self.theta_s - self.theta_r

    @field_validator("theta_s")
    @classmethod
    def _theta_s_gt_theta_r(cls, v, info):
        if "theta_r" in info.data and v <= info.data["theta_r"]:
            raise ValueError("theta_s must exceed theta_r")
        return v


class Material(Freezable):
    """Soil material definition."""

    id: int = Field(ge=1, description="Material ID (1-based)")
    name: str = ""
    K_sat: float = Field(gt=0, description="Saturated hydraulic conductivity [LT-1]")
    vanG: VanGenuchtenParams
    lam: float = Field(default=0.5, description="Mualem pore-connectivity parameter [-]")


class ETStressParams(Freezable):
    """Laio (2001) piecewise-linear ET stress-function parameters for a vegetation type.

    The stress function s(theta) is fully defined by four saturation thresholds
    (ordered s_h <= s_w <= s_e <= s_star <= 1):
      - below s_h: no ET (hygroscopic limit)
      - s_h to s_w: evaporation only, linearly increasing
      - s_w to s_star: transpiration linearly increasing to full rate
      - s_star to 1: transpiration at potential rate (no stress)
      - s_e controls the lower threshold for capillary-connectivity evaporation
    """

    s_star: float = Field(
        default=0.5, ge=0, le=1, description="Incipient stomatal closure; transpiration at full rate above this [-]"
    )
    s_w: float = Field(default=0.1, ge=0, le=1, description="Wilting point; transpiration zero below this [-]")
    s_h: float = Field(default=0.05, ge=0, le=1, description="Hygroscopic point; evaporation zero below this [-]")
    s_e: float = Field(default=0.5, ge=0, le=1, description="Capillary-continuity threshold for evaporation [-]")


class RootParams(Freezable):
    """Root structure parameters for a vegetation type.

    Defines the maximum rooting depth of a vegetation type. The kernel uses
    this at init to compute `is_root(nl, ne)` from the element's layer
    geometry. Transpiration demand is partitioned uniformly across the rooted
    layers at solve time (1/n_root weighting); no per-layer density profile
    is stored.
    """

    depth: float = Field(gt=0, description="Maximum rooting depth below surface [L]")


class RootGrowthModel(Freezable):  # stub
    """Root growth model parameters.

    Currently only the 'static' model is implemented (root structure fixed at
    the values in RootParams for the full simulation). Logistic and seasonal
    models are reserved for future implementation.
    """

    model: str = Field(
        default="static", description="Growth model: 'static' (only implemented) | 'linear' | 'exponential' (stubs)"
    )
    # --- linear growth stub ---
    # growth_rate: float = Field(default=0.01, gt=0, description="Intrinsic growth rate [1/T] (linear stub)")
    # max_depth: float = Field(default=2.0, gt=0, description="Asymptotic maximum rooting depth [L] (linear stub)")
    # --- seasonal stub ---
    # amplitude: float = Field(default=0.0, ge=0, description="Seasonal depth amplitude [L] (seasonal stub)")
    # phase_offset: float = Field(default=0.0, description="Phase offset from Julian day 1 [T] (seasonal stub)")


class Vegetation(Freezable):
    """Vegetation type definition.

    A vegetation type carries a maximum rooting depth (either a static value
    via `root` or a pair of dynamic-growth endpoints via
    `root_depth_initial` / `root_depth_final`) together with the four Laio (2001)
    ET stress thresholds. Transpiration demand is distributed uniformly across
    the rooted layers of each element at solve time (1/n_root weighting); no
    per-layer density profile is exposed.
    """

    id: int = Field(ge=1, description="Vegetation type ID (1-based)")
    name: str = ""
    et_stress: ETStressParams = Field(default_factory=ETStressParams)
    root: Optional[RootParams] = Field(
        default=None,
        description="Static root geometry (single maximum rooting depth). "
        "Required for static vegetation; ignored when root_depth_initial / "
        "root_depth_final are supplied.",
    )
    root_growth_model: str = Field(
        default="static",
        description="Root growth model: 'static' | 'linear' | 'exponential' (stub)",
    )
    root_depth_initial: Optional[float] = Field(
        default=None,
        gt=0,
        description="Root depth at t=0 as distance below the surface [L]. "
        "Layers with midpoint depth <= root_depth_initial are rooted at t=0. "
        "Used with root_growth_model='linear' or 'exponential'.",
    )
    root_depth_final: Optional[float] = Field(
        default=None,
        gt=0,
        description="Root depth at end of simulation as distance below the surface [L]. "
        "Layers with midpoint depth <= root_depth_final are rooted at t=end. "
        "Used with root_growth_model='linear' or 'exponential'.",
    )

    @field_validator("root_growth_model")
    @classmethod
    def _valid_growth_model(cls, v):
        if v not in ("static", "linear", "exponential"):
            raise ValueError("root_growth_model must be 'static', 'linear', or 'exponential'")
        return v


class SpatialDomain(Freezable):
    """Spatial discretisation and material/root assignment.

    Layer boundaries may be supplied either as a pre-assembled *bnds* array
    ``(ne, nl+1)`` or via separate *top* and *bot* arrays (more natural when
    building the domain interactively).

    *sID* and *uz* ICs may be provided in ``(ne, nl)`` orientation (row = element,
    column = layer) — the natural Python convention — and are transposed internally
    to the ``(nl, ne)`` layout required by the Fortran kernel.
    """

    # Pydantic model config — extra fields (top, bot) accepted and used during
    # the before-validator, then silently discarded.
    model_config = ConfigDict(arbitrary_types_allowed=True, extra="ignore")

    ne: int = Field(ge=1, description="Number of elements")
    nl: int = Field(ge=1, description="Number of layers per element")
    bnds: Optional[np.ndarray] = Field(
        default=None,
        description="Layer boundary elevations, shape (ne, nl+1), strictly decreasing per element. "
        "Provide either this *or* top+bot.",
    )
    sID: Any = Field(
        description="Soil material ID per layer-element. "
        "Accepted shapes: (ne, nl) or (nl, ne) or 1-D length nl/ne; "
        "stored internally as (nl, ne) int32."
    )
    vID: Any = Field(
        description="Vegetation type ID per element. " "Accepted shapes: (ne,), (ne,1), scalar; stored as (ne,) int32."
    )
    is_root: Optional[np.ndarray] = Field(
        default=None,
        description="Root mask per layer-element, shape (nl, ne); derived during init()",
    )
    _materials: dict[int, Material] = {}

    def add_material(self, mat: Material) -> None:
        if self._is_frozen:
            raise AttributeError("Cannot add material to frozen SpatialDomain")
        self._materials[mat.id] = mat

    @model_validator(mode="before")
    @classmethod
    def _preprocess(cls, data):
        """Normalise flexible inputs before field assignment."""
        if not isinstance(data, dict):
            return data

        ne_raw = data.get("ne")
        nl_raw = data.get("nl")
        if ne_raw is None or nl_raw is None:
            return data  # let required-field validators report missing ne/nl

        ne, nl = int(ne_raw), int(nl_raw)

        # ---- bnds from top + bot -----------------------------------------------
        if data.get("bnds") is None:
            top = data.get("top")
            bot = data.get("bot")
            if top is not None and bot is not None:
                top_arr = _broadcast_to_ne(top, ne, "top")
                bot_raw = np.asarray(bot, dtype=float).squeeze()
                if bot_raw.ndim == 0:
                    bot_arr = np.full((ne, nl), float(bot_raw))
                elif bot_raw.ndim == 1:
                    if bot_raw.size == nl:
                        bot_arr = np.broadcast_to(bot_raw[np.newaxis, :], (ne, nl)).copy()
                    elif bot_raw.size == ne:
                        bot_arr = np.broadcast_to(bot_raw[:, np.newaxis], (ne, nl)).copy()
                    else:
                        raise ValueError(f"'bot': 1-D size {bot_raw.size} must equal nl={nl}")
                elif bot_raw.ndim == 2:
                    if bot_raw.shape == (ne, nl):
                        bot_arr = bot_raw.copy()
                    elif bot_raw.shape == (nl, ne):
                        bot_arr = bot_raw.T.copy()
                    else:
                        raise ValueError(f"'bot': shape {bot_raw.shape} doesn't match (ne={ne}, nl={nl})")
                else:
                    raise ValueError(f"'bot': unexpected ndim={bot_raw.ndim}")
                data["bnds"] = np.concatenate([top_arr[:, np.newaxis], bot_arr], axis=1)

        if data.get("bnds") is not None:
            data["bnds"] = np.asarray(data["bnds"], dtype=float)

        # ---- sID → (nl, ne) int32 ----------------------------------------------
        if "sID" in data:
            sid_raw = np.asarray(data["sID"], dtype=np.int32)
            # collapse size-1 leading/trailing dims
            if sid_raw.ndim > 2:
                sid_raw = sid_raw.squeeze()
            if sid_raw.ndim == 0:
                data["sID"] = np.full((nl, ne), int(sid_raw), dtype=np.int32)
            elif sid_raw.ndim == 1:
                if sid_raw.size == nl:
                    data["sID"] = np.broadcast_to(sid_raw[:, np.newaxis], (nl, ne)).copy().astype(np.int32)
                elif sid_raw.size == ne:
                    data["sID"] = np.broadcast_to(sid_raw[np.newaxis, :], (nl, ne)).copy().astype(np.int32)
                else:
                    raise ValueError(f"'sID': 1-D size {sid_raw.size} can't broadcast to ({nl}, {ne})")
            elif sid_raw.ndim == 2:
                if sid_raw.shape == (nl, ne):
                    data["sID"] = sid_raw.astype(np.int32)
                elif sid_raw.shape == (ne, nl):
                    data["sID"] = sid_raw.T.astype(np.int32)
                else:
                    raise ValueError(f"'sID': shape {sid_raw.shape} doesn't match (nl={nl}, ne={ne})")

        # ---- vID → (ne,) int32 -------------------------------------------------
        if "vID" in data:
            vid_raw = np.asarray(data["vID"], dtype=np.int32).squeeze()
            if vid_raw.ndim == 0:
                data["vID"] = np.full(ne, int(vid_raw), dtype=np.int32)
            else:
                vid_flat = vid_raw.ravel()
                if vid_flat.size == 1:
                    data["vID"] = np.full(ne, int(vid_flat[0]), dtype=np.int32)
                elif vid_flat.size == ne:
                    data["vID"] = vid_flat.astype(np.int32)
                else:
                    raise ValueError(f"'vID': size {vid_flat.size} must be 1 or ne={ne}")

        return data

    @model_validator(mode="after")
    def _check_bnds(self):
        if self.bnds is None:
            raise ValueError("Provide either 'bnds' or both 'top' and 'bot'")
        if self.bnds.shape != (self.ne, self.nl + 1):
            raise ValueError(f"bnds must have shape ({self.ne}, {self.nl + 1})")
        for ex in range(self.ne):
            if not all(self.bnds[ex, i] > self.bnds[ex, i + 1] for i in range(self.nl)):
                raise ValueError(f"bnds for element {ex} must be strictly decreasing (surface to bottom)")
        return self

    def _to_fortran(self) -> dict:
        nmat = len(self._materials)
        K_sat = np.array([self._materials[i].K_sat for i in sorted(self._materials)])
        theta_s = np.array([self._materials[i].vanG.theta_s for i in sorted(self._materials)])
        theta_r = np.array([self._materials[i].vanG.theta_r for i in sorted(self._materials)])
        alpha = np.array([self._materials[i].vanG.alpha for i in sorted(self._materials)])
        vg_n = np.array([self._materials[i].vanG.n for i in sorted(self._materials)])
        lam = np.array([self._materials[i].lam for i in sorted(self._materials)])

        # PHASE 6a STUB: Vegetation data marshalling
        nveg = 0
        vID_fortran = self.vID.astype(np.int32)

        return dict(
            ne=self.ne,
            nl=self.nl,
            nmat=nmat,
            nveg=nveg,
            bnds=np.asfortranarray(self.bnds.T),
            sID=np.asfortranarray(self.sID.astype(np.int32)),
            vID=vID_fortran,
            K_sat=K_sat,
            theta_s=theta_s,
            theta_r=theta_r,
            alpha=alpha,
            vg_n=vg_n,
            lam=lam,
            is_root=np.asfortranarray(self.is_root.astype(np.int32)),
        )


class TemporalDomain(Freezable):
    """Time-stepping configuration.

    dt and dt_min are stored in user T units (whatever was passed to GWSWEXmodel).
    n_steps is derived from (stop - start) / dt when start/stop are provided.
    """

    dt: float = Field(gt=0, description="Macro time-step duration [T]")
    n_steps: int = Field(ge=1, description="Number of macro time-steps")
    adaptive: bool = Field(default=True, description="Enable CFL-based adaptive sub-stepping")
    dt_min: float = Field(default=1e-6, description="Minimum sub-step duration [T]")

    @model_validator(mode="before")
    @classmethod
    def _convert_timedeltas(cls, data):
        from datetime import timedelta as _td

        if isinstance(data, dict):
            for key in ("dt", "dt_min"):
                if isinstance(data.get(key), _td):
                    data[key] = data[key].total_seconds()
        return data

    @property
    def steps(self) -> range:
        """Range of step indices for use in manual step loops."""
        return range(self.n_steps)


class ModelParams(Freezable):
    """Green-Ampt infiltration and connectivity parameters.

    ET stress thresholds (s_star, s_w, s_h, s_e) are per-vegetation-type
    and belong to ``ETStressParams`` inside each ``Vegetation`` object.
    """

    psi_f: float = Field(default=0.1, gt=0, description="Green-Ampt suction head [L]")
    F_min: float = Field(default=0.01, gt=0, description="GA minimum cumulative infiltration [L]")
    ICratio_min: float = Field(default=0.05, ge=0, le=1, description="Minimum IC ratio [-]")


class SolverConfig(Freezable):
    """Numerical solver configuration.

    method='explicit': operator-split cascade with CFL-adaptive sub-stepping.
    method='implicit': mixed-form Richards equation with Picard iteration and TDMA.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    solver: str = Field(default="implicit", description="Solver: 'explicit' | 'implicit'")
    omp_threads: int = Field(default=1, ge=1, description="OpenMP thread count for the element-parallel loop.")
    # Explicit solver parameters
    courant_number: float = Field(default=0.9, gt=0, le=1)
    n_trapz: int = Field(default=20, ge=4, description="Quadrature points for UZ_eq")
    beta_hyst: float = Field(
        default=1.0,
        gt=0,
        le=1,
        description="Capillary hysteresis damping",
    )
    # Implicit solver parameters
    picard_tol: float = Field(default=1e-6, gt=0, description="Convergence tolerance on max |Δh| [L]")
    picard_max_iter: int = Field(default=100, ge=1, description="Maximum Picard iterations per macro-step")
    h_min: float = Field(default=-1e6, description="Lower bound on matric head [L]; numerical safety clamp")

    @field_validator("omp_threads")
    @classmethod
    def _valid_omp_threads(cls, v: int) -> int:
        physical, logical = _cpu_counts()
        if v > logical:
            raise ValueError(
                f"omp_threads={v} exceeds the number of logical CPU threads available "
                f"on this machine ({logical}). The OpenMP runtime cannot utilise more "
                f"threads than the OS exposes; reduce omp_threads to <= {logical}."
            )
        if v > physical:
            warnings.warn(
                f"omp_threads={v} exceeds the physical core count ({physical}). "
                f"Hyper-threaded / SMT threads share execution resources; setting "
                f"omp_threads above the physical core count is unlikely to improve "
                f"performance and may degrade it for compute-bound kernels.",
                UserWarning,
                stacklevel=4,
            )
        return v

    @field_validator("solver")
    @classmethod
    def _valid_solver(cls, v):
        if v not in ("explicit", "implicit"):
            raise ValueError("solver must be 'explicit' or 'implicit'")
        return v

    @property
    def solver_type_id(self) -> int:
        """Integer solver type for the Fortran kernel (1=explicit, 2=implicit)."""
        return 1 if self.solver == "explicit" else 2


class InitialConditions(Freezable):
    """Initial state specification."""

    gw: np.ndarray = Field(description="Initial GW elevations, shape (ne,)")
    sw: np.ndarray = Field(description="Initial SW depths, shape (ne,)")
    uz: np.ndarray = Field(description="Initial UZ storages, shape (nl, ne). " "Use -999 to initialise at UZ_eq.")


class LateralFluxes(Freezable):
    """Per-step lateral flux rates."""

    gw: np.ndarray = Field(description="Lateral GW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
    sw: np.ndarray = Field(description="Lateral SW flux rates [LT-1], shape (ne,)")
