"""GWSWEX coupling and BMI integration scaffolding.

This module is a forward-looking placeholder for two interoperability surfaces
that GWSWEX will eventually expose:

1. **Native couplers** to upstream / downstream hydrological models — namely a
   3-D groundwater model (e.g. MODFLOW 6, OpenGeoSys) and a 2-D overland /
   channel surface-water model (e.g. SWMM, ParFlow surface, custom kinematic
   wave). GWSWEX provides the vadose-zone "bridge" between these two domains
   and so will sit in the middle of any coupled stack: it accepts lateral GW
   exchange and lateral SW exchange as forcings (already exposed through
   ``GWSWEXmodel.set_lateral`` / ``set_forcing``) and returns recharge,
   runoff, and water-table state to the surrounding models each macro-step.

2. **BMI** (Basic Model Interface, CSDMS) — a community standard set of
   ``initialize`` / ``update`` / ``finalize`` and getter/setter methods that
   makes a model pluggable into framework-agnostic coupling tooling such as
   BMI-Tester, NextGen, PyMT, and the WRF-Hydro coupler. See the spec at
   https://bmi.csdms.io/en/stable/bmi.spec.html.

Nothing in this module is wired up yet. Every public class and method here is
a stub that ``raise``s ``NotImplementedError`` with a description of what the
production version is expected to do. The stubs are kept here (rather than as
free-standing scribbles) so that:

* downstream users can already import ``gwswex.coupler`` and discover the
  intended coupling API surface;
* the type signatures act as a design contract that future work refines
  in-place rather than reinventing;
* IDE auto-complete and ``help(...)`` already surface the planned interfaces
  to anyone evaluating GWSWEX for use in a coupled framework.

The Fortran-side BMI shim lives at ``gwswex/src/kernel_bmi.f08``; the Python
``BmiGwswex`` class below is intended to wrap that shim once it exists, but
will function as a pure-Python adapter over ``GWSWEXmodel`` until the kernel
shim is implemented.

References
----------
- Hutton, E. W. H., Piper, M. D., & Tucker, G. E. (2020). The Basic Model
  Interface 2.0: A standard interface for coupling numerical models in the
  geosciences. *Journal of Open Source Software*, 5(51), 2317.
- Peckham, S. D., Hutton, E. W. H., & Norris, B. (2013). A component-based
  approach to integrated modeling in the geosciences: The design of CSDMS.
  *Computers & Geosciences*, 53, 3–12.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from .model import GWSWEXmodel


__all__ = [
    "GWCoupler",
    "SWCoupler",
    "CoupledExchange",
    "BmiGwswex",
]


# ─────────────────────────────────────────────────────────────────────────────
# Native coupling surfaces
# ─────────────────────────────────────────────────────────────────────────────


class CoupledExchange:
    """Container for one macro-step of bidirectional flux exchange.

    All quantities are per-element 1-D arrays of length ``ne`` in the model's
    native unit system (``L/T``).

    Attributes
    ----------
    lat_gw_in : np.ndarray
        Lateral groundwater inflow rate **into** the GWSWEX column from the
        external GW model (positive = inflow). Will be passed to
        ``GWSWEXmodel.set_lateral(gw=...)`` for the next macro-step.
    lat_sw_in : np.ndarray
        Lateral surface-water inflow rate **into** the GWSWEX column from the
        external SW model (positive = inflow).
    recharge_out : np.ndarray
        Net internal recharge to the saturated zone produced by GWSWEX during
        the previous macro-step. Read from
        ``GWSWEXmodel.get_mass_balance()['recharge']``. Sign: positive = water
        added to GW.
    runoff_out : np.ndarray
        Net internal runoff produced by GWSWEX during the previous macro-step.
        Read from ``GWSWEXmodel.get_mass_balance()['runoff']``. Sign:
        positive = water added to SW network.
    gwh_out : np.ndarray
        Water-table elevation at end of the previous macro-step
        (``GWSWEXmodel.get_state()['GWH']``). Used by the external GW model
        as an updated head boundary or constraint.
    sw_out : np.ndarray
        Surface-water (ponding) depth at end of the previous macro-step.
    """

    def __init__(
        self,
        ne: int,
        lat_gw_in: NDArray[np.float64] | None = None,
        lat_sw_in: NDArray[np.float64] | None = None,
        recharge_out: NDArray[np.float64] | None = None,
        runoff_out: NDArray[np.float64] | None = None,
        gwh_out: NDArray[np.float64] | None = None,
        sw_out: NDArray[np.float64] | None = None,
    ) -> None:
        z = lambda: np.zeros(ne, dtype=np.float64)  # noqa: E731
        self.ne = ne
        self.lat_gw_in = lat_gw_in if lat_gw_in is not None else z()
        self.lat_sw_in = lat_sw_in if lat_sw_in is not None else z()
        self.recharge_out = recharge_out if recharge_out is not None else z()
        self.runoff_out = runoff_out if runoff_out is not None else z()
        self.gwh_out = gwh_out if gwh_out is not None else z()
        self.sw_out = sw_out if sw_out is not None else z()


class GWCoupler:
    """Adapter to a 3-D groundwater model (e.g. MODFLOW 6, OpenGeoSys).

    The coupler is responsible for two things each macro-step:

    1. Querying the external GW model for the lateral exchange rate to push
       into each GWSWEX column (e.g. by solving the 3-D saturated head field
       given the previous GWSWEX recharge and head boundary).
    2. Pushing GWSWEX's diagnosed recharge volume and updated water-table
       head back into the GW model as a top-boundary forcing for its own
       next solve.

    Implementation patterns expected:

    - **Sequential / loose coupling.** Most common: each model solves its own
      step, exchanges fluxes once, advances. Stable for sub-daily macro-steps
      when fluxes change slowly.
    - **Iterative / strong coupling.** GWSWEX and the GW model iterate within
      one macro-step until lateral flux and head agree to a tolerance.
      Required when the two models can disagree by more than (e.g.) 1% on
      head per step. Will require a "rewind" entry point in the GWSWEX kernel
      (currently absent — see ``todo.md``).
    """

    def __init__(self, model: GWSWEXmodel, gw_model: Any) -> None:
        """Construct a GW coupler.

        Parameters
        ----------
        model : GWSWEXmodel
            The configured (and ``init()``-ed) GWSWEX model instance.
        gw_model : Any
            Handle to the external GW model. Concrete type is framework-
            specific (e.g. a ``flopy`` simulation object for MODFLOW 6, an
            ``ogs`` Project for OpenGeoSys, or a BMI handle).
        """
        self.model = model
        self.gw_model = gw_model
        raise NotImplementedError(
            "GWCoupler is a forward-looking stub; production " "GW coupling is not yet wired. See todo.md."
        )

    def pull_lateral(self) -> NDArray[np.float64]:
        """Query the external GW model for lateral inflow rate per column.

        Returns
        -------
        np.ndarray, shape ``(ne,)``
            Lateral GW flux rate to be set on the GWSWEX side via
            ``model.set_lateral(gw=...)``. Positive = inflow into the column.
        """
        raise NotImplementedError

    def push_recharge(self, recharge: NDArray[np.float64], gwh: NDArray[np.float64]) -> None:
        """Send GWSWEX-diagnosed recharge and head back to the GW model.

        Parameters
        ----------
        recharge : np.ndarray, shape ``(ne,)``
            Per-element net recharge volume from
            ``model.get_mass_balance()['recharge']``.
        gwh : np.ndarray, shape ``(ne,)``
            Updated water-table elevation from ``model.get_state()['GWH']``.
        """
        raise NotImplementedError


class SWCoupler:
    """Adapter to an overland / channel surface-water model.

    Mirror of :class:`GWCoupler` for the surface-water domain. Exchanges:

    - **In:** lateral SW inflow rate per column (e.g. from upstream channel
      cells in a 2-D overland model, or boundary inflow). Pushed via
      ``model.set_lateral(sw=...)``.
    - **Out:** runoff volume produced by GWSWEX during the macro-step
      (``model.get_mass_balance()['runoff']``) and current SW depth
      (``model.get_state()['SW']``), both consumed by the SW model as its
      top-of-column source term.

    Routing topology (which GWSWEX columns drain into which SW reaches) is
    held by the SW-model side; this coupler only pushes per-column volumes.
    """

    def __init__(self, model: GWSWEXmodel, sw_model: Any) -> None:
        self.model = model
        self.sw_model = sw_model
        raise NotImplementedError(
            "SWCoupler is a forward-looking stub; production " "SW coupling is not yet wired. See todo.md."
        )

    def pull_lateral(self) -> NDArray[np.float64]:
        raise NotImplementedError

    def push_runoff(self, runoff: NDArray[np.float64], sw_depth: NDArray[np.float64]) -> None:
        raise NotImplementedError


# ─────────────────────────────────────────────────────────────────────────────
# BMI 2.0 adapter (CSDMS Basic Model Interface)
# ─────────────────────────────────────────────────────────────────────────────


class BmiGwswex:
    """BMI 2.0-compliant adapter for ``GWSWEXmodel``.

    Implements the CSDMS Basic Model Interface 2.0 specification
    (https://bmi.csdms.io/en/stable/bmi.spec.html). The class is a thin
    wrapper over a configured ``GWSWEXmodel``: ``initialize`` parses a YAML /
    TOML config file, calls all ``init_*``/``register_*`` setup, and invokes
    ``model.init()``; ``update`` advances one macro-step via ``run_step``;
    ``finalize`` calls ``model.deinit()``.

    Variable name conventions follow CSDMS Standard Names where possible:

    +-----------------------------------------+------------------------+--------+
    | BMI variable name                       | GWSWEX state           | Units  |
    +=========================================+========================+========+
    | ``soil_water__matric_head``             | ``h`` (implicit only)  | m      |
    | ``soil_water__volume_fraction``         | ``theta``              | -      |
    | ``land_surface_water__depth``           | ``SW``                 | m      |
    | ``soil_water__storage``                 | ``UZ`` summed          | m      |
    | ``land_water_table__elevation``         | ``GWH``                | m      |
    | ``atmosphere_water__precipitation_rate``| forcing ``precip``     | m s-1  |
    | ``land_surface_water__evaporation_rate``| forcing ``pet``        | m s-1  |
    | ``land_vegetation__transpiration_rate`` | forcing ``ptt``        | m s-1  |
    | ``land_water__lateral_groundwater_rate``| forcing ``lat_gw``     | m s-1  |
    | ``land_water__lateral_surface_rate``    | forcing ``lat_sw``     | m s-1  |
    | ``land_water__recharge_rate``           | MB ``recharge`` (out)  | m s-1  |
    | ``land_surface_water__runoff_rate``     | MB ``runoff`` (out)    | m s-1  |
    +-----------------------------------------+------------------------+--------+

    A reference Fortran-side BMI shim lives at
    ``gwswex/src/kernel_bmi.f08`` (also stubbed); this Python class will
    delegate to that shim once the f2py wrapper exposes BMI entry points.
    Until then, the production path will be pure-Python over ``GWSWEXmodel``.
    """

    # CSDMS standard fixed strings --------------------------------------------
    _name = "GWSWEX"
    _input_var_names: tuple[str, ...] = (
        "atmosphere_water__precipitation_rate",
        "land_surface_water__evaporation_rate",
        "land_vegetation__transpiration_rate",
        "land_water__lateral_groundwater_rate",
        "land_water__lateral_surface_rate",
    )
    _output_var_names: tuple[str, ...] = (
        "soil_water__volume_fraction",
        "land_surface_water__depth",
        "land_water_table__elevation",
        "soil_water__storage",
        "land_water__recharge_rate",
        "land_surface_water__runoff_rate",
    )

    def __init__(self) -> None:
        self._model: GWSWEXmodel | None = None
        self._t_index: int = -1

    # ── Control functions ───────────────────────────────────────────────────
    def initialize(self, config_file: str) -> None:
        """Read a config file and bring the model to a runnable state.

        Parameters
        ----------
        config_file : str
            Path to a YAML / TOML file describing geometry, materials,
            vegetation, IC, forcings, and solver. Schema TBD; will likely
            mirror the field names of ``gwswex.config`` Pydantic models.
        """
        raise NotImplementedError

    def update(self) -> None:
        """Advance the model by one macro-step (``dt``)."""
        raise NotImplementedError

    def update_until(self, time: float) -> None:
        """Advance until model time reaches ``time`` (seconds)."""
        raise NotImplementedError

    def finalize(self) -> None:
        """Release resources (calls ``model.deinit()``)."""
        raise NotImplementedError

    # ── Info functions ──────────────────────────────────────────────────────
    def get_component_name(self) -> str:
        return self._name

    def get_input_item_count(self) -> int:
        return len(self._input_var_names)

    def get_output_item_count(self) -> int:
        return len(self._output_var_names)

    def get_input_var_names(self) -> tuple[str, ...]:
        return self._input_var_names

    def get_output_var_names(self) -> tuple[str, ...]:
        return self._output_var_names

    # ── Variable info functions ─────────────────────────────────────────────
    def get_var_grid(self, name: str) -> int:
        """Return the grid id that variable ``name`` lives on.

        GWSWEX uses two grids: id 0 = element-only (1-D length ``ne``);
        id 1 = layered (2-D shape ``(nl, ne)``).
        """
        raise NotImplementedError

    def get_var_type(self, name: str) -> str:
        raise NotImplementedError

    def get_var_units(self, name: str) -> str:
        raise NotImplementedError

    def get_var_itemsize(self, name: str) -> int:
        raise NotImplementedError

    def get_var_nbytes(self, name: str) -> int:
        raise NotImplementedError

    def get_var_location(self, name: str) -> str:
        """Return ``'node'``, ``'edge'``, or ``'face'``. GWSWEX uses ``'node'``."""
        return "node"

    # ── Time functions ──────────────────────────────────────────────────────
    def get_current_time(self) -> float:
        raise NotImplementedError

    def get_start_time(self) -> float:
        raise NotImplementedError

    def get_end_time(self) -> float:
        raise NotImplementedError

    def get_time_units(self) -> str:
        return "s"

    def get_time_step(self) -> float:
        raise NotImplementedError

    # ── Getter / setter functions ───────────────────────────────────────────
    def get_value(self, name: str, dest: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_value_ptr(self, name: str) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_value_at_indices(self, name: str, dest: NDArray[np.float64], inds: NDArray[np.int_]) -> NDArray[np.float64]:
        raise NotImplementedError

    def set_value(self, name: str, src: NDArray[np.float64]) -> None:
        raise NotImplementedError

    def set_value_at_indices(self, name: str, inds: NDArray[np.int_], src: NDArray[np.float64]) -> None:
        raise NotImplementedError

    # ── Grid functions ──────────────────────────────────────────────────────
    def get_grid_rank(self, grid: int) -> int:
        raise NotImplementedError

    def get_grid_size(self, grid: int) -> int:
        raise NotImplementedError

    def get_grid_type(self, grid: int) -> str:
        """GWSWEX grids are ``'rectilinear'`` (uniform per-column layering)."""
        return "rectilinear"

    def get_grid_shape(self, grid: int, shape: NDArray[np.int_]) -> NDArray[np.int_]:
        raise NotImplementedError

    def get_grid_spacing(self, grid: int, spacing: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_grid_origin(self, grid: int, origin: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_grid_x(self, grid: int, x: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_grid_y(self, grid: int, y: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_grid_z(self, grid: int, z: NDArray[np.float64]) -> NDArray[np.float64]:
        raise NotImplementedError

    def get_grid_node_count(self, grid: int) -> int:
        raise NotImplementedError

    def get_grid_edge_count(self, grid: int) -> int:
        raise NotImplementedError

    def get_grid_face_count(self, grid: int) -> int:
        raise NotImplementedError

    def get_grid_edge_nodes(self, grid: int, edge_nodes: NDArray[np.int_]) -> NDArray[np.int_]:
        raise NotImplementedError

    def get_grid_face_edges(self, grid: int, face_edges: NDArray[np.int_]) -> NDArray[np.int_]:
        raise NotImplementedError

    def get_grid_face_nodes(self, grid: int, face_nodes: NDArray[np.int_]) -> NDArray[np.int_]:
        raise NotImplementedError

    def get_grid_nodes_per_face(self, grid: int, nodes_per_face: NDArray[np.int_]) -> NDArray[np.int_]:
        raise NotImplementedError
