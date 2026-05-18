# GWSWEX — Groundwater–Vadose–Surface-Water Exchange Model

[![CI](https://github.com/veethahavya-CU-cz/gwswex/actions/workflows/ci.yml/badge.svg)](https://github.com/veethahavya-CU-cz/gwswex/actions/workflows/ci.yml)

A field-scale, column-based vadose zone model with two numerical solvers,
mass-conserving numerics, and a clean Python API.

---

## Overview

GWSWEX simulates vertical water movement through three coupled reservoirs:

- **Groundwater (GW):** unconfined water table, tracked as head elevation and
  drainable volume.
- **Unsaturated zone (UZ):** soil moisture dynamics per layer, driven by gravity
  drainage, capillary redistribution, ET, and infiltration.
- **Surface water (SW):** ponding depth, receiving precipitation excess and
  discharging to runoff.

Two numerical solvers are available behind a unified Python API:

| Solver | Physics | When to use |
|---|---|---|
| `explicit` | Operator-split cascade with CFL-adaptive sub-stepping (gravity drainage, Green-Ampt infiltration, Gauss-Seidel capillary) | Exploratory runs, explicit benchmarking |
| `implicit` | Mixed-form Richards equation (Celia 1990) with Picard iteration and TDMA | Production runs; accurate water table tracking; recommended default |

The implicit solver is ~25× more accurate in water table tracking than the best
explicit variant and typically faster as well (see
[`docs/model-arch.md` Appendix B](docs/model-arch.md) for the benchmarking summary).

---

## Quick Start

```bash
pip install -e . --no-build-isolation
```

Requires Python ≥ 3.11, numpy, pydantic ≥ 2.6, netCDF4, and a Fortran compiler
(gfortran ≥ 10 or ifort). Build uses Meson + meson-python.

```python
from datetime import datetime, timedelta
import numpy as np
from gwswex import GWSWEXmodel

model = GWSWEXmodel(name="demo", T="h", L="m")

# Spatial domain: 1 element, 5 layers
model.init_space(
    ne=1, nl=5,
    top=[[1.0]],
    bot=[[0.8, 0.5, 0.0, -1.5, -2.0]],
    sID=[[1, 1, 1, 2, 2]],
    vID=[[1]],
)
model.add_material(id=1, name="topsoil", K_sat=1e-5,
                   vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
model.add_material(id=2, name="subsoil", K_sat=5e-6,
                   vanG={"alpha": 1.9, "n": 1.31, "theta_r": 0.065, "theta_s": 0.41})
model.add_vegetation(
    id=1, name="crop",
    et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
    root_depth_initial=0.8,   # rooted layers 1–3 at t=0
    root_depth_final=2.8,     # rooted layers 1–5 by end of run
    root_growth_model="linear",
)

model.init_time(
    start=datetime(2024, 1, 1), stop=datetime(2024, 1, 30),
    dt=timedelta(hours=1), dt_min=timedelta(minutes=1), adaptive=True,
)

model.set_initial_conditions(
    gw=[-1.5], sw=[0.0],
    uz=[[-999, -999, -999, -999, -999]],  # -999 → initialise at UZ equilibrium
)

model.set_solver(solver="implicit", picard_tol=1e-6)

# Forcing: 29-day run with a precipitation burst on hours 100–110
nts = model.time.n_steps
precip = np.zeros((nts, 1))
precip[100:110] = 5e-3   # 5 mm/h
model.set_forcing(precip=precip, pet=1e-4, ptt=2e-4)

model.init()
model.run()

state = model.get_state()
print(f"Final GW head: {state['GWH'][0]:.3f} m")

model.deinit()
```

See [`examples/demo.py`](examples/demo.py) for the full annotated example.

---

## Building from Source

### Environment Setup

The project uses [pixi](https://pixi.sh) for reproducible environment and task management.

**1. Install pixi** (once per machine):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

**2. Install the environment:**

```bash
pixi install
```

Resolves and installs all dependencies — Fortran compiler, Meson, Python, and all runtime and dev packages — into `.pixi/envs/default/`.

**3. Build and install the extension:**

```bash
pixi run build
```

All development commands are available as pixi tasks:

| Task | Description |
|------|-------------|
| `pixi run build` | `meson setup` + compile + editable install |
| `pixi run quick-build` | Recompile + reinstall (no re-setup) |
| `pixi run debug-build` | Build with `-g -O0 -fbacktrace` |
| `pixi run install` | Non-editable release install |
| `pixi run test` | `pytest tests/ -v` |
| `pixi run clean` | Remove `build/`, `build-debug/`, `*.egg-info` |

Use `pixi shell` for an interactive shell with all tools on `PATH`.

---

## User-Facing API

The `GWSWEXmodel` class is the only public entry point.  Everything below is
exposed on the model instance.  All quantities passed in or returned cross
through user units configured at construction (`T="s|min|h|d"`,
`L="m|cm|mm"`); the kernel itself runs in SI internally.

### Construction

```python
GWSWEXmodel(name: str = "gwswex", T: str = "s", L: str = "m")
```

`name` is used as the prefix for output and checkpoint files.  `T` and `L`
fix the user-side time and length units for **every** subsequent argument
(forcings, geometry, ICs, lateral fluxes, root depths, returned state).

### Lifecycle: configure → init → run → deinit

Configuration is built up incrementally through `init_*`, `add_*` and `set_*` calls and committed by a single `init()` call, which performs all cross-component validation, derives the per-element root mask, allocates the Fortran kernel and freezes the Pydantic configuration objects. Failures inside `init()` leave the configuration mutable so you can correct the offending input and retry without rebuilding the whole model.

| Stage                       | Method                                                                | Notes                                                                                              |
| --------------------------- | --------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| 1. Spatial domain           | `init_space(ne, nl, top, bot, sID, vID)` or `init_space(ne, nl, bnds, sID, vID)` | Create grid; `top` and `bot` are convenience inputs; layer boundaries `bnds(ne, nl+1)` are derived |
| 2. Material library         | `add_material(id, name, K_sat, vanG, lam=0.5)`                        | Registered against `sID`; `vanG = {alpha, n, theta_r, theta_s}` (alpha in `1/L`)                   |
| 3. Vegetation library       | `add_vegetation(id, name, et_stress, root..., root_growth_model)`     | Registered against `vID`; see "Vegetation and root growth" below                                   |
| 4. Temporal domain          | `init_time(start, stop, dt, dt_min, adaptive=True, n_steps=None)`     | `dt`, `dt_min` accept `timedelta`; if `n_steps` is omitted, derived from `(stop-start)/dt`         |
| 5. Initial conditions       | `set_initial_conditions(gw, sw, uz)`                                  | `uz=-999` (per layer or scalar) ⇒ kernel initialises that layer at hydrostatic equilibrium         |
| 6. Solver and model params  | `set_solver(**kwargs)` and (optional) `set_model_params(**kwargs)`    | See "Solver configuration" below; must be called before `init()`                                   |
| 7. Forcing                  | `set_forcing(precip, pet, ptt, lat_gw=None, lat_sw=None)`             | Each input broadcasts to `(n_steps, ne)` from scalar / `(ne,)` / `(n_steps,)` / `(n_steps, ne)`    |
| 8. Allocate kernel          | `init()`                                                              | Validates configuration, derives root mask, pushes everything to Fortran and freezes config objects |
| 9. Run                      | `run(...)` or `run_step(t)` in a loop                                 | See "Stepping the model" below                                                                     |
| 10. Release kernel memory   | `deinit()`                                                            | Always call when done; safe to call multiple times                                                 |

### Solver configuration

`set_solver(**kwargs)` accepts any field of `SolverConfig`.  All defaults
live in the Pydantic model; the Fortran kernel never sees an undeclared
value.

| Argument          | Default      | Used by             | Meaning                                                              |
| ----------------- | ------------ | ------------------- | -------------------------------------------------------------------- |
| `solver`          | `"implicit"` | both                | `"explicit"` or `"implicit"`                                         |
| `omp_threads`     | `1`          | both                | OpenMP thread count for the per-element loop                         |
| `courant_number`  | `0.9`        | explicit            | CFL safety factor for adaptive sub-stepping                          |
| `n_trapz`         | `20`         | explicit            | Trapezoidal-rule sub-intervals for VG ePV integral                   |
| `beta_hyst`       | `1.0`        | explicit            | Hysteresis blend factor for $K_{\text{unsat}}$                       |
| `h_min`           | `−10`        | both                | Pressure-head clamp for VG/Mualem evaluation `[L]`                   |
| `picard_tol`      | `1e-6`       | implicit            | Picard convergence tolerance on $\|\Delta h\|_\infty$ `[L]`          |
| `picard_max_iter` | `100`        | implicit            | Maximum Picard iterations per macro-step                             |

`set_model_params(psi_f=0.1, F_min=0.01, ICratio_min=0.05)` controls the
explicit-solver Green–Ampt infiltration and inter-layer connectivity (see
`docs/model-physics.md` for definitions).  The implicit solver ignores all
three.

### Vegetation and root growth

Every element references one vegetation type via `vID`.  Each
`add_vegetation(...)` call registers a type with three logical pieces:

1. **ET stress thresholds** via `et_stress={"s_star", "s_w", "s_h", "s_e"}`
   — Laio (2001) piecewise-linear stress function, applied to both
   evaporation and transpiration limbs.
2. **Root geometry** — which soil layers are active (binary `is_root` mask
   derived from rooting depth at `init()`; transpiration demand is split
   uniformly across rooted layers, i.e. 1/n_root per layer).
3. **Root growth model** — whether the rooted-layer set stays fixed for the
   whole run or expands as a linearly interpolated depth crosses each layer
   midpoint.

**Choosing how to specify roots.**  Pick exactly one of the two input
patterns below; the model raises a clear error if neither is supplied or if
both are supplied simultaneously.

| Pattern                                                   | Use when                                                                                                                    |
| --------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `root={"depth": d}`                                       | Static rooting: layers whose midpoint lies within `d` of the surface are rooted for the entire run; uniform 1/n_root sink   |
| `root_depth_initial=d0, root_depth_final=d1`              | Dynamic rooting: the rooted-layer set expands as the linearly interpolated depth `d(t)` crosses each layer midpoint; use with `root_growth_model="linear"` |

**Root-growth model selector.**  `root_growth_model` accepts
`"static" | "linear" | "exponential"`:

- `"static"` (default): the `root` depth geometry is used unchanged for the
  entire run.
- `"linear"`: with `root_depth_initial`/`root_depth_final` the rooted-layer
  set is re-evaluated at each macro-step at fraction `t / (n_steps − 1)` and
  the Fortran kernel is updated only when the mask changes.
- `"exponential"`: accepted as a Pydantic value today and follows the same
  linear interpolation path as `"linear"`; a true exponential growth law is
  reserved for a future release.

The current rooted-layer mask is accessible as `model.space.is_root`
(shape `(nl, ne)`, dtype int32, 1 = rooted) at any point after `init()`.
To drive the mask from an external source (e.g. a phenology model) call
`update_is_root(is_root)` or its alias `update_root_mask(is_root)`; both
are no-ops when the supplied mask equals the current one.

### Stepping the model

After `init()` two equivalent run modes are available:

- **Single call.** `run(n_steps=None, output_file=None, callback=None)` —
  sweeps `0..n_steps-1` using the forcing stored by `set_forcing(...)`.
  `output_file` opens a CF-1.8 NetCDF writer; `callback(t, state)` fires
  after every macro-step with the post-step `get_state()` dictionary.
- **Manual loop.** `for t in model.Time.steps: model.run_step(t)` —
  identical semantics, lets you interleave per-step diagnostics, on-the-fly
  forcing updates, or checkpoint writes.

Inside (or before) a manual loop you can override forcings without
rebuilding the whole forcing block:

| Method                                          | What it changes                                                                            |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `set_lateral(gw, sw)`                           | Lateral GW and SW flux rates for the **next** `step()`/`run_step()`; consumed once         |
| `update_lateral_fluxes(gw, sw)`                 | Alias for `set_lateral`, provided for symmetry with the other `update_*` methods           |
| `update_forcing(t, **kwargs)`                   | Overwrites stored `precip`/`pet`/`ptt`/`lat_gw`/`lat_sw` at step index `t` in place        |
| `update_is_root(is_root)`                       | Pushes a new `(nl, ne)` binary root mask (int32) to the kernel; no-op if mask unchanged    |
| `step(dt, precip, pet, ptt)`                    | Low-level: advance one macro-step bypassing the stored forcing entirely                    |

### State and diagnostics

| Method                          | Returns                                                                                                                   |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `get_state()`                   | `{ "GWH": (ne,), "GWV": (ne,), "SW": (ne,), "UZ": (nl, ne), "theta": (nl, ne) }` in user units                            |
| `get_mass_balance()`            | Per-step accumulators in SI: `precip, infiltration, evap, transp, recharge, runoff, lat_gw, lat_sw, delta_gw/sw/uz, n_substeps` |
| `mass_balance_history`          | `list[dict]` of every per-step `get_mass_balance()` snapshot accumulated during `run`/`run_step`; reset at each new `run`   |
| `Time` / `time`                 | The frozen `TemporalDomain` (provides `.dt`, `.n_steps`, `.steps` range)                                                  |
| `space`                         | The frozen `SpatialDomain` (provides `.ne`, `.nl`, `.bnds`, `.is_root`)                                                   |
| `solver`                        | The frozen `SolverConfig` (provides `.solver`, `.solver_type_id`)                                                         |

### Checkpointing and restart

The checkpoint API persists every kernel state field needed to resume an
interrupted run with arbitrary new forcings under either solver.

```python
model.save_checkpoint("ckpt_t100.nc", t=100)   # called any time after init()
...
fresh = GWSWEXmodel(...)                        # rebuild with same domain/solver
... configure identically up to init() ...
last_t = fresh.load_checkpoint("ckpt_t100.nc")  # returns the stored t (-1 if absent)
fresh.set_forcing(precip=new_p, pet=new_pe, ptt=new_pt)  # optionally redefine forcing
for t in range(last_t + 1, fresh.time.n_steps):
    fresh.run_step(t)
```

| Method                                                       | Description                                                                                                                                                                                                                                |
| ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `save_checkpoint(filepath, t=None)`                          | Writes a NetCDF file containing GWH/GWV/SW/UZ/theta state, explicit-solver `IC`/`ICratio`/`F_GA`, and (for `solver="implicit"`) the Picard matric-head profile `h(nl, ne)`.  Records solver identity, timestep `t`, dt, and units as global attributes. |
| `load_checkpoint(filepath) -> int`                           | Restores all state into the live kernel.  Refuses to load a checkpoint produced by a different solver (raises `RuntimeError`).  Returns the saved timestep so the caller can resume the loop at `returned_t + 1`.                              |
| `GWSWEXmodel.list_checkpoints(directory, pattern="*.nc")`    | Static method.  Scans a directory for GWSWEX checkpoints, returns one dict per file with `path, filename, timestep, solver, ne, nl, dt_seconds, n_steps_total, T_unit, L_unit`, sorted by timestep.  Non-checkpoint NetCDFs are skipped silently. |

A typical "run with new fluxes from any timestep" workflow is therefore:

```python
# 1. List existing checkpoints to choose a restart point
ckpts = GWSWEXmodel.list_checkpoints("./outputs")
print(ckpts)   # [{'filename': 'ckpt_t100.nc', 'timestep': 100, 'solver': 'implicit', ...}, ...]

# 2. Rebuild the model identically (same domain, materials, vegetation, solver)
m = build_model(solver="implicit")
m.init()
last_t = m.load_checkpoint(ckpts[-1]["path"])

# 3. Override the forcings going forward (e.g. a new precip series)
m.set_forcing(precip=alt_precip, pet=pe, ptt=pt, lat_gw=alt_lat_gw, lat_sw=alt_lat_sw)

# 4. Continue stepping
for t in range(last_t + 1, m.time.n_steps):
    m.run_step(t)
```

For per-step injection of *individual* fluxes (e.g. a one-off lateral
spike) without rewriting the whole forcing block, combine `set_lateral`
and/or `update_forcing(t, ...)` immediately before each `run_step(t)` call.

---

## Repository Structure

```
gwswex/                Python package + Fortran sources
  __init__.py          Public re-exports
  config.py            Pydantic configuration models
  model.py             GWSWEXmodel class (lifecycle, stepping, state access)
  io.py                NetCDF writer/reader
  wrapper.f90          f2py interface
  src/
    kernel.f08         Model singleton, lifecycle, stepping loop
    shared/            constants, types, physics, geometry, lateral,
                       mass_balance
    explicit/          processes, timestep, solver (operator-split cascade)
    implicit/          solver (Richards / Picard / TDMA)
docs/
  model-physics.md     Governing equations, numerical formulation
  model-arch.md        Software architecture
  issues.md            Known limitations and open issues
examples/
  demo.py              Annotated single-column demo
  gwswex-vs-hydrus1d/  HYDRUS-1D comparison notebook
tests/
  test_api.py          Python API unit tests (Pydantic validators, broadcasting)
  test_kernel.py       Fortran kernel integration tests (both solvers)
  test_physics.py      Constitutive-function unit tests
.archive/
  model-diagnosis/     Archived V1–V9 diagnostic workspace
```

---

## Physics Summary

- **Constitutive relations:** Van Genuchten retention curve, Mualem conductivity.
- **Storativity at the water table:** drainable porosity $S_y = \theta_s - \theta_r$
  (not elastic storativity $S_s$, which is appropriate for confined conditions only).
- **ET partitioning:** Laio (2001) piecewise-linear stress function; uniform
  transpiration sink across rooted layers (1/n_root per layer).
- **Implicit solver:** Celia (1990) mixed-form Richards, warm-started between
  macro-steps; adaptive macro-step subdivisions via Picard count proxy.
- **Explicit solver:** gravity drainage → infiltration (Green-Ampt) → capillary
  redistribution (Gauss-Seidel, configurable sweeps) → ET cascade; CFL-adaptive
  sub-stepping.
- **Mass balance:** per-step flux accumulators (precip, infiltration, AE, AT,
  recharge, runoff, lateral GW/SW, ΔGW, ΔSW, ΔUZ) are exposed via
  `model.mass_balance_history`.
- **OMP parallelism:** element-level `!$omp parallel do schedule(dynamic)` in
  the kernel step loop. Thread count set via `set_solver(omp_threads=N)`.

---

## Continuous Integration

The GWSWEX test suite runs automatically on every push to `main` and on all pull
requests via GitHub Actions. Builds test against Python 3.11, 3.12, and 3.13 on
Ubuntu with the latest gfortran.

**Build status:** See the
[Actions tab](https://github.com/veethahavya-CU-cz/gwswex/actions/workflows/ci.yml)
for the latest test results.

**Testing locally:**

```bash
source .env.d/dev.env
python -m pytest tests/ -v
```

All 150 tests must pass before changes are merged to `main`.

---

## Dependencies

| Package | Purpose |
|---|---|
| numpy | Array operations and f2py interface |
| pydantic ≥ 2.6 | Input validation, default management |
| netCDF4 | Output file writing and checkpoint I/O |
| gfortran ≥ 10 | Fortran 2008+ compilation (OMP, coarrays not used) |
| meson + meson-python | Build system |

---

---

## Licence

GWSWEX is released under the **GNU General Public License v3.0** (GPLv3).
This ensures that all improvements, extensions, and derivative works remain
open source and available to the hydrological modelling community.

See [LICENSE](LICENSE) for the full text.

---

## Development and GenAI Attribution

**Conceptualization and physical derivation:** All foundational model concepts, governing equations, and physical process derivations were performed without the aid of generative AI. The development of GWSWEX has been presented incrementally at the European Geosciences Union General Assembly for three consecutive years and forms part of an MSc thesis:

- **2024:** Kootanoor Sheshadrivasan, V. — EGU General Assembly 2024. 
  [https://doi.org/10.5194/egusphere-egu24-8263](https://doi.org/10.5194/egusphere-egu24-8263)

- **2023:** Kootanoor Sheshadrivasan, V. — EGU General Assembly 2023.
  [https://doi.org/10.5194/egusphere-egu23-6432](https://doi.org/10.5194/egusphere-egu23-6432)

- **2022:** Kootanoor Sheshadrivasan, V. — EGU General Assembly 2022.
  [https://doi.org/10.5194/egusphere-egu22-4567](https://doi.org/10.5194/egusphere-egu22-4567)

- **2021:** Kootanoor Sheshadrivasan, V. [*A 0-dimensional Conceptual Model to facilitate Coupling of Groundwater and Surface-Water Numerical Models and its application to a bog-wetland study area.*](https://github.com/Veethahavya/masters-thesis/blob/main/main.pdf) Master of Science thesis, University of Stuttgart.

**Code generation:** Parts of the codebase (primarily software engineering, implementation, and testing infrastructure) were developed with the assistance of generative AI tools, specifically Claude Haiku, Sonnet, and Opus from GitHub Copilot. These tools were used to accelerate development and improve the code quality but did not influence the underlying model physics or scientific methodology.
