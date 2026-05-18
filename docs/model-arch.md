# GWSWEX Model Architecture
This document specifies the GWSWEX implementation architecture. All physics and solver logic are defined in [model-physics.md](model-physics.md); this document maps them to a concrete, layered software architecture without re-deriving the physics.

The package offers two numerical solvers behind a unified API:

- **Explicit** (`solver='explicit'`): operator-split cascade with CFL-adaptive sub-stepping (model-physics.md §3).
- **Implicit** (`solver='implicit'`): mixed-form Richards equation with Picard iteration and TDMA (model-physics.md §4).

The implementation comprises four layers:

1. **Model Kernel** (modern Fortran 2008–2018): the GWSWEX model singleton that stores the structure of the model and dispatches to the selected solver. The explicit path handles CFL evaluation and sub-stepping internally; the implicit path solves the full-column Richards equation per macro-step.
2. **Model Processes** (modern Fortran 2008–2018): all model physics and state management, structured as a set of small, self-contained modules. Explicit-only modules (cascade processes, CFL time-stepper) and implicit-only modules (Richards assembly, TDMA) coexist without coupling.
3. **Fortran f2py wrapper** (F90): a thin interface that exposes kernel data structures, setters, getters, and entry points to Python via f2py, using only primitive-type arguments (scalars, 1-D/2-D arrays of `real(dp)` / `integer`) with no deferred shapes.
4. **Python API** (`gwswex` package): Pydantic-validated end-user library that handles configuration, input pre-processing, macro time-step orchestration, NetCDF I/O, checkpointing / restart, post-processing, and plotting. Designed so that even a non-programming hydrologist can set up and run a simulation from a short, self-documenting script.

---

## 0. Responsibility boundary

| Concern | Lives in |
|---|---|
| Numerical physics (VG, Mualem, Green-Ampt, Laio stress, specific moisture capacity) | Fortran kernel |
| Solver dispatch (explicit vs implicit) | Fortran kernel (`kernel_init`) |
| Explicit: CFL evaluation and adaptive sub-stepping loop | Fortran kernel (explicit path) |
| Explicit: solver cascade (precipitation partitioning through geometry resolution) | Fortran kernel (explicit path) |
| Implicit: Richards matrix assembly, TDMA solve, Picard iteration | Fortran kernel (implicit path) |
| Implicit: water table detection from converged head profile | Fortran kernel (implicit path) |
| State storage (prev/curr pairs) and end-of-step copy | Fortran kernel |
| Explicit: geometry computation ($l^*$, $d_a$, ePV, $\text{UZ}_{\text{eq}}$, $V_{\text{cum}}$) | Fortran kernel (`gwswex_geometry.update_geometry`) |
| Implicit: inter-node harmonic conductivity ($K_{½}$) | Fortran kernel (`gwswex_geometry.compute_K_half`) |
| Mass-balance accumulator bookkeeping within a macro-step | Fortran kernel |
| OpenMP parallelism over elements | Fortran kernel |
| Marshalling derived types to/from flat arrays | F90 wrapper |
| **Parameter default values** (e.g., `courant_number=0.9`, `picard_tol=1e-6`) | Python (Pydantic only) |
| Input validation, unit conversion, config freezing | Python (Pydantic) |
| Macro time-step loop, forcing lookup, recording schedule | Python |
| Lateral flux injection, root-weight updates per macro-step | Python → kernel (setters) |
| NetCDF I/O, checkpointing, restart from snapshot | Python |
| Time-series storage across macro-steps | Python (numpy / NetCDF) |
| Post-processing, analysis, plotting, PET estimation | Python |

The model is **temporally agnostic**: it knows nothing about calendar dates. All durations are in seconds; all lengths in metres. The Python layer handles unit conversion before passing data to the kernel.

---

## 1. Dimension conventions and indexing

| Symbol | Meaning                                                | Scope                 |
| ------ | ------------------------------------------------------ | --------------------- |
| `ne`   | Number of independent model elements (soil columns)    | Global, fixed at init |
| `nl`   | Number of layers per element (uniform across elements) | Global, fixed at init |
| `nmat` | Number of distinct soil materials                      | Global, fixed at init |
| `nveg` | Number of distinct vegetation cover types              | Global, fixed at init |

Layer indexing is **1-based, top-down** (layer 1 = surface, layer `nl` = deepest). Layer boundaries are stored in `bnds(nl+1, ne)`: `bnds(1, ex)` is the surface elevation of element `ex` and `bnds(nl+1, ex)` its domain bottom. The boundary-index-first layout keeps `bnds(:, ex)` contiguous in Fortran column-major order.

All per-layer arrays are `(nl, ne)`, addressed as `arr(lx, ex)` (`lx` = 1-based layer index, `ex` = element index). The `arr(:, ex)` column slice is contiguous — efficient for per-element loops and for assumed-shape `(:)` arguments.

State variables use `*_prev`/`*_curr` pairs (e.g. `UZ_prev(nl,ne)` / `UZ_curr(nl,ne)`) rather than a buffer dimension; at macro-step end, `*_curr` is copied to `*_prev`. The GW state is dual-tracked: `GWH_*` stores head elevation and `GWV_*` stores drainable volume $V_\text{GW}(\text{GW}_t)$, kept in sync throughout the solver; head is recovered via $V_\text{GW}^{-1}$ only when volume changes.

Soil material assignments are `sID(nl, ne)`, indexing into the material library `(1:nmat)`. Each element-layer pair is assigned one soil material type; soil material types are registered on the model itself via `model.add_material(id=..., ...)` and carry the user-supplied Van Genuchten / Mualem parameters (`alpha`, `vg_n`, `lam`, `K_sat`) together with the water-content limits (`theta_s`, `theta_r`). Two further per-material quantities are *derived* inside the kernel and are not user inputs: the Mualem closure exponent `vg_m = 1 − 1/vg_n` and the unconfined specific yield `Sy = theta_s − theta_r` (drainable porosity at the water table; see [model-physics.md §4.1.3](model-physics.md#413-specific-moisture-capacity-ch)). Vegetation assignments are `vID(ne)`, indexing into the vegetation library `(1:nveg)`. Each element is assigned one vegetation type ID; vegetation types are registered on the model itself via `model.add_vegetation(id=..., ...)` and carry the ET stress parameters, rooting depth, and (optionally) a root growth model for that type. The per-layer root mask `is_root(nl,ne)` together with the per-element rooted-layer count `n_root(ne)` is derived from the assigned vegetation type and layer geometry at kernel init, and refreshed between macro-steps when a dynamic growth model is in use; transpiration demand is partitioned uniformly across the rooted layers (1/`n_root` per rooted layer) at solve time.

---

## 2. Data classification

All model data falls into one of six categories. Correct classification drives memory layout, persistence, and API exposure.

### 2.1. State variables (prev/curr pairs, copied at macro-step end)

These hold the previous (`*_prev`) and updated (`*_curr`) model states. At the end of each macro-step the `*_curr` arrays are copied to `*_prev`, unless the step was flagged `dry_run`.

**Shared state (both solvers):**

| Variable | Dimension | Units | Description |
|---|---|---|---|
| `GWH_prev`, `GWH_curr` | `(ne)` | L | GW table elevation per element |
| `GWV_prev`, `GWV_curr` | `(ne)` | L | Drainable GW volume $V_\text{GW}(\text{GW}_t)$ |
| `SW_prev`, `SW_curr` | `(ne)` | L | Surface water (ponding) depth |
| `UZ_prev`, `UZ_curr` | `(nl, ne)` | L | Unsaturated zone storage per layer |
| `theta_prev`, `theta_curr` | `(nl, ne)` | - (VWC) | Derived: volumetric water content. Explicit solver: $\theta = \theta_r + \text{UZ}/d_a$ (clamped to $\theta_s$). Implicit solver: written directly from the converged head profile via van Genuchten. Populated from the hydrostatic profile at `kernel_set_ic` so that `get_state()` is valid immediately after `init()`. |

`theta` is derived but stored for output convenience; it is not independently prognostic. `theta_i`, the blended VWC at the GW boundary layer, is computed on read-back for plotting and is not stored in the kernel.

**Implicit solver state (allocated only when `solver_type == SOLVER_IMPLICIT`):**

| Variable | Dimension | Units | Description |
|---|---|---|---|
| `h_prev`, `h_curr` | `(nl, ne)` | L | Matric head profile per layer. Primary prognostic for the implicit solver; $h < 0$ in unsaturated zone, $h \geq 0$ below water table. Warm-started between steps. |

### 2.2. Persistent parameters (stored, updatable via setters)

Updated between steps (or between macro-steps from Python). They carry memory across time-steps.

**Explicit solver only (not allocated when `solver_type == SOLVER_IMPLICIT`):**

| Variable | Dimension | Description |
|---|---|---|
| `IC` | `(nl, ne)` | Infiltration-front depth tracker per layer |
| `ICratio` | `(nl, ne)` | Inter-layer connectivity ratio (derived from IC) |
| `F_GA` | `(ne)` | Green-Ampt cumulative infiltration depth per element |

### 2.3. Quasi-static parameters (stored, set at init, not updated during run)

| Variable      | Dimension    | Description                                                                                                            |
| ------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------- |
| `bnds`        | `(nl+1, ne)` | Layer boundary elevations (top-down) per element; `bnds(1,ex)` = surface elevation, `bnds(nl+1,ex)` = domain bottom    |
| `vID`         | `(ne)`       | Vegetation type ID per element; `vID(ex)`. Maps each element to one entry in the vegetation library `(1:nveg)`.         |
| `K_sat`       | `(nl, ne)`   | Saturated hydraulic conductivity per element-layer (expanded from per-material input at `kernel_init`)                 |
| `theta_s`     | `(nl, ne)`   | Saturated VWC per element-layer                                                                                        |
| `theta_r`     | `(nl, ne)`   | Residual VWC per element-layer                                                                                         |
| `alpha`       | `(nl, ne)`   | VG alpha parameter per element-layer                                                                                   |
| `vg_n`        | `(nl, ne)`   | VG n parameter per element-layer                                                                                       |
| `vg_m`        | `(nl, ne)`   | Derived: $1 - 1/n$, per element-layer                                                                                  |
| `lambda`      | `(nl, ne)`   | Mualem pore-connectivity parameter per element-layer                                                                   |
| `Sy`          | `(nl, ne)`   | Specific yield $= \theta_s - \theta_r$, per element-layer                                                              |
| `is_root`     | `(nl, ne)`   | Boolean root mask per layer; derived from `Vegetation.root.depth` (or, for time-varying vegetation, from the current interpolated rooting depth) and layer `bnds`; updatable via setter between macro-steps |
| `n_root`      | `(ne)`       | Number of currently rooted layers per element (`sum(is_root(:, ex))`); recomputed automatically whenever `is_root` is (re)assigned and used to weight transpiration demand uniformly (`1 / n_root`) across the rooted layers at solve time |

**Note on material expansion.** The user supplies material properties as `(nmat)` vectors and a soil material ID map `sID(nl, ne)` to `kernel_init`. The kernel immediately expands these to `(nl, ne)` per-layer arrays via the `sID` lookup and stores only the expanded arrays. Neither `sID` nor the original per-material vectors are retained in the model singleton after initialisation. This layout eliminates runtime indirect addressing and enables efficient OpenMP per-element vectorisation.

### 2.4. Ephemeral variables (computed per step, not persisted)

Recalculated at the start of each effective time-step (or sub-step) from the current state.

**Shared:**

| Variable | Dimension | Description |
|---|---|---|
| `gw_bnd_idx` | `(ne)` | GW boundary layer index $l^*$ |
| `d_a` | `(nl, ne)` | Active (unsaturated) thickness |
| `ePV` | `(nl, ne)` | Effective pore volume |
| `V_cum` | `(nl+1, ne)` | Cumulative drainable volume at each layer top (for $V_{\text{GW}}^{-1}$) |

**Explicit solver only:**

| Variable | Dimension | Description |
|---|---|---|
| `UZ_eq` | `(nl, ne)` | Equilibrium UZ storage (VG integral) |
| `K_unsat` | `(nl, ne)` | Unsaturated hydraulic conductivity |
| `tc` | `(nl, ne)` | Transfer capacity $= K_{\text{unsat}} \cdot \delta t$ |

**Implicit solver only:**

| Variable | Dimension | Description |
|---|---|---|
| `K_half` | `(nl-1, ne)` | Upstream-weighted (donor-cell) inter-node conductivity at layer interfaces |
| `sink` | `(nl, ne)` | ET extraction rate per layer [s$^{-1}$] (Laio stress-reduced, distributed by root weight) |

**Implicit solver work arrays** (per-element, reused each Picard iteration):

| Variable | Dimension | Description |
|---|---|---|
| `a`, `b`, `c`, `d_rhs` | `(nl)` | TDMA tridiagonal coefficients and RHS vector |
| `h_k` | `(nl)` | Picard iterate of matric head |
| `theta_k` | `(nl)` | VWC at current iterate |
| `K_k` | `(nl)` | K at current iterate |
| `C_k` | `(nl)` | Specific moisture capacity at current iterate |

### 2.5. Flux accumulators (accumulated across sub-steps, reported per macro-step)

| Variable           | Dimension | Description                        |
| ------------------ | --------- | ---------------------------------- |
| `acc_precip`       | `(ne)`    | Accumulated precipitation volume   |
| `acc_infiltration` | `(ne)`    | Accumulated infiltration volume    |
| `acc_evap`         | `(ne)`    | Accumulated actual evaporation     |
| `acc_transp`       | `(ne)`    | Accumulated actual transpiration   |
| `acc_recharge`     | `(ne)`    | Net internal GW recharge ($\Delta V_\text{GW}^\text{int} - \text{acc\_lat\_gw}$) |
| `acc_runoff`       | `(ne)`    | Net internal SW runoff             |
| `acc_lat_gw`       | `(ne)`    | Accumulated lateral GW flux volume |
| `acc_lat_sw`       | `(ne)`    | Accumulated lateral SW flux volume |
| `acc_delta_gw`     | `(ne)`    | Internal GW volume change $\Delta V_\text{GW}^\text{int}$ (without lateral) |
| `acc_delta_sw`     | `(ne)`    | Internal SW change (without lateral) |
| `acc_delta_uz`     | `(ne)`    | UZ storage change including $\theta_r \cdot \Delta z_\text{sat}$ residual correction |
| `n_substeps`       | `(ne)`    | Number of sub-steps taken          |

### 2.6. Model forcings (set per macro-step from Python)

| Variable      | Dimension | Units | Description                         |
| ------------- | --------- | ----- | ----------------------------------- |
| `precip_rate` | `(ne)`    | LT-1  | Precipitation rate                  |
| `pet_rate`    | `(ne)`    | LT-1  | Potential soil evaporation rate     |
| `ptt_rate`    | `(ne)`    | LT-1  | Potential transpiration rate        |
| `lat_gw_rate` | `(ne)`    | LT-1  | Lateral GW flux rate (+ve = inflow) |
| `lat_sw_rate` | `(ne)`    | LT-1  | Lateral SW flux rate (+ve = inflow) |

---

## 3. Fortran kernel architecture

The kernel is split into small, loosely coupled modules. Each module has a single concern and communicates exclusively through the model singleton. The module dependency graph is strictly acyclic: lower modules never `use` higher ones.

### 3.1. Module dependency graph

```
gwswex_constants          (dp kind, pi, epsilon, status codes)
       |
gwswex_types              (derived types: model singleton, solver config, material library)
       |
gwswex_physics            (pure VG/Mualem/Laio/C(h) functions; no state dependency)
       |
   +---+---+---+---+------+-------------------+
   |       |       |      |                    |
gwswex_  gwswex_   gwswex_ gwswex_          gwswex_
geometry explicit_  time   mass_balance     solver_implicit
   |     processes  |        |             (Richards TDMA,
   |       |        |        |              Picard iteration,
   |       |        |        |              ET sink, WT detection,
   |       |        |        |              compute_K_half_upstream)
   +---+---+---+----+--------+                |
           |              gwswex_lateral       |
           |         (SW & GW lateral flux;    |
           |          used by both solvers)    |
           +---+----------+-------------------+
               |
     gwswex_solver_explicit
     (cascade + sub-stepping)
               |
               +---------------+-------------------+
                               |
                         gwswex_kernel
                         (singleton, solver dispatch, init/deinit)
                               |
                         wrapper.f90
                         (f2py interface)
```

The solver dispatch is resolved at `kernel_init` time: an integer `solver_type` field (1 = explicit, 2 = implicit) selects which solver path `kernel_step` calls. No runtime branching overhead occurs inside the per-element solve loop.

### 3.2. `gwswex_constants`

→ [`gwswex/src/shared/constants.f08`](../gwswex/src/shared/constants.f08)

### 3.3. `gwswex_types`

Defines all kernel derived types; `gwswex_model` is the global singleton holding complete model state.

Key fields of `gwswex_model`:
- `solver_type`: integer (1 = explicit, 2 = implicit), set at `kernel_init`.
- Implicit state arrays: `h_prev(nl, ne)`, `h_curr(nl, ne)` — allocated only when `solver_type == 2`.
- Picard configuration fields: `picard_tol`, `picard_max_iter`. The saturated-zone capacity is the per-material `Sy = theta_s − theta_r` (drainable porosity), derived in `kernel_init`.
- ET stress parameters per vegetation type: `s_star`, `s_w`, `s_h`, `s_e`.

**IMPORTANT:** No default values are defined in Fortran types (see [§8.4](#84-default-parameter-handling)). All parameters are set explicitly by the Python layer at kernel initialisation.

→ [`gwswex/src/shared/types.f08`](../gwswex/src/shared/types.f08)

### 3.4. `gwswex_physics`

Pure, stateless functions for constitutive relationships. Every function is `elemental` or `pure` where possible to enable vectorisation and OpenMP worksharing.

Additional functions for the implicit solver:
- `vg_C(h, ...)`: specific moisture capacity $C(h) = \partial\theta/\partial h$ (model-physics.md §4.1.3).
- `vg_theta_from_h(h, ...)`: water content from matric head (model-physics.md §4.1.1).
- `vg_mualem_kusat(h, ...)`: conductivity from matric head (model-physics.md §4.1.2).

→ [`gwswex/src/shared/physics.f08`](../gwswex/src/shared/physics.f08)

### 3.5. `gwswex_geometry`

Geometric pre-computations structured as two solver-specific subroutines:

**Explicit-solver only:**

- `update_geometry(M, ex)`: translates a GW elevation into boundary-layer index $l^*$, active layer thicknesses $d_a$, effective pore volumes ePV, VG-equilibrium UZ storage, and the cumulative drainable-volume array $V_\text{cum}$. Called at geometry-update points ([model-physics.md §3.4](model-physics.md#34-geometry-update)) and geometry-resolution ([§3.10](model-physics.md#310-geometry-resolution)) in the explicit cascade, and once at initial-condition setup (`kernel_set_ic`). NOT called during implicit solver stepping; the implicit solver derives all geometry from `h_prev` at each Picard iteration.

**Both solvers — exported but not used in the implicit Picard loop:**

- `compute_K_half(K, dz, nl, K_half)`: assembles the inter-node thickness-weighted harmonic-mean conductivity array $K_{½}(1:nl-1)$. This function is exported from `gwswex_geometry` but is **not** called inside `picard_solve`. The implicit Picard loop uses `compute_K_half_upstream` (donor-cell upstream weighting; private to `gwswex_solver_implicit`) instead, which is essential for robust wetting-front advance into dry media (Forsyth & Kropinski, 1997). `compute_K_half` is retained for potential use in diagnostics or alternative solver variants where the harmonic mean is preferred.

→ [`gwswex/src/shared/geometry.f08`](../gwswex/src/shared/geometry.f08)

### 3.6. `gwswex_explicit_processes`

One subroutine per physical process in the explicit solver cascade. Each operates on a single element. The subroutines are called in the order specified by [model-physics.md §3.1](model-physics.md#31-control-flow). All routines mutate the model singleton in-place.

→ [`gwswex/src/explicit/processes.f08`](../gwswex/src/explicit/processes.f08)

### 3.7. `gwswex_time`

CFL evaluation and sub-step loop logic (explicit solver only).

→ [`gwswex/src/explicit/timestep.f08`](../gwswex/src/explicit/timestep.f08)

### 3.8. `gwswex_solver_explicit`

Orchestrates the full explicit solver cascade for a single element within one sub-step, and manages the sub-stepping loop.

→ [`gwswex/src/explicit/solver.f08`](../gwswex/src/explicit/solver.f08)

### 3.9. `gwswex_solver_implicit`

Implements the mixed-form Richards solver with Picard iteration for a single element over one macro-step (model-physics.md §4). Stateless with respect to the module: all inputs and outputs pass through the model singleton or explicit arguments.

Contains the following procedures:

- `solve_element_implicit(M, ex, dt)`: top-level entry point. Evaluates ET sinks, runs Picard iteration, updates GW/SW/UZ from converged head profile.
- `build_picard_system(...)`: assembles the Celia (1990) mixed-form tridiagonal system for one Picard iteration (model-physics.md §4.5).
- `solve_tdma(n, a, b, c, d, x)`: Thomas algorithm, O(n), in-place.
- `compute_et_sink(M, ex, ...)`: evaluates Laio (2001) stress functions per layer and distributes ET as a sink vector (model-physics.md §4.8).
- `locate_water_table(h, z_c, n)`: linear interpolation of the $h = 0$ contour (model-physics.md §4.9).
- `h_to_uz(h, ...)`: converts converged head profile back to UZ storage array for state consistency with the shared output pathway.

Lateral exchange in the implicit solver. The `apply_lateral` shared subroutine is invoked with `apply_gw=.false.`, so only the SW source/sink branch runs at the start of the step. The GW lateral rate `lat_gw_rate(ex)` is instead injected as a distributed Picard source `src_lat_gw(1:nl)` (uniform per saturated thickness; routed to the bottom layer when the column is fully unsaturated) and folded into the per-iterate sink vector inside `picard_solve`, so that the converged head profile is mechanically aware of the lateral inflow rather than being overwritten by `h_to_state` after a direct GWH mutation. The accumulator `acc_lat_gw` is updated separately to `lat_gw_rate * dt`.

Diagnostic accounting details in `solve_element_implicit`:
- `acc_delta_gw` and `acc_delta_uz` are accumulated directly from state updates.
- `acc_runoff` is derived centrally in `gwswex_mass_balance` from $\max(P - I, 0)$; the solver itself only populates `acc_precip` (prescribed) and `acc_infiltration` (recovered from a layer-1 mass balance at the converged Picard iterate).
- **Surface water (ponding):** at the start of the step `SW_curr` is initialised to `SW_prev` and the Neumann top flux is augmented with $\text{SW}_\text{prev}/\Delta t$ so that ponded water can re-infiltrate. After the Picard solve, `SW_curr = \max(P + \text{SW}_\text{prev} - I_\text{actual},\,0)$, matching the explicit-solver convention. `h_to_state` does not touch `SW_curr`; the update is performed by the caller using the layer-1 mass-balance estimate of actual infiltration.
- **Water-table location.** The water table is located as the top of contiguous saturation propagating up from the bottom of the column: the shallowest layer $k$ such that $h_k, h_{k+1}, \dots, h_{n_l} \ge 0$. A positive $h_1$ without saturation of the layer immediately below represents surface ponding (imposed by the Dirichlet top BC, or produced by excess Neumann inflow on a sub-step that has not yet conducted the water downward) and is *not* a rising water table. This distinction is essential: using $h_1 \ge 0$ as a unique indicator of full saturation snaps GWH to the surface elevation the first time the top layer ponds and destroys the UZ/GWV partition.

→ [`gwswex/src/implicit/solver.f08`](../gwswex/src/implicit/solver.f08)

### 3.10. `gwswex_mass_balance`

Per-macro-step output diagnostics, structured in three solver-aware sections:

- `calc_boundary_fluxes(M, ex)` (§A — both solvers): documentation hook confirming that boundary flux accumulators (`acc_evap`, `acc_transp`, `acc_lat_gw`, `acc_lat_sw`, `acc_precip`, `acc_infiltration`) have been populated in-place by the solver. No recomputation; reserved for future post-step flux corrections.
- `calc_storage_deltas(M, ex)` (§B — both solvers): computes $\Delta V_\text{GW}$, $\Delta \text{SW}$, $\Delta \text{UZ}$ from the `*_prev`/`*_curr` state pairs. $\Delta \text{UZ}$ includes an explicit residual-saturated-zone correction $\sum_l \theta_r^{[l]} \cdot [z_{\text{sat,curr}}^{[l]} - z_{\text{sat,prev}}^{[l]}]$, which folds the $\theta_r \cdot z_\text{sat}$ water that both solvers keep outside $\text{UZ}_\text{curr}$ (UZ is zero for fully saturated layers) and outside $\text{GWV}$ (drainable part only, $S_y \cdot z_\text{sat}$) into the storage-change accounting. Without this correction, a vertical migration of the water table silently creates (rise) or destroys (fall) $\theta_r \cdot \Delta z_\text{sat}$ of water in the state-pair difference; with it, the column closure $\text{total} = \text{GWV} + \sum \text{UZ} + \text{SW} + \theta_r \cdot z_\text{sat}$ is exact at every macro-step.
- `calc_discharge_outputs(M, ex)` (§C — explicit solver only): derives net internal GW recharge ($\Delta V_\text{GW} - \text{acc\_lat\_gw}$) and net positive SW change (runoff proxy) from the storage deltas and applied lateral fluxes. Not called for the implicit solver, which accumulates recharge directly in `solve_element_implicit`.
- `output_calc(M, ex, dt_macro)`: entry point called by `kernel_step`; dispatches §A → §B → §C (conditional on solver type).

→ [`gwswex/src/shared/mass_balance.f08`](../gwswex/src/shared/mass_balance.f08)

### 3.11. `gwswex_kernel`

Top-level module that holds the global singleton and provides the macro-step entry point. This is the only module that the f2py wrapper interacts with.

Init is intentionally **two-stage**: `kernel_init` allocates the singleton and loads geometry / material / IC arrays; `kernel_set_vegetation` is then called separately to install the vegetation library and (for vegetation types with a static rooting depth) recompute `is_root` / `n_root`. Splitting the two lets vegetation be reconfigured mid-run (e.g. from a coupled crop model) without reallocating kernel state.

Key entry points:
- `kernel_init(ne, nl, nmat, solver_type, bnds, sID, K_sat, theta_s, theta_r, alpha, vg_n, lam, is_root, ierr)`: allocates all per-element arrays (including `n_root(ne) = sum(is_root(:, ex))`), expands material parameters per layer, derives `vg_m = 1 − 1/vg_n` and `Sy = theta_s − theta_r` per layer, and conditionally allocates the implicit head-state arrays (`h_prev`, `h_curr`) when `solver_type == SOLVER_IMPLICIT`. The `solver_type` integer (1 = explicit, 2 = implicit) is stored on `Model%solver%solver_type` and selects the per-element dispatch in `kernel_step`; there is no separate `set_solver_type` setter.
- `kernel_set_vegetation(nveg, ne, vID, root_depth, s_star, s_w, s_h, s_e, ierr)`: stores the vegetation library and per-element `vID`. For each type with `root_depth > 0` the per-element `is_root` mask is recomputed by tagging every layer whose midpoint lies within `root_depth`; `n_root` is then refreshed accordingly. For types with `root_depth == 0` (the convention used by dynamic-growth vegetation, where the Python layer drives the per-step mask through `set_is_root`) the mask supplied to `kernel_init` is preserved.
- `kernel_set_is_root(rmask)`: replaces the per-element rooting mask with the supplied `(nl, ne)` integer array and refreshes `n_root`. Used between macro-steps when the Python orchestration layer pushes a freshly interpolated root depth for time-varying vegetation.
- `kernel_step(dt, precip, pet, ptt, lat_gw, lat_sw, ierr)`: dispatches to explicit or implicit solver per element:

```
kernel_step(dt, ...)
    !$omp parallel do schedule(static)
    do ex = 1, M%ne
        if (M%solver%solver_type == SOLVER_EXPLICIT) then
            call solve_element_explicit(M, ex, dt)
        else
            call solve_element_implicit(M, ex, dt)
        end if
    end do
```

- `kernel_get_h(h_out)`, `kernel_set_h(h_in)`: head-state getters/setters for the implicit solver.
- `kernel_switch_solver(new_solver, ierr)`: in-place change of the active solver between macro-steps. Translates the live state into the convention required by the destination solver: explicit→implicit warm-starts `h_prev`/`h_curr` from `theta_curr` via `vg_h_inv`, with saturated layers (those whose midpoint lies below `GWH_curr`) assigned a hydrostatic positive head; implicit→explicit resets the persistent Green-Ampt / connectivity fields (`IC = 0`, `ICratio = ICratio_min`, `F_GA = F_min`), refreshes geometry via `update_geometry`, and re-derives `UZ_curr` as `theta_curr * d_a` for active layers (zero for saturated layers). Each layer's `UZ_curr` is then clamped to $(1 - 10^{-6})\cdot\text{ePV}$ and any trim is routed to `SW_curr`, which keeps the column closure exact and prevents the explicit precipitation CFL from collapsing on a column that the implicit solver has left at $\theta = \theta_s$ in one or more layers (see [model-physics.md §3.2.2](model-physics.md#322-precipitation-cfl)). The `h_prev` / `h_curr` arrays are kept allocated even after a switch to the explicit solver so that subsequent switches back to the implicit solver are O(1).
- `kernel_warm_start_explicit(ICratio_in, F_GA_in, ierr)`: optional override of the explicit solver's GA / connectivity persistent state immediately after `kernel_switch_solver(SOLVER_EXPLICIT, ...)`. Both arguments are scalars uniform across layers and elements; a negative value is the sentinel for "leave at the cold-start default". `ICratio_in` is clamped into $[\text{ICratio\_min}, 1]$ and the per-layer wetting-front depth is set to `IC(l, ex) = ICratio_in * d_a(l, ex)` (the inverse of the cascade's own update rule). `F_GA_in` is clamped from below by `F_min`.
- `kernel_warm_start_explicit_proxy(ierr)`: physics-aware seed of the same persistent state, derived from the converged Picard `theta_curr` profile that `kernel_switch_solver` has just translated. Per active layer the effective saturation $S_e = (\theta - \theta_r)/(\theta_s - \theta_r)$ sets `IC = S_e * d_a` and `ICratio = max(S_e, ICratio_min)`; per element the cumulative Green-Ampt infiltration is set to the column integral of $(\theta - \theta_r)\Delta z$ over active layers, capped at $5\,\psi_f$ and floored at `F_min`. See [model-physics.md §3.12](model-physics.md#312-warm-start-of-ga--connectivity-state-on-switch-from-implicit) for the physics rationale and the bias direction.

→ [`gwswex/src/kernel.f08`](../gwswex/src/kernel.f08)

---

## 4. F90 f2py wrapper

The wrapper module is a thin translation layer. It `use`s only `gwswex_kernel` and exposes subroutines with `intent(in/out)` scalar and array arguments that f2py can marshal. No derived types cross the boundary.

The wrapper preserves the kernel's two-stage init: `wrapper.init(...)` calls `kernel_init`; `wrapper.set_vegetation(...)` calls `kernel_set_vegetation`. The Python `GWSWEXmodel.init()` method invokes them in that order (see [`gwswex/model.py`](../gwswex/model.py) `init` and `_call_set_vegetation`).

Key wrapper subroutines (full f2py signatures in [`gwswex/wrapper.f90`](../gwswex/wrapper.f90)):
- `init(ne, nl, nmat, solver_type, bnds, sID, K_sat, theta_s, theta_r, alpha, vg_n, lam, is_root, ierr)`: only the user-supplied VG / Mualem parameters cross the boundary. The kernel derives `vg_m = 1 − 1/vg_n` and `Sy = theta_s − theta_r` internally; these are *never* user inputs. f2py hides the integer dimension arguments (`ne`, `nl`, `nmat`) so the Python signature is dimension-free.
- `set_vegetation(nveg, ne, vID, root_depth, s_star, s_w, s_h, s_e, ierr)`: second stage of init; installs the vegetation library and per-element `vID`. `root_depth(nveg)` carries the per-type rooting depth in SI metres (zero meaning "keep the Python-supplied mask").
- `set_solver_params(courant, dt_min, beta_h, n_trapz, h_min)`: explicit-solver controls plus the shared `h_min` clamp.
- `set_picard_params(picard_tol, picard_max_iter)`: implicit-solver Picard controls.
- `set_model_params(psi_f, F_min, ICratio_min)`: explicit-solver Green–Ampt and connectivity controls.
- `set_omp_threads(n)`, `set_is_root(rmask)`, `set_ic(gw, sw, uz)`, `set_ic_state(IC, ICratio, F_GA)`.
- `switch_solver(new_solver, ierr)`: f2py shim around `kernel_switch_solver`; takes the destination solver-type integer (1 = explicit, 2 = implicit) and returns a status code.
- `warm_start_explicit(icratio_in, f_ga_in, ierr)`: f2py shim around `kernel_warm_start_explicit`; called optionally by Python after a switch to the explicit solver to seed `ICratio` and `F_GA` from user-supplied scalar values.
- `warm_start_explicit_proxy(ierr)`: f2py shim around `kernel_warm_start_explicit_proxy`; called by the Python `switch_solver` default (`warm_start="proxy"`) to seed the GA state from the converged implicit profile.
- `get_h(h_out)`, `set_h(h_in)`: matric head round-trip for implicit-solver checkpoints.
- `get_gw`, `get_gwv`, `get_sw`, `get_uz`, `get_theta`, `get_accumulators`, `get_ic_state`: state and diagnostic accessors.
- `deinit()`: deallocate all kernel arrays.

→ [`gwswex/wrapper.f90`](../gwswex/wrapper.f90)

---

## 5. Python API

### 5.1. Package structure

```
gwswex/
  __init__.py          # public re-exports, __version__
  config.py            # Pydantic configuration models (Freezable base)
  model.py             # GWSWEXmodel class (lifecycle, stepping, state access)
  io.py                # GwswexNCWriter, GwswexNCReader (CF-1.8 NetCDF)
  utils.py             # PETEstimator, AnalyseResults, plot_storage_dynamics
  wrapper.f90          # f2py interface (kept at top of gwswex/)
  src/
    shared/            # constants, types, physics, geometry, lateral,
                       #   mass_balance, kernel
    explicit/          # processes, timestep, solver
    implicit/          # solver
```

### 5.2. Pydantic configuration models (`config.py`)

All configuration types inherit from `Freezable` (Pydantic `BaseModel` subclass with a `.freeze()` method preventing further mutation after registration).

`SolverConfig` covers both solvers. For unconfined conditions the correct water-table storage capacity is `Sy = θ_s − θ_r` (drainable porosity); the compressibility-based specific storativity `Ss` is not applicable and is not included. `Sy` is derived per-material inside the kernel from `theta_s` and `theta_r`; see [model-physics.md §4.1.3](model-physics.md#413-specific-moisture-capacity-ch).

```python
class SolverConfig(Freezable):
    model_config = ConfigDict(extra="forbid")   # rejects unknown kwargs

    solver: Literal["explicit", "implicit"] = "implicit"
    omp_threads: int = Field(default=1, ge=1)
    # Explicit solver parameters
    courant_number: float = Field(default=0.9, gt=0, le=1)
    n_trapz: int = Field(default=20, ge=4)
    beta_hyst: float = Field(default=1.0, gt=0, le=1)
    # Implicit solver parameters
    picard_tol: float = Field(default=1e-6, gt=0)
    picard_max_iter: int = Field(default=100, ge=1)

    @property
    def solver_type_id(self) -> int:
        return 1 if self.solver == "explicit" else 2
```

→ [`gwswex/config.py`](../gwswex/config.py)

### 5.3. Model class (`model.py`)

The `GWSWEXmodel` class manages the full lifecycle: registration of configuration, Fortran kernel init/deinit, macro time-step loop, state retrieval, checkpointing, and restart.

#### 5.3.1. Construction and unit handling

```python
GWSWEXmodel(
    name: str = "gwswex",
    T: str = "s",
    L: str = "m",
    write_output: bool = True,
    output_fpath: str | pathlib.Path | None = None,
)
```

The two unit selectors fix the user-facing time and length units for *every* subsequent argument and return value. The kernel itself is unit-agnostic and runs in SI; the Python layer converts on the way in (multiplying by `_T_scale` / `_L_scale`) and on the way out (dividing).

`write_output` and `output_fpath` control the built-in NetCDF writer. When `write_output=True` (the default), `init()` opens a CF-1.8 `GwswexNCWriter` at `output_fpath`; every `run_step()` writes a timestep automatically, and `deinit()` closes the file. If `output_fpath` is `None`, a default path is derived from `name` in the current working directory using `GWSWEXmodel._default_output_filename(name)` — whitespace is replaced with `_`, the characters ``/ \ : * ? " < > | ' "`` are stripped, consecutive underscores are collapsed, the stem is lower-cased and trimmed to 64 characters, and `.nc` is appended. Passing `write_output=False` disables the built-in writer entirely, in which case the caller is responsible for instantiating `GwswexNCWriter` manually (or using `run(output_file=...)` which opens its own writer for the duration of the sweep).

| Selector | Accepted values        | SI factor stored     |
| -------- | ---------------------- | -------------------- |
| `T`      | `"s" \| "min" \| "h" \| "d"` | seconds-per-user-T   |
| `L`      | `"m" \| "cm" \| "mm"`        | metres-per-user-L    |

#### 5.3.2. Configuration and init pattern

Configuration is built up incrementally through `init_*`, `add_*` and `set_*` calls and committed by a single `init()` call. `init()` performs all cross-component validation (checks `sID`/`vID` references, validates IC dimensions, etc.), derives the per-element root mask from each element's vegetation type and layer geometry, allocates the Fortran kernel, and freezes the Pydantic objects only on full success. Failures inside `init()` leave the configuration mutable so the user can correct the offending input and retry without rebuilding the model.

| Method                                                                  | Stage                | What it does                                                                                                                                                                                                            |
| ----------------------------------------------------------------------- | -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `init_space(ne, nl, top=..., bot=..., bnds=..., sID, vID)`              | spatial              | Construct an unfrozen `SpatialDomain`. Either supply `bnds(ne, nl+1)` directly or the convenience pair `top + bot` (broadcast to per-element columns). `sID` and `vID` accept `(ne, nl)`/`(nl, ne)`/scalar.                |
| `add_material(id, name, K_sat, vanG, lam=0.5)`                          | spatial              | Register a soil material into the in-progress `SpatialDomain.materials`. `vanG = {alpha [1/L], n, theta_r, theta_s}` is a `VanGenuchtenParams` payload; `lam` is the Mualem pore-connectivity exponent.                  |
| `add_vegetation(id, name, et_stress, root..., root_growth_model)`       | vegetation           | Register a vegetation type into `model._vegetation`. See §5.3.3.                                                                                                                                                        |
| `init_time(start, stop, dt, dt_min, adaptive=True, n_steps=None)`       | temporal             | Construct an unfrozen `TemporalDomain`. `dt`/`dt_min` accept `datetime.timedelta` or numeric (in user T units). When `n_steps` is omitted it is derived from `(stop − start)/dt`.                                          |
| `set_initial_conditions(gw, sw, uz)`                                    | ICs                  | Broadcasts `gw, sw` to `(ne,)` and `uz` to `(nl, ne)` (accepts `(ne, nl)` row-major from the user). The sentinel value `−999` for any UZ cell tells the kernel to initialise that cell at hydrostatic equilibrium.            |
| `set_solver(**kwargs)`                                                  | solver               | Build and freeze a `SolverConfig`. Must be called *before* `init()`. To change solver *after* `init()`, use `switch_solver(**kwargs)` (in-place, state-preserving) or `deinit()` + `set_solver(...)` + `init()` + `set_initial_conditions(...)` (full reconfigure).                                                |
| `switch_solver(**kwargs)`                                               | solver (post-init)   | In-place change of the active solver on a live kernel. Builds a new frozen `SolverConfig`, calls the kernel translator (see [§3.11](#311-gwswex_kernel)), re-pushes `set_solver_params` (and `set_picard_params` if the destination is implicit) and preserves the previous `omp_threads`. State (`GWH`, `SW`, `UZ`, `theta`, `h`) is preserved and translated; mass-balance accumulators continue across the boundary. Requires `_is_initialised`. The keyword-only `warm_start` argument selects how the explicit GA / connectivity state is seeded after a switch to explicit: `"proxy"` (default) derives `IC`, `ICratio` and `F_GA` from the converged Picard $\theta$ profile; `"cold"` keeps the kernel-applied cold defaults (`IC = 0`, `ICratio = ICratio_min`, `F_GA = F_min`); `"manual"` accepts the additional scalars `icratio_init` and `f_ga_init` for a uniform user-supplied seed (`IC = icratio_init * d_a` per layer, `F_GA = f_ga_init` per element).
| `set_model_params(psi_f=0.1, F_min=0.01, ICratio_min=0.05)`             | model params         | Explicit-solver Green–Ampt and inter-layer connectivity controls. The implicit solver ignores all three.                                                                                                                  |
| `set_forcing(precip, pet, ptt, lat_gw=None, lat_sw=None)`               | forcing              | Coerces every input to shape `(n_steps, ne)` in user L/T units. Scalars, `(ne,)`, `(n_steps,)`, and `(n_steps, ne)` shapes all broadcast cleanly; anything else raises.                                                  |
| `init()`                                                                | kernel allocation    | Validates `sID`/`vID` references, derives `is_root(nl, ne)` and `n_root(ne)` from each element's vegetation type and layer geometry, freezes all Pydantic objects, then marshals geometry/material/vegetation arrays to SI and calls `wrapper.init`, `set_solver_params`, `set_picard_params` (implicit only), `set_omp_threads`, `set_model_params`, `set_ic`, `set_vegetation`. |
| `deinit()`                                                              | kernel teardown      | Calls `wrapper.deinit()` and closes any open NetCDF writer; safe to call repeatedly.                                                                                                                                     |

`SolverConfig` arguments (all defaults Pydantic-side):

| Field             | Default      | Used by   | Meaning                                                          |
| ----------------- | ------------ | --------- | ---------------------------------------------------------------- |
| `solver`          | `"implicit"` | both      | `"explicit"` or `"implicit"`                                     |
| `omp_threads`     | `1`          | both      | OpenMP thread count for the per-element loop                     |
| `courant_number`  | `0.9`        | explicit  | CFL safety factor for adaptive sub-stepping                      |
| `n_trapz`         | `20`         | explicit  | Trapezoidal-rule sub-intervals for the VG ePV integral           |
| `beta_hyst`       | `1.0`        | explicit  | Hysteresis blend factor for $K_{\text{unsat}}$                   |
| `h_min`           | `−1e6`       | both      | Pressure-head clamp for VG/Mualem evaluation `[L]`. Numerical safety floor only; not a physical limit. The default of $-10^6$ m is intentionally far below any realistic soil suction (residual-water suctions are typically $-10^2$ to $-10^4$ m) so that the clamp activates only on Picard blow-up. Calibrate downward only when a physically tighter bound is required for diagnosis. |
| `picard_tol`      | `1e-6`       | implicit  | Picard convergence tolerance on $\|\Delta h\|_\infty$ `[L]`      |
| `picard_max_iter` | `100`        | implicit  |                                                                  |

#### 5.3.3. Vegetation types and root-growth models

Each vegetation type is a frozen `Vegetation` object with three logical pieces: ET stress thresholds, a root-distribution specification, and a growth model selector.

**ET stress (`et_stress`).** A `ETStressParams` dict `{s_star, s_w, s_h, s_e}` carrying Laio (2001) thresholds in saturation units; thresholds must satisfy `s_h ≤ s_w ≤ s_e ≤ s_star ≤ 1`. Used by both solvers (`compute_et_sink` for implicit, `evaporation`/`transpiration` cascade subroutines for explicit).

**Root distribution.** Pick exactly one of the two mutually exclusive input patterns:

| Pattern                                                                    | Behaviour at `init()`                                                                                                                                                                                                                                                                          |
| -------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `root = RootParams(depth)`                                                 | Static rooting: `is_root(:, ex)` is set for every layer whose midpoint is within `depth` of the surface. The mask is pushed to the kernel once at `init()` and never updated.                                                                                                                  |
| `root_depth_initial = d0`, `root_depth_final = d1` (with `root_growth_model != "static"`) | Dynamic rooting: at each macro-step the current rooting depth is `d0 + (d1 − d0) · t/(n_steps − 1)`; layers whose midpoint depth ≤ current depth are flagged. The Python helper `_recompute_is_root_at_frac(frac)` rebuilds the mask and pushes it via `update_is_root` (a no-op when the mask has not changed since the previous step). |

Transpiration demand is partitioned uniformly over the rooted layers (`1 / n_root` per rooted layer) at solve time — there is no per-layer root-density profile.

**Root-growth model (`root_growth_model`).** Selector with three accepted values:

- `"static"` (default): the rooting mask derived from `RootParams(depth)` is pushed once at `init()` and never updated. The kernel-side `is_root` array is read-only thereafter unless the user explicitly calls `update_is_root(...)`.
- `"linear"`: requires `root_depth_initial`/`root_depth_final`. At each step the Python loop in `run`/`run_step` calls `_recompute_is_root_at_frac(t / (n_steps − 1))` and pushes the new mask to the kernel via `update_is_root`. The kernel-side comparison short-circuits to a no-op when the mask has not changed (e.g. between two layer-midpoint crossings).
- `"exponential"`: accepted as a Pydantic value; current implementation reuses the same linear-interpolation update path. Reserved for a future genuine exponential growth law.

The `_recompute_is_root_at_frac` helper is fully internal; user code that needs externally-driven phenology (e.g. a coupled crop model) should call `update_is_root(mask)` directly.

#### 5.3.4. Stepping the kernel

After `init()` two equivalent stepping modes are available; both consume the forcing block stored by `set_forcing(...)`.

| Method                                                            | Use                                                                                                                                                                                                              |
| ----------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run(n_steps=None, output_file=None, callback=None)`              | Sweeps `0..n_steps-1`. `output_file` opens a CF-1.8 `GwswexNCWriter`; `callback(t, state)` fires after every macro-step with the post-step `get_state()` dict. Resets `mass_balance_history` at entry.            |
| `run_step(t)`                                                     | Executes step `t` exactly as `run()` would (including dynamic root-weight refresh); the step counter must remain in `[0, n_steps)`. Raises if `set_forcing` has not been called.                                  |
| `step(dt, precip, pet, ptt)`                                      | Low-level: bypasses stored forcing entirely. Use to drive the kernel from outside the model's forcing book-keeping (e.g. a coupled framework). `set_lateral(...)` must be called first if non-zero lateral flux is wanted; lateral is consumed once. |

Per-step / pre-step modifiers usable inside a manual `run_step(t)` loop (or before a `step()` call):

| Method                                                | Effect                                                                                                                                                |
| ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `set_lateral(gw, sw)`                                 | Pre-loads lateral GW and SW flux *rates* (user `L/T`) for the next macro-step. Consumed and cleared inside `step()`.                                   |
| `update_lateral_fluxes(gw, sw)`                       | Alias for `set_lateral`.                                                                                                                              |
| `update_forcing(t, **kwargs)`                         | Mutates stored forcing in place at index `t`. Only keys present in `self._forcing` (`precip`, `pet`, `ptt`, `lat_gw`, `lat_sw`) are accepted.            |
| `update_is_root(mask)` / `update_root_mask(mask)`     | Pushes `(nl, ne)` integer rooting mask to kernel; short-circuits when values are unchanged. The kernel automatically refreshes `n_root(ex) = sum(is_root(:, ex))` on every assignment.   |

#### 5.3.5. State and diagnostic accessors

| Method                            | Returns                                                                                                                                                                                                                   |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `get_state()`                     | `dict` of `GWH(ne)`, `GWV(ne)`, `SW(ne)`, `UZ(nl, ne)`, `theta(nl, ne)` — all in user units.                                                                                                                                |
| `get_mass_balance()`              | `dict` of accumulator arrays, all `(ne,)` in SI: `precip`, `infiltration`, `evap`, `transp`, `recharge`, `runoff`, `lat_gw`, `lat_sw`, `delta_gw`, `delta_sw`, `delta_uz`, `n_substeps`. The accumulators are populated identically by both solvers (see §3.10). |
| `mass_balance_history` (property) | `list[dict]` of every per-step `get_mass_balance()` snapshot collected during `run` / `run_step`. Reset at the start of each `run`; manual `run_step` calls append.                                                          |
| `Time` (alias for `time`)         | The frozen `TemporalDomain`; provides `.dt`, `.n_steps`, `.steps` (range).                                                                                                                                                  |
| `space`                           | The frozen `SpatialDomain`; provides `.ne`, `.nl`, `.bnds`, `.is_root`.                                                                                                                                                     |
| `solver`                          | The frozen `SolverConfig`; provides `.solver` (string) and `.solver_type_id` (1 = explicit, 2 = implicit).                                                                                                                  |

#### 5.3.6. Checkpointing and restart

The checkpoint API persists every kernel state field needed to resume a run with arbitrary new forcings under either solver.

| Method                                                                  | Purpose                                                                                                                                                                                                                                                                                |
| ----------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `save_checkpoint(filepath, t=None)`                                     | Writes a NetCDF file with the full kernel state in SI: state pair (`GWH`, `GWV`, `SW`, `UZ`, `theta`); explicit-solver persistent fields (`IC`, `ICratio`, `F_GA`); for `solver="implicit"` also the matric-head profile `h(nl, ne)`. Records `solver`, `timestep`, `dt_seconds`, `n_steps_total`, `T_unit`, `L_unit`, and a `gwswex_checkpoint_version` global attribute. |
| `load_checkpoint(filepath) → int`                                       | Restores all state into the live kernel via `set_ic`, `set_ic_state`, and (for implicit) `set_h`. Refuses cross-solver loads with a `RuntimeError`. Returns the saved `timestep` index (`-1` if absent) so the caller can resume the loop at `returned_t + 1`.                            |
| `GWSWEXmodel.list_checkpoints(directory, pattern="*.nc")` (static)      | Scans `directory`, returns one dict per file with `path`, `filename`, `timestep`, `solver`, `ne`, `nl`, optional `dt_seconds`, `n_steps_total`, `T_unit`, `L_unit`. Sorted by `timestep` ascending. Files lacking the `gwswex_checkpoint_version` attribute are skipped silently, so a run-output directory can be passed without false hits. |

The Picard warm-start in the implicit solver depends on `h_prev` being the converged head profile from the previous macro-step. Without `h` round-tripping, a restart would silently fall back to the all-zero sentinel branch in `kernel_set_h`, which substitutes a hydrostatic profile and distorts the first post-restart step. `save_checkpoint` therefore writes `h` whenever the active solver is implicit, and `load_checkpoint` requires it to be present in that case.

A typical "restart from any timestep with new fluxes" workflow:

```python
ckpts = GWSWEXmodel.list_checkpoints("./outputs")
m = build_model(solver="implicit")     # rebuild with same domain/solver
m.init()
last_t = m.load_checkpoint(ckpts[-1]["path"])
m.set_forcing(precip=alt_p, pet=alt_pe, ptt=alt_pt,
              lat_gw=alt_lat_gw, lat_sw=alt_lat_sw)
for t in range(last_t + 1, m.time.n_steps):
    m.run_step(t)
```

For a one-off per-step override (e.g. injecting an observed lateral spike at one timestep without rewriting the whole forcing block), call `set_lateral(...)` and/or `update_forcing(t, ...)` immediately before each `run_step(t)`.

→ [`gwswex/model.py`](../gwswex/model.py)

### 5.4. NetCDF I/O (`io.py`)

→ [`gwswex/io.py`](../gwswex/io.py)

### 5.5. End-user demo script

This is the canonical user experience: a complete simulation in ~50 lines of self-documenting code. Solver selection is a single keyword argument to `set_solver()`.

```python
from datetime import datetime, timedelta
import numpy as np
from gwswex import GWSWEXmodel

model = GWSWEXmodel(name="demo", T="h", L="m")

# Spatial domain
model.init_space(ne=1, nl=5, top=[[1.0]], bot=[[0.8, 0.5, 0.0, -1.5, -2.0]],
                 sID=[[1, 1, 1, 2, 2]], vID=[[1]])
model.add_material(id=1, name="topsoil", K_sat=1e-5,
                   vanG={"alpha": 3.6, "n": 1.56, "theta_r": 0.078, "theta_s": 0.43})
model.add_material(id=2, name="subsoil", K_sat=5e-6,
                   vanG={"alpha": 1.9, "n": 1.31, "theta_r": 0.065, "theta_s": 0.41})
model.add_vegetation(id=1, name="crop1",
                     et_stress={"s_star": 0.5, "s_w": 0.1, "s_h": 0.05, "s_e": 0.5},
                     root_growth_model="linear",
                     root_depth_initial=0.3,
                     root_depth_final=0.8)

# Temporal domain
model.init_time(start=datetime(2024, 1, 1), stop=datetime(2024, 1, 30),
                dt=timedelta(hours=1), dt_min=timedelta(minutes=1), adaptive=True)

# Initial conditions (uz=-999 → UZ_eq in kernel)
model.set_initial_conditions(gw=[-1.5], sw=[0.0],
                             uz=[[-999, -999, -999, -999, -999]])

# Solver and forcing
model.set_solver(solver="implicit", picard_tol=1e-6)  # or solver="explicit"
precip = np.zeros((model.time.n_steps, 1))
model.set_forcing(precip=precip, pet=1e-4, ptt=2e-4)

model.init()
model.run()   # or: for t in model.Time.steps: model.run_step(t)

state  = model.get_state()
print(f"Final GW head: {state['GWH'][0]:.3f} m")
print(f"Total ET:      {sum(mb['evap'][0] + mb['transp'][0] for mb in model.mass_balance_history):.4f} m")

model.deinit()
```

→ [`examples/demo.py`](../examples/demo.py)

---

## 6. Build system

The project uses **Meson** as the build system, with **meson-python** for PEP-517 wheel building and **f2py** for the Fortran-Python bridge.

### 6.1. Compilation order

```
constants.f08
types.f08
physics.f08
geometry.f08
explicit_processes.f08
time.f08
solver_explicit.f08
solver_implicit.f08
mass_balance.f08
kernel.f08
    → compiled as static library: libgwswex_core.a

wrapper.f90 + libgwswex_core.a
    → f2py generates: f_gwswex.cpython-3XX.so / .pyd
```

The module dependency graph is strictly acyclic. `solver_implicit.f08` depends on `gwswex_physics` and `gwswex_geometry` but has no dependency on any explicit-solver module, and vice versa. Both solver modules are linked into the same static library; unused code paths incur no runtime cost.

### 6.2. `meson.build` (sketch)

→ [`meson.build`](../meson.build)

---

## 7. Testing strategy

The in-tree pytest suite is fully standalone (no HYDRUS / phydrus / external
binaries required).  Coverage is split across three files:

| File | What is tested |
|---|---|
| `tests/test_physics.py` | Van Genuchten retention/conductivity, Mualem $K(h)$, $C(h) = \mathrm{d}\theta/\mathrm{d}h$, TDMA solver against `numpy.linalg.solve`, Laio (2001) ET stress functions |
| `tests/test_api.py` | Pydantic validators (`SpatialDomain`, `TemporalDomain`, `SolverConfig`, `Material`, `VanGenuchtenParams`, `Vegetation`, `RootParams`, `ModelParams`, `InitialConditions`, `LateralFluxes`), broadcast utilities, `Freezable.freeze`, `GWSWEXmodel` lifecycle guards and forcing broadcast |
| `tests/test_kernel.py` | End-to-end Fortran kernel runs (both solvers): zero-forcing smoke tests, multi-element broadcasting, mass-balance accumulator closure, lateral-flux accumulators, ET extraction, runoff/ponding under saturated profile, checkpoint round-trip, NetCDF output (CF-1.8), error paths, `run()` vs `run_step()` consistency |

The HYDRUS-1D comparison lives in
`examples/gwswex-vs-hydrus1d/comparison-basic-loam.ipynb` and is *not* part of the pytest
suite (it is treated as an external benchmark, not a regression test).

### 7.1. Reference implementation (`reference.py`)

A pure-Python reference implementation of the complete solver cascade: unoptimised, written for clarity, with every equation implemented directly from [model-physics.md](model-physics.md) with minimal abstraction. Ground-truth for unit and integration tests; the Fortran kernel must produce identical results to floating-point tolerance for all test configurations.

---

## 8. Design principles and extensibility

### 8.1. Adding a new process module

1. Create `gwswex/new_process.f08` containing module `gwswex_new_process`. `use` only `gwswex_constants`, `gwswex_types`, and optionally `gwswex_physics`.
2. Define a subroutine operating on a single element via the model singleton.
3. Add any new state variables to `gwswex_types::gwswex_model`.
4. Insert the call at the appropriate point in the relevant solver module (`gwswex_solver_explicit` for the explicit path, `gwswex_solver_implicit` for the implicit path, or both if shared).
5. Add the new source to `meson.build` (before the solver module(s) that `use` it).
6. Expose any new state/setters through the f2py wrapper.
7. Add a Pydantic model for any new parameters in `core.py`.
8. No existing module's code changes, except the insertion in the solver cascade.

### 8.2. Replacing the ET stress function

1. Add an alternative function (e.g. `feddes_stress`) in `physics.f08`.
2. Add a selector field to `et_params` (e.g. `stress_model = 'laio' | 'feddes'`).
3. In the explicit solver: dispatch in `gwswex_explicit_processes::evaporation` and `transpiration`.
4. In the implicit solver: dispatch in `gwswex_solver_implicit::compute_et_sink`.
5. No other module changes.

### 8.3. Fortran code style

- Vectorise array operations where loop-free expressions are possible (`where`, `merge`, `min`, `max` on arrays).
- Avoid branching inside hot loops; prefer masked operations or branchless `min`/`max`/`merge`.
- All per-element work arrays are pre-allocated in the singleton; no allocation inside the stepping loop.
- OpenMP parallelism is over elements (outer loop) with `schedule(static)` (the kernel uses static chunking on the per-element dispatch loop in `kernel_step`; the initially considered `dynamic` schedule was rejected after benchmarking because per-element work is dominated by Picard/sub-step counts whose cross-element variance is small relative to the dynamic-scheduling overhead at typical $n_e$). For very heterogeneous element loads (e.g. spatially clustered ponding), revisiting `schedule(dynamic, chunk)` is a worthwhile sensitivity experiment.
- Each element's data is contiguous via column-major `(nl, ne)` layout — `arr(:, ex)` is a contiguous column slice, ensuring cache-friendly access in per-element loops.

### 8.4. Default parameter handling (CRITICAL PRINCIPLE)

**RULE: NO default values shall ever be defined on the Fortran side. ALL defaults must be defined exclusively in the Python Pydantic models.**

**Documented exception.** `solver_config%h_min = -1.0e6_dp` in [`gwswex/src/shared/types.f08`](../gwswex/src/shared/types.f08) carries a Fortran-side default that mirrors the Python default in [`gwswex/config.py`](../gwswex/config.py). This exists as a defensive guard so that any code path that touches `vg_C`/`vg_theta_from_h`/`vg_mualem_kusat` before `set_solver_params` has executed (e.g. tests that exercise the kernel without a full Python init, or future hot-reload paths) cannot dereference an uninitialised pressure-head clamp and produce silent NaNs in the constitutive evaluations. The two values are kept identical by code review; any change to the Python default must be mirrored in the Fortran type and a regression test added. New parameters must not follow this pattern.

Rationale: (1) **Transparency + auditability** — all defaults are visible and traceable in Python `Field(default=...)` declarations, never hidden in compiled Fortran. (2) **Single source of truth** — no dual-definition drift or maintenance burden from defaults specified in two places. (3) **Kernel purity + interoperability** — the Fortran kernel is a stateless executor accepting only explicitly set values; future C++/Julia bindings can swap the Python layer without Fortran-side defaults creating inconsistencies.

**Implementation:**
- Fortran derived types (`solver_config`, `et_params`, `material`, etc.) declare fields with **no** initialisation values.
- The f2py wrapper subroutines accept all parameter values as mandatory arguments.
- `SolverConfig`, `ETParams`, etc. in `core.py` define all defaults via Pydantic `Field(default=...)`.
- The Python layer always passes explicit values to the kernel; Fortran never receives sentinels or undeclared values.

**Example violation to avoid:**

```fortran
! DO NOT DO THIS:
type :: solver_config
  real(dp) :: courant_number = 0.9_dp   ! Fortran-side default — FORBIDDEN
end type solver_config
```

**Correct pattern:**

```fortran
! DO THIS:
type :: solver_config
  real(dp) :: courant_number   ! No default; value set by Python layer
end type solver_config
```

```python
# And in Python:
class SolverConfig(Freezable):
    courant_number: float = Field(default=0.9, gt=0, le=1, ...)
```

When `kernel_set_solver_params()` is called, it always supplies an explicit value: `kernel_set_solver_params(0.9, ...)`.

---

## 9. Implementation notes

### 9.1. Unit handling

The Fortran kernel operates exclusively in SI units (metres, seconds). The Python layer converts all user-facing inputs to SI before passing to the kernel and converts outputs back on retrieval, ensuring the kernel is unit-agnostic and physics are dimensionally consistent.

### 9.2. Coupling protocol (future)

For coupling to external SW and GW models at different intervals:

1. SW fluxes are fast, so couple every macro-step: set `lat_sw` per step.
2. GW fluxes are slow, so couple at coarser intervals (e.g. every $N$ steps): accumulate `acc_recharge` over $N$ steps and feed to the GW model, then receive updated `lat_gw` for the next $N$ steps.
3. Convergence iteration: run GWSWEX with lateral fluxes, compare delta-GW and delta-SW between GWSWEX and the external models, iterate until below a convergence threshold.

### 9.3. `theta_i` (blended VWC at boundary layer)

For plotting the continuous VWC profile through the boundary layer, the Python layer computes on read-back:

$$\theta_i^{[l^*]} = \frac{\theta_s^{[l^*]} \cdot (t_{l^*} - \text{GW}) + \theta^{[l^*]} \cdot d_a^{[l^*]}}{t_{l^*} - b_{l^*}}$$

where $\theta^{[l^*]}$ is the UZ-derived VWC of the boundary layer's unsaturated portion. This is a display quantity only; it does not participate in the physics.

### 9.4. Pre-processing of the rooting mask

The rooting mask is derived inside `init()` from each element's vegetation type before the configuration is frozen: for static vegetation, layers whose midpoint depth is within `RootParams.depth` of the surface are flagged; for dynamic vegetation, the mask is initialised at the t=0 rooting depth and refreshed every macro-step from `_recompute_is_root_at_frac`. Vegetation types must supply either a static `RootParams(depth=...)` or both `root_depth_initial` and `root_depth_final` together with a non-static `root_growth_model`; missing or contradictory specifications raise inside `init()`.

### 9.5. Caching of time-varying parameters

When `update_is_root()` is called, it compares the new mask against the currently loaded array and makes the f2py call only if values have changed, avoiding unnecessary Fortran-side copies for infrequently updated parameters (e.g. seasonal root growth that crosses a layer boundary only a handful of times per year). The same pattern applies to any future time-varying parameters.

---

## 10. Documentation standards

### 10.1. Fortran kernel — FORD

All Fortran source files are documented for [FORD](https://github.com/Fortran-FOSS-Programmers/ford). The primary docmark is `!>` (FORD predocmark), placed **before** the documented entity; every line of the block uses `!>`. The `!!` post-placed docmark is reserved exclusively for derived-type component descriptions (field-trailing comments).

#### 10.1.1. Comment markers

| Marker | Placement     | Role                                                                                             |
| ------ | ------------- | ------------------------------------------------------------------------------------------------ |
| `!>`   | Before entity | **Primary docmark** — pre-placed doc-comment (all lines of the block)                            |
| `!!`   | After entity  | Post-placed docmark — used for derived-type component descriptions only                          |
| `!*`   | Before entity | Alternative pre-placed block-starter (`predocmark_alt`): subsequent plain `!` lines are absorbed |
| `!`    | Anywhere      | Ordinary (non-documentation) comment — ignored by FORD                                           |

Use `!>` for all module, subroutine, function, and derived-type doc-blocks. Use `!!` only for trailing component descriptions inside `type` blocks. Use plain `!` freely for inline implementation notes.

#### 10.1.2. Document formatting: No emoji characters

This document uses only plain text and standard Markdown formatting. **Never use emoji characters** (Unicode emoji like ❌, ✅, 🔧, 📝, etc.) in this Markdown file **or** in Fortran source documentation.

**ASCII symbols are permitted only in this Markdown documentation** for brevity when they improve readability. Examples of acceptable ASCII symbols in Markdown only:
- `→` (right arrow) for "leads to" or "produces"
- `←` (left arrow) for "receives from"
- `✓` or `•` for bullet lists or emphasis (in addition to Markdown bullets)
- `↑` / `↓` for "increase" / "decrease"
- `◇` or `○` for cross-references or markers

**CRITICAL: ASCII symbols are strictly forbidden in Fortran source files.** All Fortran comments and documentation must use only basic ASCII letters, numbers, and standard punctuation (no special Unicode characters at all).

For text alternatives in documentation, use plain-text instead of emoji:
- For "correct" or "preferred": Use bold text or explicit language (e.g., "**Correct pattern:**", "**Recommended approach:**")
- For "incorrect" or "avoid": Use bold text or explicit language (e.g., "**Example violation to avoid:**", "**Do not use:**")
- For emphasis or callouts: Use `**bold**`, `_italic_`, or Markdown inline code formatting

This ensures the document remains:
- **Accessible**: no platform rendering or encoding issues.
- **Professional**: consistent with formal technical documentation standards.
- **Machine-readable**: compatible with automated tooling and version control.
- **Portable**: no dependency on special character rendering libraries.

#### 10.1.3. Module-level doc-block

Place the doc-block **before** the `module` statement using `!>` on every line:

```fortran
!> Process-level subroutines implementing the solver cascade defined in
!> [model-physics.md S3](model-physics.md).  Each subroutine operates on a
!> single element index and mutates the model singleton in-place.
!>
!> author: <author>
!> version: 1.0
module gwswex_processes
  use gwswex_constants, only: dp, EPS
  ...
```

Required meta-data keys (placed at the very top of the doc-block, before any blank `!>` line):

| Key       | Required for      | Notes                               |
| --------- | ----------------- | ----------------------------------- |
| `author`  | all modules       | contributor name(s)                 |
| `version` | all modules       | matches `pyproject.toml` version    |
| `summary` | optional override | FORD uses first paragraph if absent |

#### 10.1.4. Subroutine and function doc-blocks

Place the entire doc-block **before** the `subroutine`/`function` statement using `!>` on every line. Use `@param` and `@result` tags inside the doc-block — do not add trailing `!!` to individual dummy argument declarations:

```fortran
!> Apply prescribed lateral GW and SW flux rates to element `ex` for
!> sub-step duration `dt` (S3.3).  GW is updated in volume-space and
!> GWH recovered via \(V_\text{GW}^{-1}\) only if the volume changes.
!>
!> @param M   Model singleton (mutated in-place).
!> @param ex  Element index (1-based).
!> @param dt  Sub-step duration [T].
subroutine apply_lateral(M, ex, dt)
  type(gwswex_model), intent(inout) :: M
  integer,            intent(in)    :: ex
  real(dp),           intent(in)    :: dt
```

For functions, document the return value with `@result`:

```fortran
!> Recover GW table elevation from drainable volume (S2.4).
!>
!> @param  V_gw  Drainable GW volume [L3].
!> @param  geo   Element geometry descriptor.
!> @result gwh   GW table elevation above datum [L].
pure function V_gw_inv(V_gw, geo) result(gwh)
  real(dp),          intent(in) :: V_gw
  type(geo_element), intent(in) :: geo
  real(dp)                      :: gwh
```

Rules:
- The summary (first `!>` line) must be a complete sentence ending with a period.
- Cross-reference the physics spec section in parentheses, e.g. `(S3.3)`.
- Use `@param <name>  <description>` for every dummy argument; two spaces between name and description.
- Use `@result <name>  <description>` for the function result variable.
- Do **not** add trailing `!!` lines to dummy argument declarations — `@param`/`@result` in the preceding doc-block is sufficient.
- Use `\( … \)` for inline LaTeX, `$$ … $$` or `\[ … \]` for display math.
- For `pure`/`elemental` functions, include units in square brackets in the `@param`/`@result` description: `@param Ks  Saturated hydraulic conductivity [LT-1].`

#### 10.1.5. Derived-type doc-blocks

Place the type's doc-block **before** the `type` statement with `!>`. Component descriptions use trailing `!!` (post-placed), since FORD has no `@param` equivalent for type fields:

```fortran
!> author: <author>
!> summary: Top-level model singleton holding all kernel state.
!>
!> All allocatable arrays are indexed column-major `(nl, ne)`.
!> State arrays use `*_prev`/`*_curr` pairs; no buffer-index dimension.
type :: gwswex_model
  integer :: ne = 0
    !! Number of model elements.
  integer :: nl = 0
    !! Number of layers per element.
  ...
```

Every public component must have a trailing `!!` description. Private components may omit it.

#### 10.1.6. Special environments

Use FORD's block tags where appropriate — they are rendered as styled boxes in the generated HTML:

```fortran
!> Resolve layer state transitions and enforce physical bounds after a
!> sub-step (S3.10).
!>
!> @param M   Model singleton.
!> @param ex  Element index (1-based).
!>
!> @note
!> Must be called *after* all process subroutines and *before* the
!> CFL check for the next sub-step.
!> @endnote
!>
!> @warning
!> Calls `update_geometry` internally; do not call it again immediately
!> afterwards without an intervening state change.
!> @endwarning
subroutine geometry_resolution(M, ex)
```

Available tags: `@note` / `@endnote`, `@warning` / `@endwarning`, `@todo` / `@endtodo`, `@bug` / `@endbug`. Tags are case-insensitive.

Argument and result tags (no closing tag needed):

| Tag                      | Use                                                      |
| ------------------------ | -------------------------------------------------------- |
| `@param <name>  <text>`  | Documents a dummy argument by name                       |
| `@result <name>  <text>` | Documents the function result variable                   |
| `@see <reference>`       | Points to related routines, types, or external resources |

#### 10.1.7. Cross-linking

Use FORD's double-bracket link syntax to cross-reference other documented entities inline:

```fortran
!> See [[update_geometry]] and [[V_gw_inv(function)]].
!> This subroutine is called from [[solve_element]].
```

Syntax: `[[name]]` for an unambiguous name; `[[name(type)]]` to disambiguate (e.g. `[[V_gw(function)]]`).

For referencing related routines, types, or external material use `@see`:

```fortran
!> Recompute layer geometry after a state change (S3.10).
!>
!> @param M   Model singleton.
!> @param ex  Element index (1-based).
!>
!> @see [[geometry_resolution]], [[gwswex_geometry(module)]]
!> @see model-physics.md Section 3.10
subroutine update_geometry(M, ex)
```

`@see` entries are rendered as a *See also* list in the HTML output. Multiple `@see` lines are allowed.

#### 10.1.8. What NOT to document inline

- Implementation comments (`! local loop index`, `! clamp to domain bottom`) — use plain `!`.
- Temporary variables local to a subroutine body — do not add `!!` to local variable declarations unless they carry non-obvious physics meaning.
- The f2py wrapper module and its subroutines — these are auto-generated interface plumbing; a single module-level `!!` block noting "f2py interface layer" is sufficient.

---

### 10.2. Python — Black + mypy

All Python source files (`core.py`, `model.py`, `io.py`, `utils.py`, `reference.py`) must conform to **Black** formatting and pass **mypy** strict type-checking.

#### 10.2.1. Black formatting

Enforce with:

```bash
black --line-length 100 gwswex/
```

Key conventions that Black enforces (do not fight them):
- Double quotes for all strings.
- Trailing commas in multi-line argument lists and collections.
- Magic trailing comma in function signatures to force per-argument line breaks.
- 100-character line limit (configured in `pyproject.toml`: `[tool.black] line-length = 100`).

#### 10.2.2. mypy type checking

Run with:

```bash
mypy --strict gwswex/
```

Required configuration in `pyproject.toml`:

```toml
[tool.mypy]
strict = true
ignore_missing_imports = true    # for f2py compiled module (f_gwswex)
```

All public functions and methods must carry fully annotated signatures:

```python
def step(
    self,
    dt: float,
    precip: np.ndarray,
    pet: np.ndarray,
    ptt: np.ndarray,
) -> None: ...
```

Use `np.ndarray` for all NumPy arrays. Where shape constraints matter, document them in the docstring rather than using runtime shape-check types.

#### 10.2.3. Docstring format

Use **Google-style** docstrings (compatible with both Sphinx and readable plain-text):

```python
def get_state(self) -> dict[str, np.ndarray]:
    """Retrieve the post-step model state in user units.

    Returns the state arrays after the most recently completed macro-step.
    Arrays are copies; mutating them does not affect kernel state.

    Returns:
        dict with keys:
            ``GWH``: GW table elevation, shape ``(ne,)``.
            ``GWV``: Drainable GW volume, shape ``(ne,)``.
            ``SW``:  Surface water depth, shape ``(ne,)``.
            ``UZ``:  UZ storage per layer, shape ``(nl, ne)``.
            ``theta``: Volumetric water content, shape ``(nl, ne)``.

    Raises:
        RuntimeError: If the kernel has not been initialised via ``init()``.
    """
```

Rules:
- Every public class, method, and module-level function must have a docstring.
- First line: single-sentence summary, no period for one-liners, period for multi-liners.
- Sections used: `Args:`, `Returns:`, `Raises:`, `Note:`, `Example:` — only include sections that are non-trivial (skip `Args:` if the function has no arguments).
- For Pydantic `Field(...)` descriptors: the `description=` string in the field is the authoritative doc; the class docstring describes the type's purpose and invariants only.
- Private methods (leading `_`) do not require docstrings unless the logic is non-obvious.

#### 10.2.4. What NOT to document

- Pydantic validator methods (`@field_validator`, `@model_validator`) — these are implementation detail; no docstring required.
- The `__init__` of Pydantic models — Pydantic generates this; document field semantics via `Field(description=...)` instead.
- Trivial property aliases (e.g. `@property` that simply returns `self._x`).
- Type alias definitions (`MyType = dict[str, np.ndarray]`) — inline comment is sufficient.

---

## 11. Tuning parameters

User-facing knobs that calibrate a GWSWEX simulation against a reference
without changing the underlying physics. The physical role of each empirical
parameter is documented in
[model-physics.md](model-physics.md#empirical-parameters-checklist-sensitivity-analysis-candidates);
this section lists where they are surfaced in the API, their typical ranges,
and the direction of their effect.

They split into three families:

- **Discretisation** (quasi-tuning): set numerical accuracy and cost. Changing
  these does not alter the governing equations, but truncation error is not
  negligible for coarse choices, so they must be swept alongside the
  solver-specific knobs.
- **Model parameters** (`set_model_params`): shape the Green-Ampt infiltration
  closure and the explicit cascade's connectivity floor. Physically motivated,
  but their numerical values are effectively regularisers and can be calibrated.
- **Solver parameters** (`set_solver`): control the internal sub-stepping,
  quadrature, hysteresis blending, and Picard convergence. Different between
  the two solvers.

A worked sensitivity analysis against HYDRUS-1D, including recommended
calibrated configurations for a deep-WT and a shallow-WT setup, is in
[docs/sensitivity-analysis-from-comparison.md](sensitivity-analysis-from-comparison.md).

### 11.1 Discretisation

| Parameter | Where | Typical range | Effect direction | Notes |
|-----------|-------|---------------|------------------|-------|
| `DZ` (layer thickness) | `init_space` | 0.01–0.05 m | Finer → lower truncation error; cost scales $\propto 1/DZ$ | 0.02 m is a good Pareto point for loam at field scale; below 0.01 m no further accuracy gain. |
| macro-step `dt` | `init_time` | explicit: forcing-cadence; implicit: 1 h–several h | Shorter → lower time-splitting error; cost scales $\propto 1/dt$ | Implicit tolerates larger macro-steps than the forcing cadence because its Picard iteration smooths internal sub-dynamics. |
| `dt_min` | `init_time` | seconds to a few minutes | Sets adaptive-step floor (explicit) and Picard fallback (implicit) | Only relevant when adaptive sub-stepping triggers. |

### 11.2 Model parameters (`set_model_params`)

| Parameter | Range | Effect direction | Notes |
|-----------|-------|------------------|-------|
| `psi_f` (Green-Ampt suction head) | 0.005–0.20 m | Higher → larger infiltration capacity at front | Only active when $P > K_\text{sat}$; inert on loam under moderate forcing. |
| `F_min` (min cumulative infiltration) | $10^{-8}$ to $10^{-4}$ m | Regulariser for the Green-Ampt capacity denominator | Only nudges behaviour at the onset of a storm; inert on loam. |
| `ICratio_min` (min inter-layer connectivity) | 0.01–0.50 | Higher → faster vertical redistribution in the explicit cascade | Noticeably affects explicit WT drawdown rate; implicit solver is insensitive. |

### 11.3 Explicit-solver parameters (`set_solver(solver="explicit", …)`)

| Parameter | Range | Effect direction | Notes |
|-----------|-------|------------------|-------|
| `courant_number` | 0.3–0.95 | Lower → more sub-steps, less splitting error, higher cost | Default 0.9 is fine for loam; very stiff soils may need 0.5. |
| `n_trapz` | 5–40 | Trapezoidal quadrature nodes for $UZ_\text{eq}$ | 10 is sufficient for loam; 20 is default; over-resolving a near-linear retention profile gives no benefit. |
| `beta_hyst` | 0.7–1.0 | Higher → more-imbibition-like retention in draining zones | 1.0 recovers the van Genuchten main-drainage curve; 0.85 mimics a hysteretic wetting branch. |

### 11.4 Implicit-solver parameters (`set_solver(solver="implicit", …)`)

| Parameter | Range | Effect direction | Notes |
|-----------|-------|------------------|-------|
| `picard_tol` | $10^{-7}$ to $10^{-4}$ m | Lower → tighter head convergence per macro-step | Default $10^{-5}$ under-converges the capillary fringe for deep WT cases; $10^{-6}$ is safe. |
| `picard_max_iter` | 50–500 | Headroom for the tolerance | Raise together with a tightened `picard_tol`. |
| `beta_hyst` | 0.7–1.0 | Same role as explicit | Only applied when computing diagnostic $\theta$ from the head solution. |
| `n_trapz` | 10–40 | Only used for IC construction from hydrostatic equilibrium | Does not affect the transient. |

### 11.5 Parameters **not** to calibrate

These change the physics and must match the reference experiment, not the solver:

- van Genuchten–Mualem retention: `alpha`, `n`, `theta_r`, `theta_s`, `lam`
- Saturated hydraulic conductivity `K_sat`
- Root-growth schedule

### 11.6 Per-soil OAT calibration: lessons learned

A one-at-a-time (OAT) sweep across all twelve soil × setup × solver combinations
of the multi-soil HYDRUS-1D benchmark
([sensitivity-analysis-from-comparison.md](sensitivity-analysis-from-comparison.md))
yields the following practical guidance:

- **`ICratio_min` is the dominant lever for the explicit cascade.** Its optimum
  is highly soil-dependent: clay-rich profiles benefit from low values
  ($0.05$); near-pure loam wants high values ($0.5$); sand and layered profiles
  sit in the middle ($0.2$–0.42). Implicit results are insensitive — leave at
  default.
- **`courant_number` matters most for low-$K_\text{sat}$ soils.** In the
  basic-setup clay benchmark, raising `courant_number` from 0.9 to 0.95
  removes a residual splitting bias in the WT envelope.
- **`n_trapz` saturates by $\sim 20$ for the explicit cascade.** Going from
  10 to 20 visibly improves clay and loam-clay; further refinement is
  numerical noise.
- **`picard_tol` is the only meaningful implicit knob in field-scale runs.**
  Tightening from $10^{-5}$ to $10^{-6}$ or $10^{-7}$ closes a small
  late-recession bias in basic-setup loam; pair with `picard_max_iter` $\geq 200$
  to avoid non-convergence.
- **ET-stress thresholds are calibratable within their physical bounds.** The Laio (2001) thresholds `s_star`, `s_w`, `s_h`, `s_e` carry physical meaning but their optimal values are sensitive to the soil–plant system. Multi-soil benchmarks against HYDRUS-1D show that joint sweeps over these four parameters (bracketed by `s_h ≤ s_w ≤ s_e ≤ s_star ≤ 1`) reduce WT-RMSE by 5–30 % across both solvers. They should be included in any OAT calibration; per-case calibrated values are tabulated in [sensitivity-analysis-from-comparison.md](sensitivity-analysis-from-comparison.md).
- **`beta_hyst` and `psi_f` are essentially inert** at the OAT scale tested here; leave at defaults unless the reference experiment imposes a specific hysteretic branch.
- **The GWSWEX–HYDRUS-1D pre-canopy evaporation accounting differs structurally.** GWSWEX applies potential soil evaporation and transpiration as separate prescribed rates (wet-phase fluxes pass through to the soil surface unchanged), whereas HYDRUS-1D partitions canopy interception before entering Richards. This is a known structural difference rather than a calibration target; see [sensitivity-analysis-from-comparison.md §3.2](sensitivity-analysis-from-comparison.md) for details.
- **Structural ceilings persist.** For the intensive-setup clay and loam
  cases, no combination of tuning knobs closes the WT trajectory to within
  $\sim 30$ cm of HYDRUS-1D; this reflects the discretisation gap between
  the GWSWEX layered-bucket cascade and the HYDRUS-1D Richards solver under
  rapidly-cycling deep-WT conditions, not a calibration deficit. See
  [sensitivity-analysis-from-comparison.md §3](sensitivity-analysis-from-comparison.md)
  for the full numerical record.

---

## Appendix B: Historical V9 architecture record

> Consolidated from `docs/model-arch-v9.md`.
> The V9 decisions described here are the basis of the current implementation
> and are documented for traceability.

### B.1 Variant selection rationale

The V1–V9 diagnostic campaign (archived under `.archive/model-diagnosis/`)
evaluated 26 scheme variants. Headline results:

| Variant | GW-RMSE | R² | MB_max | Runtime |
|---|---|---|---|---|
| **V9 implicit Picard** | **1.1 cm** | **0.998** | 3.52 mm | **8.5 s** |
| V4.St (Strang + CFL-B, best explicit) | 28.2 cm | −0.62 | 1.49 mm | 268 s |
| V_best (best explicit ensemble) | 27.0 cm | −0.48 | 1.54 mm | 434 s |
| V1.St (Strang + CFL-A) | 35.2 cm | −1.52 | 14.95 mm | 7.3 s |

The implicit V9 solver is approximately 25× more accurate in water-table tracking
than the best explicit variant and typically faster.

### B.2 Physics extensions over the V9 barebone

The Python barebone V9 (in `.archive/model-diagnosis/barebone/`) established the
core physics. Four extensions are present in the production Fortran kernel:

| Extension | Status | Rationale |
|---|---|---|
| Adaptive time-stepping (Picard-count proxy) | Included | Robustness during sharp wetting fronts; larger steps during dry periods. |
| ET stress function (Laio 2001) | Included | Realistic transpiration under water stress. |
| Root water uptake (distributed layer sink) | Included | Column-scale ET partitioning. |
| Newton-Raphson linearisation option | Deferred | Picard converges robustly; Newton was 2× slower in V1.5 experiments. Solver interface remains abstract for future Newton add. |
| Hysteresis in retention curve | Excluded | Mualem-VG with warm-start head profile captures dominant wetting–drying history. |

### B.3 Third-party Fortran libraries

**FOODIE** (Fortran Object-Oriented Differential-equations Integration
Environment) was evaluated. After finite-volume spatial discretisation the
Richards system is a tridiagonal DAE; the optimal solver is a custom Thomas
(TDMA) algorithm in $O(n)$. FOODIE's implicit integrators do not exploit
tridiagonal structure and would require user-supplied Jacobian-vector products.
**Not used.**

**LAPACK** is not required: LAPACK's `dgtsv` for tridiagonal systems is
asymptotically equivalent to the custom TDMA; no external dependency is
warranted. **Not used.**

The kernel has no external numerical-library dependencies. All constitutive
relations (VG, Mualem, Laio), the TDMA solver, and the Picard iteration are
self-contained.

### B.4 Solver-abstraction design

The current `gwswex_solver_implicit` module implements the Celia (1990)
mixed-form Picard iteration directly. The original V9 design called for an
abstract `abstract_solver_t` type with deferred bindings and parallel
`gwswex_solver_picard` / `gwswex_solver_newton` modules. This was deferred:
Newton offered no demonstrable convergence-rate advantage in the V1.5 Python
experiments (often 2× slower with Jacobian issues at harmonic-mean K
interfaces), and a single-implementation module is simpler. The dispatcher
field `solver_config%solver` (`'explicit'` vs `'implicit'`) lives in
`gwswex_kernel` and selects the correct solver subroutine at the macro-step
boundary; adding a third solver in future would only require a new module and
a dispatch branch.
