# GWSWEX — Consolidated TODO

Single source of truth for forward-looking work. Items here come from `.agent/implementation_plan.md`, in-code `! reserved for ...` comments, the documentation audit performed alongside Phase 15, and the BMI / coupling scaffolding files just added.

For *historical* notes on completed phases see [`.agent/implementation_notes.md`](.agent/implementation_notes.md).

Conventions: `[ ]` not started · `[~]` in progress · `[x]` done · `[!]` blocked.

---

## §1 — BMI integration (forward-looking scaffolding)

Origin: [`gwswex/coupler.py`](gwswex/coupler.py) (`BmiGwswex`) and [`gwswex/src/kernel_bmi.f08`](gwswex/src/kernel_bmi.f08). Both files are pure stubs — every non-trivial method raises `NotImplementedError` (Python) or `error stop` (Fortran). Reference: [BMI 2.0 spec](https://bmi.csdms.io/en/stable/bmi.spec.html).

Implementation order (each step depends on the previous):

- [ ] **1.1 Control group.** Implement `initialize(config_file)`,
`update()`, `update_until(time)`, `finalize()` in `BmiGwswex`. Map onto the existing `GWSWEXmodel.init` / `step` / `deinit` lifecycle. Decide on YAML vs TOML for the config file (recommend TOML; it is in the stdlib from Python 3.11).
- [ ] **1.2 Time group.** Implement `get_current_time`, `get_start_time`,
`get_end_time`, `get_time_step`. The trivial `get_time_units → "s"` is already wired.
- [ ] **1.3 Variable info group.** For each entry in `_input_var_names` /
`_output_var_names` populate `get_var_grid`, `get_var_type`, `get_var_units`, `get_var_itemsize`, `get_var_nbytes`. CSDMS Standard Names already chosen.
- [ ] **1.4 Getters / setters.** `get_value`, `get_value_ptr` (with care —
zero-copy view requires NumPy `np.ascontiguousarray` round-trip from the kernel), `get_value_at_indices`, `set_value`, `set_value_at_indices`. Map to the existing `GWSWEXmodel.update_forcing` and `get_state` paths.
- [ ] **1.5 Grid functions.** Decide grid type per output: ponding/ET as
`scalar` per element, profile θ/UZ as `rectilinear` over (z, ne). Implement `get_grid_rank`, `get_grid_size`, `get_grid_shape`, `get_grid_x` (= z midpoints from `bnds`), `get_grid_node_count`. The non-rectilinear surface-mesh functions (`get_grid_face_*`, `get_grid_edge_*`) can keep raising `NotImplementedError` indefinitely.
- [ ] **1.6 Fortran BMI.** Implement the matching `bmi_*` procedures in
[`gwswex/src/kernel_bmi.f08`](gwswex/src/kernel_bmi.f08). Wire into [`meson.build`](meson.build) only after the Python side is stable, since every Fortran rebuild adds CI cost.
- [ ] **1.7 Add `tests/test_bmi.py`.** Cover at least: initialize/update/finalize
round-trip; `get_value` / `set_value` for each variable in `_input_var_names` and `_output_var_names`; grid accessor consistency.
- [ ] **1.8 Wire `gwswex.coupler.BmiGwswex` into `gwswex/__init__.py`** once
the control + time + getter groups are functional.

## §2 — GW / SW coupling adapters

Origin: [`gwswex/coupler.py`](gwswex/coupler.py) (`GWCoupler`, `SWCoupler`, `CoupledExchange`). Forward-looking; no concrete partner model yet.

- [ ] **2.0 Optional spatial awareness layer (shared GW/SW grid).** GWSWEX is,
and shall remain, *physically* spatially agnostic — every element is solved in isolation with no inter-element fluxes. The optional layer is purely a **bookkeeping / addressing** facility: accept at `init_space` an optional `grid` descriptor (the *same* horizontal grid used by both the partner GW and SW models) so that each GWSWEX element carries an unambiguous spatial identity. This makes GW ↔ GWSWEX ↔ SW exchange native and intrinsic (mapping is by grid index, no external reindexing needed) without changing any kernel physics. Sub-tasks:
    - [ ] **2.0.1 Grid descriptor data model.** Add a `Grid` Pydantic type in
[`gwswex/config.py`](gwswex/config.py) supporting at least the BMI grid types `scalar`, `points`, `uniform_rectilinear`, `rectilinear`, and `unstructured` (face-centred). Carry rank, shape, origin/spacing or x/y/z arrays, and optional `crs` (EPSG / WKT). Frozen at `register_space()`.
    - [ ] **2.0.2 Element ↔ grid mapping.** Add a required `cell_id(ne)`
array (1-based) when a `grid` is supplied; validate that `1 ≤ cell_id ≤ grid.size` and that ids are unique. Optionally accept `cell_id` as a structured `(i, j)` / `(i, j, k)` index for rectilinear grids.
    - [ ] **2.0.3 Round-trip through I/O.** Persist the grid descriptor and
`cell_id` mapping in NetCDF outputs ([`gwswex/io.py`](gwswex/io.py)) and in checkpoints, with CF-1.8 grid conventions. Restart must reconstruct the mapping.
    - [ ] **2.0.4 Coupler integration.** `GWCoupler` / `SWCoupler` become
thin index-mappers: `pull_lateral` reads the partner's per-cell flux at `cell_id[ex]`; `push_recharge` / `push_runoff` write to the same. Eliminates the per-coupling-pair boilerplate that would otherwise be needed for every partner combination.
    - [ ] **2.0.5 BMI exposure.** When a grid is registered, `BmiGwswex`
`get_grid_*` returns the *partner* grid (not the synthetic `(z, ne)` profile grid), so couplers see GWSWEX as a co-located, grid-native component. The profile / column grid remains available under a second BMI grid id for vertical-state outputs.
    - [ ] **2.0.6 Backward compatibility.** When no `grid` is supplied,
behaviour is unchanged: `cell_id` defaults to `1..ne` and the synthetic `points` grid is reported via BMI. All existing tests and examples must continue to pass without modification.
    - [ ] **2.0.7 Worked example.** Add a short `examples/spatial_grid.py`
showing a 4 × 3 uniform-rectilinear grid mapped to 12 GWSWEX elements, with one or two cells deactivated (no GWSWEX element) to demonstrate sparse occupancy.

Rationale: the kernel stays physically agnostic (no neighbour stencils, no inter-element coupling, no horizontal fluxes inside GWSWEX), but *coupling becomes native* because both partners and GWSWEX share an externally agreed grid topology. This is the same pattern used by PRMS-MODFLOW6 (GSFLOW), CSDMS BMI grid contracts, and the OpenGeoSys surface-coupling layer.

- [ ] **2.1 Pick a first concrete GW partner.** Recommended: MODFLOW 6 via
`flopy` (BMI-aware, well-documented, used by many similar coupling efforts).
- [ ] **2.2 Pick a first concrete SW partner.** Recommended: SWMM via
`pyswmm` for plot/network-scale, or a kinematic-wave reference for purely overland coupling.
- [ ] **2.3 Sequential ("loose") coupling.** Implement `pull_lateral` →
`step` → `push_recharge` / `push_runoff` for one macro-step. Exchange the full `CoupledExchange` payload at every macro-step.
- [ ] **2.4 Iterative ("strong") coupling.** Add a "rewind" entry point in the
kernel that restores `*_prev` arrays without recomputing geometry, so the coupling outer iteration can re-issue a step with updated lateral fluxes. Required only when partner models also need to iterate to convergence; not on the critical path for §2.3.
- [ ] **2.5 Document coupling patterns** in a new `docs/coupling.md` once at
least one partner is wired up end-to-end.

## §3 — Kernel reserved hooks

Two `! reserved for ...` placeholders are in the current source. Neither blocks any current functionality.

- [ ] **3.1 Rate-form mass-balance outputs.**
[`gwswex/src/shared/mass_balance.f08`](gwswex/src/shared/mass_balance.f08) line 131 (`associate(unused_dt => dt_macro)`). Currently all MB diagnostics are stored as cumulative volumes per macro-step. Add an optional rate-form output (volume / dt_macro) for partners that consume per-step fluxes natively (most BMI couplings; MODFLOW 6 stress periods).
- [ ] **3.2 Geometry-correction diagnostics.**
[`gwswex/src/explicit/solver.f08`](gwswex/src/explicit/solver.f08) line 116. Currently the explicit cascade silently corrects layer storages when geometry refresh detects a water-table crossing; expose the per-substep correction magnitude through `get_accumulators` for diagnostics.

## §4 — Configuration refinements

- [ ] **4.1 Real exponential root-growth law.**
[`gwswex/config.py`](gwswex/config.py) line ~153 — `RootGrowthModel` currently accepts `"exponential"` as a label but the Python loop (`_recompute_root_weights_at_frac`) silently falls back to linear interpolation. Implement the genuine exponential law `d(t) = d_max · (1 − exp(−k · t/T))` with `k` exposed as a vegetation parameter. Currently documented as a known stub in [`docs/model-arch.md` §5.3.3](docs/model-arch.md#533-vegetation-types-and-root-growth-models); remove the "reserved for a future genuine exponential growth law" caveat once shipped.

## §5 — Documentation

- [ ] **5.1 Write `docs/coupling.md`** documenting sequential and iterative
coupling patterns once §2.3 / §2.4 are functional.
- [ ] **5.2 Write `docs/bmi.md`** with `BmiGwswex` usage notes, CSDMS Standard
Names used, and the grid id scheme once §2.1–2.5 are functional.

## §6 — Testing

- [ ] **6.1 Cross-solver checkpoint guard test.** Save a checkpoint under one
solver, attempt to load it under the other in a fresh process; confirm the `RuntimeError` from `load_checkpoint` fires before any kernel mutation. (The guard is implemented; the test currently lives only in the Phase 13 integration notes.)
- [ ] **6.2 Re-baseline `tests/test_kernel.py`** when the rate-form MB output
(§4.1) lands.

## §7 — Priority structural extensions (from manuscript §Limitations)

Origin: manuscript section *Limitations and priority extensions*. Listed here for tracking; expansion of any one requires a dedicated design note before implementation.

- [ ] **7.1 Richer evapotranspiration partitioning.** The dry-phase RMSE
reductions of 5 to 30 percent observed under Laio-threshold tuning in the OAT campaign indicate that the dry-phase residual against HYDRUS-1D is dominated by the ET parameterisation. Candidate replacements for the current Laio four-threshold closure: (i) a sub-daily energy-balance closure driven by net radiation, ground heat flux and aerodynamic resistance; (ii) a canopy interception store with throughfall and stem drainage; (iii) a coupled energy-and-water vegetation scheme with prognostic LAI and stomatal conductance. Priority near-term target.
- [ ] **7.2 Improved capillary-rise representation.** The cascade resolves
upward Darcian flux through the inter-layer connectivity ratio plus the equilibrium-storage quadrature, both of which degrade under sharp wetting fronts in low K_sat media. Candidate fixes are enumerated in `docs/issues.md` and in the open-source repository's issue tracker. Priority structural target after 7.1.
- [ ] **7.3 Guided UZ-discretisation utility.** Recommend, per element, a
layer count and per-layer thickness from the local soil-hydraulic structure, the climatological GWH range, and the expected GW-UZ-SW exchange intensity. Goal: relieve the user of the discretisation choice that currently has to be made by inspection and that drives the accuracy-versus-cost trade-off identified in the performance benchmark. Operational target.
- [ ] **7.4 Process closures not yet represented.** Lateral UZ flow,
retention-curve hysteresis, snow accumulation and melt, canopy interception, and non-equilibrium / preferential flow. Each is an independent extension and should be scoped on its own merits.

## §8 — Per-element solver dispatch (mixed explicit/implicit within a single run)

Background: the current architecture assigns a single solver to the entire model instance for the duration of a macro-step. `switch_solver()` can change the solver between simulation periods but applies the change globally to all elements. A genuinely heterogeneous regional domain (e.g. floodplain wetlands alongside deep-GWH uplands within the same GWSWEX domain instance) would benefit from the ability to route individual elements through the explicit solver while others use the implicit solver *within the same macro-step*, so that CFL-adaptive sub-stepping is incurred only where the physics demands it.

Design questions that must be answered before any implementation:

1. **State isolation.** The explicit and implicit solvers share the same kernel
state arrays (`uz`, `gw`, `sw`, `ic`, etc.) but carry disjoint working arrays (`f_ga`, `icratio` for explicit; pressure-head work arrays and Picard counters for implicit). A mixed-dispatch step must ensure that the per-element working arrays for the inactive solver are either never written or are safely ignored on the next step, and that `switch_solver()` still produces a consistent warm-start when applied globally after a mixed step.
2. **Kernel loop restructure.** `kernel_step` in `gwswex/src/kernel.f08`
currently dispatches all `ne` elements to either `solve_element_explicit` or `solve_element_implicit` via a single conditional on the `solver` flag. Mixed dispatch requires a per-element solver mask (e.g. `integer(int8) :: solver_mask(ne)` in kernel state), the OMP loop body branching on `solver_mask(ex)`, and a guarantee that the two branches are still OpenMP-safe (no false sharing between elements running different branches).
3. **Python API.** Decide the exposure surface:
   - Option A: `set_solver(solver="mixed", mask=np.array([...]))` — explicit
mask; simplest but requires the user to manage the mask across checkpoint save/restore.
   - Option B: `set_solver(solver="auto", threshold={"ponding_depth": 0.5,
"gwh_depth": 2.0})` — auto-routing each macro-step based on per-element state diagnostics; more powerful but adds a per-step classification cost and must handle the transition consistently at the warm-start boundary.
   - Option C: both, with mask overrideable by the user and auto-routing as
the default. Recommend exploring C.
4. **Checkpoint compatibility.** The `solver_mask` must be persisted in
checkpoints and validated on load; a checkpoint saved under a mixed configuration must restore correctly and must not confuse the existing cross-solver guard in `load_checkpoint`.
5. **Mass-balance accounting.** Verify that `acc_recharge`, `acc_runoff` and
the per-element mass-balance residual remain correct when elements in the same step use different solvers.
6. **OAT / sensitivity implications.** A mixed-dispatch run is no longer
directly comparable to a single-solver run on the same forcing and soil; the sensitivity analysis infrastructure and any future GSA harness must be able to record the per-element solver assignment alongside the parameter sweep metadata.

Suggested exploration order:

- [ ] **8.0.1 Feasibility audit.** Read `gwswex/src/kernel.f08` and
`gwswex/src/explicit/solver.f08` / `gwswex/src/implicit/solver.f08` to inventory all arrays that are solver-specific vs shared; map the working arrays that would need to be allocated per-element under mixed dispatch. Write up findings as a design note in `.agent/implementation_notes.md`.
- [ ] **8.0.2 Prototype kernel change.** Add `solver_mask(ne)` to kernel
state; modify the `kernel_step` OMP loop to branch per element; confirm bit-equivalence with the existing single-solver paths when the mask is all-explicit or all-implicit. No Python API change at this stage.
- [ ] **8.0.3 Python API draft.** Implement Option A (explicit mask) in
`GWSWEXmodel.set_solver` and `switch_solver`; expose `solver_mask` through `get_state` and persist in checkpoint.
- [ ] **8.0.4 Auto-routing heuristic.** Implement Option B auto-routing
and benchmark against all-explicit and all-implicit runs on the 12 comparison cases to quantify the accuracy/cost trade-off.
- [ ] **8.0.5 Tests and documentation.** Add `tests/test_mixed_dispatch.py`;
update `docs/model-arch.md` and `docs/coupling.md`.

---

## §9 — OpenMP and Python-Fortran call-overhead optimisation

Motivation: observed CPU utilisation on the element loop is well below the theoretical maximum (saturation curve in the computational-performance benchmark turns over at ~4 to 6 threads on a 10-core M-series CPU even when ne is large). Two distinct overheads contribute: (a) work imbalance and fork-join cost inside the OMP parallel regions in `gwswex/src/kernel.f08`; (b) Python-side overhead in the per-step loop in `GWSWEXmodel.run` / `GWSWEXmodel.run_step`, which makes one f2py call per macro-step and several NumPy slice copies per call.

Each item below is independent and can be worked in any order. Expected gains are stated where they can be reasoned from the existing code; all must be confirmed by re-running `examples/gwswex-vs-hydrus1d/computational-performance-benchmark` and the 12 comparison notebooks before being declared done.

### 9.1 Element-loop scheduling

- [ ] **9.1.1 Switch the explicit-solver element loop to dynamic or guided
scheduling.** Lines: `gwswex/src/kernel.f08` ~line 360 (`!$omp parallel do schedule(static)` around `solve_element_explicit`). The CFL-adaptive sub-stepping inside `solve_element_explicit` produces per-element workloads that vary by one to two orders of magnitude (ponded / near-saturated columns take many sub-steps; dry columns take one). A static schedule therefore leaves the long-running threads waiting at the barrier. Recommend `schedule(dynamic, chunk)` with `chunk` tuned to a few elements (or `schedule(guided)` if profiling shows that dynamic overhead is non-trivial). Confirm bit-equivalent output; the per-element algorithm is unchanged.
- [ ] **9.1.2 Investigate the implicit-solver element loop.** Lines:
`gwswex/src/kernel.f08` ~line 366 (`solve_element_implicit`). The Picard iteration count varies by element: dry columns converge in one iteration, sharp-front columns may hit `picard_max_iter`. Same fix as above. Bench with at least one ponding-prone setup (intensive-clay) before merging.

### 9.2 Reduce parallel-region fork-join count

- [ ] **9.2.1 Fuse `solve_element_*` and `output_calc` into a single
parallel region.** Lines: `gwswex/src/kernel.f08` ~line 360 to ~line
  378. Currently each `kernel_step` opens, closes, and re-opens an OMP
parallel region (solver loop, then `output_calc` loop). On a 10-core CPU with small ne the fork-join overhead is several microseconds per fork; for ne ~ 100 and macro-steps ~ 10000, this accumulates. Wrap both loops in a single `!$omp parallel`, then use `!$omp do` for each, with an implicit barrier between them where required.
- [ ] **9.2.2 Audit branch divergence inside `solve_element_*`.** The
two solvers carry many regime branches (saturated column vs unsaturated, ponding vs dry, lateral inflow / outflow, ET regime). Element-level parallelism keeps these branches local to each thread, so divergence does not hurt correctness, but it does prevent SIMD inside the per-element loops where the inner loops over layers might vectorise. Identify hot inner loops in `gwswex/src/explicit/processes.f08` and `gwswex/src/implicit/solver.f08` that scan over `nl`; mark vectorisation-safe loops with `!$omp simd` (or restructure to remove branches from the layer loop body) and confirm under `-fopt-info-vec` whether the compiler accepts them.

### 9.3 Memory layout and false sharing

- [ ] **9.3.1 First-touch initialisation of per-element arrays.** Lines:
arrays allocated in `kernel_init` (`gwswex/src/kernel.f08`). On multi-socket / NUMA hosts the master thread allocates and zero-fills every element-indexed array, which pins all pages on one NUMA node. Add a `!$omp parallel do` initialisation loop immediately after each allocation to force first-touch placement onto the thread that will own that element during stepping. Modest gain on a single-socket M-series laptop; substantial on dual-socket EPYC / Xeon nodes.
- [ ] **9.3.2 False-sharing audit on accumulator arrays.** Lines:
`gwswex/src/kernel.f08` `acc_recharge`, `acc_runoff`, `acc_lat_*`, `acc_delta_*`, `n_substeps` (all length-`ne`). Each thread writes to its own element index, so semantically there is no sharing; but if `ne` is small and elements within a single 64-byte cache line are written by different threads, false sharing can occur. Verify chunk sizes and consider padding the per-element accumulators when ne is modest.

### 9.4 Reduce Python-Fortran call overhead in the per-step loop

- [ ] **9.4.1 Provide a `kernel_run` entry point that loops over
macro-steps in Fortran.** Lines: `gwswex/src/kernel.f08` (new subroutine `kernel_run(n_steps, dt, precip, pet, ptt, lat_gw, lat_sw)`); `gwswex/wrapper.f90` (matching `run` wrapper); `gwswex/model.py` (`GWSWEXmodel.run` switches to single batched f2py call when no per-step callback is supplied). Eliminates `n_steps` round-trips through f2py and the associated NumPy slice copies. Per-step callback and writer paths must continue to work via the existing `run_step` loop; `run` should auto-select the batched path only when it is safe (no callback, no dynamic root growth, no lateral fluxes that change per step). Largest expected gain for long-horizon, single-element / small-ne runs; the comparison notebooks fall in this regime.
- [ ] **9.4.2 Pre-bind forcing arrays at `init` time.** Lines:
`gwswex/wrapper.f90` (new `bind_forcing(precip, pet, ptt, lat_gw, lat_sw)` that stores pointers into the kernel state); `gwswex/model.py` (`set_forcing` writes once and then `run_step` advances an internal step index). Avoids per-step NumPy ascontiguousarray and slice copies. Compatible with 8.4.1; either both or one or the other can land.
- [ ] **9.4.3 Suppress per-step `get_state` / mass-balance accounting in
performance-sensitive runs.** Lines: `gwswex/model.py` `run` and `run_step`. `run_step(track=False)` already exists; ensure `run` exposes the same control (`track=True` default) so that the OAT harness and other long-batch workflows can disable per-step diagnostic copies, and document the behaviour in `docs/model-arch.md`.
- [ ] **9.4.4 Avoid repeated f2py module attribute lookups inside the
hot loop.** Lines: `gwswex/model.py` `step` method. The lookup `_F.gwswex_wrapper.step` is performed every macro-step; bind it once to a local in `init` (e.g. `self._step_fn = _F.gwswex_wrapper.step`) and call through the cached attribute. Negligible per call but cumulative for very long runs.

### 9.5 Profiling and verification

- [ ] **9.5.1 Add a profiling target.** A script under
`examples/gwswex-vs-hydrus1d/computational-performance-benchmark` that runs an OMP scaling sweep at ne ~ 10, 100, 1000 and reports parallel efficiency. Use the result to tune the `chunk` size in 8.1.1 and to confirm any scheduling change.
- [ ] **9.5.2 Bit-equivalence guard.** Any change in 9.1 to 9.4 must
reproduce `examples/demo-explicit/demo-explicit.nc` and `examples/demo-implicit/demo-implicit.nc` to machine precision and must pass all 12 comparison notebooks unchanged. The OAT harness outputs are a secondary check.
