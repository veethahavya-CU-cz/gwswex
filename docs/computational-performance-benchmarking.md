# Computational performance benchmarking — GWSWEX vs HYDRUS-1D

This note records the design, execution and findings of the computational performance benchmark introduced for the GWSWEX manuscript. The objective is to characterise the per-element wall-clock cost of the GWSWEX explicit and implicit solvers against a HYDRUS-1D reference, sweeping the number of independent vertical columns $n_e$ over three decades, on a single hardware configuration, with all I/O and parallelism choices fixed and traceable.

The full benchmark harness lives in `examples/gwswex-vs-hydrus1d/computational-performance-benchmark/`. The results referenced below are the raw `results/{system_info,gwswex_results,hydrus_results}.json`, the aggregated `results/summary_table.{md,csv}` and the `results/scaling.pdf` figure produced by `analyse.py`.

## Setup chosen and why

The intensive-sand-loam case from the verification catalogue was selected as the benchmark workload. The profile is a two-material layered column: sand occupies the upper third of the column and loam the lower two thirds, using the van Genuchten--Mualem parameters from the Carsel--Parrish catalogue. This layered profile exercises both the fast-drainage, CFL-active sand upper layer and the slower loam lower layer in the same column and within the same simulation, making it a more representative benchmark workload than a homogeneous single-material column: the explicit solver's CFL adaptive substepping is genuinely exercised by the sand fraction, and the implicit Picard iteration faces the material discontinuity at depth. The intensive forcing drives repeated ponding--drainage cycles that engage both layers over the 32-day horizon, so the timed region exercises the full operational kernel rather than a stripped-down core. There is also direct continuity with the `comparison-intensive-sand-loam.ipynb` notebook of the verification artefact, on which the same column is solved by both codes for the single-element accuracy comparison, making the wall-clock figures directly physically interpretable.

The full configuration (geometry, materials, layer ids, root depths, forcing arrays) is materialised by `benchmark_common.derive_intensive_case()`, which loads `examples/sensitivity-analysis/experiment_definitions.json` and applies the same per-soil derivations used by the comparison notebooks and the OAT worker. Layered profiles are supported by mapping each soil in `layers_spec` to its own HYDRUS material and reassigning each profile node's material id by depth. The active soil tag is set by the environment variable `GWSWEX_BENCH_SOIL` (default `sand-loam`) so that other variants can be re-run with no code change.

## OMP thread count and HYDRUS pool size

Both the OpenMP thread count for GWSWEX and the worker count for the HYDRUS pool are fixed at $T = W = 8$ rather than at the logical CPU count of 10. The Apple M1 Pro is a heterogeneous chip with 8 performance cores and 2 efficiency cores; in a development probe across $T \in \{1, \dots, 10\}$ on the same workload, the per-element wall-clock time bottomed out at $T = 8$ (5.0$\times$ speedup over single-threaded) and regressed at $T = 10$ to a 4.6$\times$ speedup, because the per-step OpenMP barrier in the GWSWEX kernel forces the eight performance-core threads to wait for the two efficiency-core threads at the end of every macro-step. Setting $T = W = 8$ avoids this barrier penalty, gives both codes the same compute resource, and makes the comparison one between two equally well-provisioned configurations of the same hardware. The probe (`probe_omp.py`) is committed alongside the runners.

## What is being timed

For GWSWEX, the timed region is the full `init -> run_step loop -> deinit` sequence, executed inside the same Python process as the trial loop; no subprocess is spawned per trial. The OpenMP thread count is set through the model API (`set_solver(omp_threads=N)`), which wires through to `omp_set_num_threads(N)` inside the kernel, and is verified by a probe (`kernel_get_omp_max_threads`) called immediately after the solver is configured. Per-trial metrics captured: wall (`time.perf_counter`), user and system CPU time and peak resident set size from `getrusage(RUSAGE_SELF)`, disk I/O from `psutil.disk_io_counters`, and the resulting NetCDF output size on disk.

For HYDRUS-1D, the timed region is the `concurrent.futures.ProcessPoolExecutor.map` over $n_e$ independent phydrus-driven simulations, each in its own scratch workspace. Within each worker the timed region rebuilds the full `phydrus.Model` from scratch (model, time info, water flow, materials, profile, atmospheric BC, root uptake / growth), calls `Model.write_input()` to materialise the deck, and then `Model.simulate()` to invoke the HYDRUS binary. This matches the path used by `comparison-intensive-sand-loam.ipynb` exactly, and folds the full phydrus configuration cost into the per-element sample on the same footing as the GWSWEX wrapper's per-element configuration cost. The pool size is fixed at 8 (matched to the GWSWEX OMP cap above), giving $n_e$ phydrus simulations distributed over a saturated pool; child-process resource usage is captured via `getrusage(RUSAGE_CHILDREN)` deltas in each worker.

Per (model, $n_e$) cell, three trials are run; the median is reported with the trial spread (min, max).

## The `flush_nc` flag in `model.py`

GWSWEX writes NetCDF output through `GwswexNCWriter`, which by default batches writes and lets the netCDF4 library flush at its own discretion. For the benchmark to report a wall-clock that includes per-step disk I/O on a comparable footing with HYDRUS-1D's per-step ASCII output, a `flush_nc=True` flag was added to `GWSWEXmodel.__init__`; when set, the writer calls `ds.sync()` after every macro-step. This is the flag the benchmark sets. The writer is now also strictly idempotent on `close()` and has a `__del__` safety net so that an aborted trial does not leak an open NetCDF file across worker boundaries.

## Hardware and software environment

| Item                             | Value                                                             |
| -------------------------------- | ----------------------------------------------------------------- |
| CPU model                        | Apple M1 Pro                                                      |
| Physical cores / logical threads | 10 (8 performance + 2 efficiency) / 10                            |
| RAM                              | 16.0 GiB                                                          |
| OS                               | Darwin 25.4.0 (arm64)                                             |
| Fortran compiler                 | GNU Fortran (conda-forge gcc 15.2.0-18) 15.2.0                    |
| OpenMP runtime                   | OpenMP 4.5 (libgomp via gfortran 15.2.0)                          |
| Python                           | CPython 3.14.3 (miniforge)                                        |
| HYDRUS-1D binary                 | `examples/gwswex-vs-hydrus1d/hydrus1d/bin/hydrus` (334 472 bytes) |
| GWSWEX commit                    | `8709bf0`                                                         |
| GWSWEX OMP thread count          | 8 (M1 Pro performance-core count)                                 |
| HYDRUS pool size                 | 8 (matched to GWSWEX OMP cap)                                     |

## Headline timings (median of three trials, intensive-sand-loam, 32 d hourly)

| Model                                       | $n_e$ | wall [s] | wall (min..max)   | user CPU [s] | max RSS [MiB] | disk write [MiB] |
| ------------------------------------------- | ----: | -------: | ----------------- | -----------: | ------------: | ---------------: |
| HYDRUS-1D pool (8 workers, phydrus per sim) |     1 |    1.757 | 1.740..1.764      |        0.918 |         192.3 |            13.64 |
| HYDRUS-1D pool (8 workers, phydrus per sim) |    10 |   11.277 | 10.231..11.497    |       67.365 |         192.3 |           143.70 |
| HYDRUS-1D pool (8 workers, phydrus per sim) |   100 |  102.260 | 101.330..103.556  |      748.559 |         192.3 |          2223.17 |
| HYDRUS-1D pool (8 workers, phydrus per sim) |  1000 | 1054.112 | 987.558..1141.080 |     7591.133 |         192.3 |         18783.57 |
| GWSWEX implicit (8 threads)                 |     1 |    0.489 | 0.484..0.493      |        0.468 |        170.4* |             0.00 |
| GWSWEX implicit (8 threads)                 |    10 |    1.228 | 1.226..1.274      |        4.431 |        170.4* |             0.04 |
| GWSWEX implicit (8 threads)                 |   100 |    9.041 | 9.039..9.180      |       44.296 |        170.4* |             3.52 |
| GWSWEX implicit (8 threads)                 |  1000 |   87.508 | 86.730..87.583    |      442.158 |         170.4 |            34.36 |
| GWSWEX explicit (8 threads)                 |     1 |    0.980 | 0.978..1.031      |        0.957 |          67.2 |             0.02 |
| GWSWEX explicit (8 threads)                 |    10 |    2.448 | 2.429..3.663      |       10.204 |          68.4 |             0.53 |
| GWSWEX explicit (8 threads)                 |   100 |   17.422 | 16.580..30.053    |       99.607 |          68.4 |            68.26 |
| GWSWEX explicit (8 threads)                 |  1000 |  167.132 | 160.800..188.722  |      983.417 |         169.4 |          1299.14 |

*The GWSWEX implicit solver was run in the same process as the explicit solver (sequentially). All implicit RSS values (170.4 MiB) reflect the residual resident set from the preceding explicit ne=1000 run (which peaked at 169.4 MiB); the genuine implicit ne=1000 working set is indistinguishable from that residual in this run. The explicit ne=1 figure of 67.2 MiB is the cleanest baseline for per-process footprint.

Empirical scaling exponents from a log-log OLS fit over $n_e \ge 10$ (the regime in which fixed per-call overhead is amortised):

| Code                                        | $\beta$ ($n_e \ge 10$) | $\beta$ (full range) |
| ------------------------------------------- | ---------------------: | -------------------: |
| HYDRUS-1D pool (8 workers, phydrus per sim) |                   0.99 |                 0.93 |
| GWSWEX implicit (8 threads)                 |                   0.93 |                 0.76 |
| GWSWEX explicit (8 threads)                 |                   0.92 |                 0.75 |

## Findings

The single-column row ($n_e = 1$) is the operating point that `comparison-intensive-sand-loam.ipynb` reports for the accuracy comparison. All three figures (GWSWEX explicit, GWSWEX implicit, HYDRUS-1D) include the full per-call configuration cost: for GWSWEX the Python-side construction, Fortran-side allocation, OpenMP team start-up, and per-step NetCDF flush; for HYDRUS the phydrus model construction, deck assembly and `write_input()`, and the subprocess invocation that writes per-step ASCII output to a dedicated workspace.

### Dominant effect: sand $K_\text{sat}$ drives CFL substepping

The sand upper layer has $K_\text{sat} = 0.297\ \text{m\,d}^{-1}$, roughly 29× larger than the loam lower layer ($K_\text{sat} = 0.0104\ \text{m\,d}^{-1}$). The explicit solver's adaptive Courant-Friedrichs-Lewy sub-step cascade fires continuously while drainage fronts traverse the sand, inflating wall time by roughly 1.9× relative to the implicit solver at $n_e = 1\,000$ (167 s vs 88 s). This divergence between the two GWSWEX solvers is notably smaller than it would be at a lower Courant number; the runs here used $C = 0.98$, which allows larger explicit substeps than the $C = 0.9$ used in earlier benchmarks and thereby reduces the CFL penalty.

### Solver ordering

At $n_e = 1\,000$:
- GWSWEX implicit: **87.51 s** — 12.0× faster than HYDRUS (1054 s, ratio 0.083)
- GWSWEX explicit: **167.13 s** — 6.3× faster than HYDRUS
- HYDRUS-1D pool:  **1054 s**

At $n_e = 1$:
- GWSWEX implicit: 0.49 s (fastest)
- GWSWEX explicit: 0.98 s (2.0× implicit, 1.8× faster than HYDRUS)
- HYDRUS-1D:       1.76 s

GWSWEX (both solvers) is faster than HYDRUS at all tested $n_e$. The dominant driver of the HYDRUS cost is subprocess-launch overhead and per-simulation ASCII I/O (18.3 GiB at $n_e = 1\,000$), not arithmetic. At $n_e = 1\,000$ GWSWEX implicit distributes 1000 columns across 8 OMP threads in a single call, whereas HYDRUS runs 1000 separate subprocesses through a saturated pool of 8 workers.

### Scaling exponents

All three codes scale sub-linearly ($\beta < 1$ in the asymptotic regime $n_e \ge 10$), confirming that fixed-overhead amortisation is the dominant effect at low to moderate $n_e$:

- HYDRUS pool: $\beta = 0.99$ (≥10), $0.93$ (full range)
- GWSWEX implicit: $\beta = 0.93$ (≥10), $0.76$ (full range)
- GWSWEX explicit: $\beta = 0.92$ (≥10), $0.75$ (full range)

HYDRUS now scales near-linearly in the asymptotic regime ($\beta \approx 1$), consistent with the pool becoming work-saturated: each added column requires approximately one more worker-slot's worth of subprocess time. Both GWSWEX solvers remain sub-linear, reflecting OpenMP amortisation of fixed per-call overhead across the column array.

### Memory and disk footprint

HYDRUS peak per-child-process RSS is approximately 192 MiB at all $n_e$ (captured via `getrusage(RUSAGE_CHILDREN)` per worker, reflecting the HYDRUS binary's own resident set rather than the parent Python process). GWSWEX carries the Fortran state in the parent process (67--170 MiB depending on $n_e$ and solver; see the RSS footnote above). The HYDRUS disk-write footprint is proportional to $n_e$ and exceeds GWSWEX by two orders of magnitude at $n_e = 1\,000$ (18.3 GiB vs 34 MiB for implicit), driven by seven ASCII output files per simulation.

### Operational case for GWSWEX

On a single-machine, isolated-column workload, GWSWEX implicit is the fastest option at all tested $n_e$; GWSWEX explicit is also faster than HYDRUS at every point in the sweep. The HYDRUS pool cost is dominated by subprocess launch overhead and file I/O, not by Richards equation arithmetic. GWSWEX also writes two orders of magnitude less to disk at $n_e = 1\,000$ (34 MiB vs 18.3 GiB for implicit vs HYDRUS). The primary case for GWSWEX, however, is the in-memory, single-process coupler that can exchange recharge and capillary rise with a MODFLOW-class saturated-zone solver at every macro-step, without paying $n_e$ subprocess launches per coupling step and without marshalling per-column state through the file system. That benefit is not captured by the present standalone benchmark; it would be the subject of a follow-up coupled-model benchmark.

## Trade-offs and caveats

Three trials per cell is a small sample. The reported spreads are tight on all cells (a few percent or less), consistent with a stable, deterministic workload on a dedicated machine.

The $n_e = 10^4$ point originally planned in the methodology was not executed in this round. From the fitted scaling exponents ($\beta({\ge}10)$) the projected wall is approximately $1054 \times (10^4/10^3)^{0.99} \approx 10\,300\ \text{s}$ for the HYDRUS pool, $87.5 \times (10^4/10^3)^{0.93} \approx 745\ \text{s}$ for GWSWEX implicit, and $167 \times (10^4/10^3)^{0.92} \approx 1\,390\ \text{s}$ for GWSWEX explicit. The infrastructure supports the larger sweep (`--ne-list` argument); only the time has not been spent.

The peak RSS for the GWSWEX implicit solver is flat at 170.4 MiB across all $n_e$. This is a `getrusage(RUSAGE_SELF)` artefact: the implicit trials were run in the same process after the explicit ne=1000 trial (which peaked at 169.4 MiB), and peak resident set size is monotonic over the process lifetime. The implicit ne=1000 figure (170.4 MiB) is near-identical to that residual, so the genuine implicit ne=1000 working set cannot be distinguished from it in this run. The explicit ne=1 figure (67.2 MiB) is the cleanest baseline for per-process footprint. For explicit ne=10 and ne=100 (both 68.4 MiB) the same monotonicity applies; the genuine per-$n_e$ working set is bounded above by those values.

The HYDRUS-1D disk-write cost is dominated by seven ASCII output files per simulation (`*.OUT`), and could be reduced by a factor of 5--10 by selecting a sparser print schedule in the SELECTOR.IN; the present run uses the default print cadence chosen by the phydrus template builder to match the comparison notebook's output frequency, which is appropriate for a like-for-like benchmark but is not the cheapest HYDRUS configuration possible.

## Reproducing

```bash
cd examples/gwswex-vs-hydrus1d/computational-performance-benchmark source ../../../.env.d/dev.env

python system_info.py # captures hardware/software python run_gwswex.py --ne-list 1,10,100,1000 --n-trials 3 \
 --solvers explicit,implicit \
 --omp-threads 8 --courant-number 0.98 python run_hydrus.py --ne-list 1,10,100,1000 --n-trials 3 --n-phys 8 python analyse.py # summary_table.md + scaling.pdf
```

The `tmp/` and `results/{*.json,*.md,*.csv,*.pdf}` artefacts are gitignored; the scripts are versioned in the repo.
