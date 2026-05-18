# Computational performance benchmark

Benchmarks GWSWEX (explicit & implicit solvers) and HYDRUS-1D on the
**intensive-clay** profile from the `examples/gwswex-vs-hydrus1d/` source of truth
(`experiment_definitions.json`). The case is chosen because it is the regime
where the explicit operator-split solver is expected to be most strongly
ahead of both the implicit Picard solver and HYDRUS-1D: low-Ks clay under
hourly forcing produces a ponded surface that the explicit Green-Ampt
infiltration cap and bucket cascade resolve in O(1) per step, whereas the
mixed-form Richards solver in HYDRUS-1D and in the GWSWEX implicit kernel
both pay an iterative cost per step.

## Layout

```
benchmark_common.py    -- load SoT, build intensive-clay forcing arrays
system_info.py         -- snapshot of hardware + OS + Python + library versions
run_gwswex.py          -- sweep n_e in {1,10,100,1000,10000}, both solvers
run_hydrus.py          -- run n_e independent HYDRUS-1D sims via multiprocessing
analyse.py             -- summarise results into tables and log-log plots
results/               -- gitignored output (per-trial timings + system info)
tmp/                   -- gitignored scratch (per-element NetCDFs, HYDRUS workspaces)
```

## Method

- `OMP_NUM_THREADS` is set to the number of logical CPU cores for GWSWEX.
- HYDRUS-1D is single-threaded; `multiprocessing.Pool(n_physical_cores)` is used
  to run independent simulations in parallel, giving HYDRUS the same total
  hardware budget as GWSWEX.
- Each (model, solver, n_e) cell is repeated three times; median, min, max
  reported.
- Per-trial metrics: wall time, user CPU, system CPU, peak resident set size,
  cumulative disk read/write bytes, NetCDF write time (where applicable).
- GWSWEX uses `flush_nc=True` so that every macro-step writes through to disk,
  matching the HYDRUS-1D default of writing T-Level / NOD_INF every output time.

## Usage

```sh
source .env.d/dev.env

# 1. record hardware/software environment
$PYTHON examples/gwswex-vs-hydrus1d/computational-performance-benchmark/system_info.py

# 2. GWSWEX benchmark (explicit + implicit)
$PYTHON examples/gwswex-vs-hydrus1d/computational-performance-benchmark/run_gwswex.py

# 3. HYDRUS-1D benchmark (multiprocessing pool)
$PYTHON examples/gwswex-vs-hydrus1d/computational-performance-benchmark/run_hydrus.py

# 4. summary tables + figures
$PYTHON examples/gwswex-vs-hydrus1d/computational-performance-benchmark/analyse.py
```

Both `run_*.py` scripts accept `--ne-list` (comma-separated) and `--n-trials`.
