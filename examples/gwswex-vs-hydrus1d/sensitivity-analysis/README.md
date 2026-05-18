# GWSWEX OAT sensitivity analysis — reproducibility kit

This directory contains the standalone scripts used to drive the
one-at-a-time (OAT) parameter sensitivity analysis that backs the results
in `docs/sensitivity-analysis-from-comparison.md` and the tuned
configurations baked into the twelve `examples/gwswex-vs-hydrus1d/`
comparison notebooks.

The scripts are independent of the GWSWEX test suite (they are not
collected by `pytest`); they are placed under `tests/` so they live with
the rest of the verification material and can be re-run by anyone with a
working GWSWEX + HYDRUS-1D environment.

## Contents

| File | Role |
|---|---|
| `oat_harness.py` | Orchestrator. Iterated coordinate-descent OAT across all 24 (soil $\times$ setup $\times$ solver) cells. Spawns one `oat_worker.py` subprocess per trial. Writes `oat_results/oat_results.json`. |
| `oat_worker.py` | Single-trial GWSWEX runner. Reads a JSON parameter set from `argv[1]`, executes the relevant comparison-notebook cells with those parameters substituted in, prints `RMSE_CM <value>` (or `ERROR <msg>`) on stdout. Subprocess-isolated to protect the harness from Fortran-kernel state corruption when an individual trial crashes. |
| `apply_oat_optima.py` | Reads `oat_results/oat_results.json` and produces `oat_results/tuned_params.json` — one accepted configuration per (soil, setup, solver) cell. Per-setup acceptance margins: basic $\geq 2$ %, intensive $\geq 5$ %; cells that do not clear the margin revert to the notebook baseline. |
| `gen_variants.py` | Template-substitution notebook generator. Reads the tuned-params JSON and (re)writes `set_solver(...)`, `set_model_params(...)`, and `ET_STRESS = dict(...)` lines in each comparison notebook. |
| `fill_sensitivity_tables.py` | Fills the per-case parameter tables in `docs/sensitivity-analysis-from-comparison.md` §2.3 from `tuned_params.json`. |
| `validate_all.py` | End-to-end validation. Re-executes the twelve comparison notebooks via `jupyter nbconvert --to notebook --execute --inplace`, scrapes the printed Full-period WT RMSE / NSE per solver, and writes `oat_results/validation_summary.md`. |
| `oat_results.reference.json` | Reference OAT output produced on 2026-04-22 against commit `f4b1c6b`-era notebooks. Used to sanity-check that a fresh re-run reproduces the same picks. |
| `tuned_params.reference.json` | Reference tuned-params output corresponding to `oat_results.reference.json`. |

All scripts also have a header docstring describing their inputs and outputs.

## Prerequisites

- A working GWSWEX install (`pip install -e .` from the repo root, with
  the project's `meson` + Fortran toolchain).
- HYDRUS-1D 4.x (from the vendored binary at
  `examples/gwswex-vs-hydrus1d/hydrus1d/bin/hydrus`) and `phydrus`.
- `jupyter`, `nbconvert`, `nbformat`, `numpy`, `xarray`, `matplotlib`.
- A POSIX shell. The OAT runs are CPU-only and shell out via `subprocess`.

The recommended environment is the project's `gwswex` conda env (see
`.env.d/dev.env`).

## Reproducing the full study

From the repository root:

```bash
# 0. Activate the GWSWEX environment
source .env.d/dev.env

# 1. Run the full OAT sweep over all 24 cases (~3-4 h on a recent laptop).
#    On macOS with sandboxed shells, set TMPDIR to a writeable location,
#    a non-default MPLCONFIGDIR, and disable KMP affinity.
export TMPDIR="$HOME/.tmp"          # any writeable dir is fine
export MPLCONFIGDIR="$TMPDIR/mpl"
export KMP_AFFINITY=disabled
mkdir -p "$MPLCONFIGDIR"
$PYTHON tests/sensitivity-analysis/oat_harness.py --all

# 2. Promote winning trials into a tuned-params JSON, applying the
#    per-setup acceptance margins.
$PYTHON tests/sensitivity-analysis/apply_oat_optima.py

# 3. Bake the tuned values into the comparison notebooks.
$PYTHON tests/sensitivity-analysis/gen_variants.py

# 4. Update the per-case parameter tables in the sensitivity doc.
$PYTHON tests/sensitivity-analysis/fill_sensitivity_tables.py

# 5. Re-execute every comparison notebook end-to-end and record the
#    full-pipeline WT RMSE/NSE.
$PYTHON tests/sensitivity-analysis/validate_all.py
```

After step 1, `tests/sensitivity-analysis/oat_results/oat_results.json`
should be byte-identical (modulo timestamp keys, if any) to
`oat_results.reference.json`. Small differences in the trailing decimals
of RMSEs are expected when the underlying GWSWEX or HYDRUS-1D versions
have changed; the picks themselves should remain stable as long as
neither the model nor the OAT grid has been modified.

## Re-running a single case for debugging

The harness accepts `--cases`, `--setups`, and `--solvers` filters:

```bash
$PYTHON tests/sensitivity-analysis/oat_harness.py \
    --cases loam --setups basic --solvers implicit
```

This is a useful smoke test (about 10 minutes on a recent laptop, of
which roughly one minute is the HYDRUS-1D reference run; the rest is the
GWSWEX trials).

## What the OAT does, briefly

For every (soil, setup, solver) cell, the harness performs **iterated
coordinate-descent OAT**:

1. Build a baseline by parsing the current notebook source: the
   `MODEL_PARAMS = dict(...)`, `ET_STRESS = dict(...)`, and the
   `set_solver(...)` keyword arguments. This means hand-tuned values
   carried in the notebook (e.g. `ICratio_min = 0.42` in
   `comparison-basic-clay.ipynb`) become the starting point for the next
   OAT round and are only displaced by a strict improvement.
2. Run the HYDRUS-1D reference once and cache the WT trajectory.
3. Per pass: sweep one parameter at a time across its discrete grid,
   holding the others at the current best. Accept the per-parameter
   winner if it reduces the WT-RMSE by at least the per-setup margin
   (basic: 2 %, intensive: 5 %) relative to the running RMSE. Continue
   into the next parameter from the just-accepted configuration.
4. Stop when no parameter improves the RMSE by at least 0.5 % within a
   pass, or when `MAX_PASSES = 3` is reached.

The acceptance margins, sweep grids, and stopping criterion live at the
top of `oat_harness.py` and can be adjusted. The full sweep grid is
documented in `docs/sensitivity-analysis-from-comparison.md` §1.3.

## Subprocess isolation

`oat_harness.py` deliberately runs every trial in a fresh `oat_worker.py`
subprocess. This is to work around a known issue: GWSWEX uses a
module-level Fortran `Model` variable, and a failed trial leaves it in a
half-initialised state that causes the next allocate to fail with

```
At line 44 of file gwswex/src/kernel.f08
Fortran runtime error: Attempting to allocate already allocated variable 'model'
```

A subprocess per trial is the simplest robust fix; the per-trial cost is
roughly 4 seconds of Python startup, which is negligible compared to the
GWSWEX run itself.

## Output layout

```
tests/sensitivity-analysis/oat_results/
  cache/
    h_gw_d-{soil}-{setup}.npy        # cached HYDRUS WT references
  oat_results.json                   # full per-trial record
  tuned_params.json                  # per-cell accepted configuration
  validation_summary.md              # post-execution RMSE/NSE table
  validation_logs/
    val-{setup}-{soil}.log           # nbconvert stdout per notebook
  tmp/                               # ephemeral worker outputs
```
