# type: ignore
"""GWSWEX benchmark: sweep n_e and solver, capture per-trial metrics.

Each (solver, n_e) cell is run `n_trials` times in this same process. Per
trial we capture:

  - wall      : `time.perf_counter()` around `init -> run_step loop -> deinit`
  - user_cpu  : `resource.getrusage(RUSAGE_SELF).ru_utime`
  - sys_cpu   : `resource.getrusage(RUSAGE_SELF).ru_stime`
  - max_rss   : `ru_maxrss` (kB on Linux, bytes on macOS)
  - disk_read : `psutil.disk_io_counters()` delta
  - disk_write: `psutil.disk_io_counters()` delta
  - n_steps   : number of macro-steps actually executed (sanity check)

NetCDF output is written with deferred flushing (flush_nc=False, the
default), matching the behaviour of the comparison notebook
(`comparison-intensive-sand-loam.ipynb`). Per-step state extraction
(get_state, get_mass_balance) and NetCDF per-step writes are suppressed
via run_step(track=False) — the output file is still created and its
configuration metadata written by init()/deinit(), so the file exists for
post-run inspection. Per-step flushing (flush_nc=True) would add ~768
synchronous disk-sync calls and inflate wall time by ~1.5–2.5 s at ne=1,
which is not representative of normal usage.

The `--courant-number` default is 0.98 (maximum-efficiency explicit stepping).
Note that the comparison notebook uses courant_number=0.3 for accuracy fidelity;
the relaxation here is intentional and disclosed: this is a throughput benchmark,
not an accuracy benchmark. The explicit solver is run with ``n_trapz=5`` to match
the trapezoidal quadrature setting used in the comparison notebook.
"""

from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np
import psutil
from benchmark_common import N_TRIALS_DEFAULT, NE_DEFAULT, RESULTS_DIR, TMP_DIR, derive_intensive_case, ensure_dirs

# Default OMP thread count: number of P-cores on Apple silicon (8 on M1 Pro).
# `psutil.cpu_count(logical=False)` returns 10 on M1 Pro, which includes the
# 2 E-cores; OMP threads landing on E-cores stall the per-step barrier and
# regress total throughput, so we cap at 8 by default. Override with
# `--omp-threads`.
_PHYS = psutil.cpu_count(logical=False) or 1
_DEFAULT_OMP = min(8, _PHYS)


def _run_one_trial(solver: str, ne: int, trial: int, omp_threads: int, courant_number: float, cfg: dict) -> dict:
    from gwswex import GWSWEXmodel

    out_nc = TMP_DIR / f"gwswex_{solver}_ne{ne}_t{trial}.nc"
    if out_nc.exists():
        out_nc.unlink()

    NL = cfg["NL"]
    bot = [cfg["bot_layer_elevations"]] * ne
    top = [[cfg["Z_TOP"]] for _ in range(ne)]
    sID = [cfg["layer_mat_ids"] for _ in range(ne)]
    vID = [[1] for _ in range(ne)]

    m = GWSWEXmodel(name=f"perf-{solver}-ne{ne}-t{trial}", T="h", L="m", output_fpath=str(out_nc))
    m.init_space(ne=ne, nl=NL, top=top, bot=bot, sID=sID, vID=vID)
    for mat in cfg["materials"]:
        m.add_material(**mat)
    m.add_vegetation(
        id=1,
        name="bench-veg",
        root_depth_initial=cfg["ROOT_D0"],
        root_depth_final=cfg["ROOT_D1"],
        root_growth_model="linear",
        et_stress=dict(s_star=0.5, s_w=0.1, s_h=0.05, s_e=0.5),
    )
    m.init_time(n_steps=cfg["N"], dt=1.0, dt_min=1 / 60)
    m.set_model_params(psi_f=0.09, F_min=1e-06, ICratio_min=0.05)
    solver_kwargs: dict = dict(solver=solver, omp_threads=omp_threads, courant_number=courant_number)
    if solver == "explicit":
        solver_kwargs["n_trapz"] = 5  # match comparison-notebook setting
    m.set_solver(**solver_kwargs)
    m.set_initial_conditions(gw=cfg["Z_WT"], sw=0.0, uz=-999)

    prec = np.broadcast_to(cfg["prec_h"][:, None], (cfg["N"], ne)).copy()
    pet = np.broadcast_to(cfg["pet_h"][:, None], (cfg["N"], ne)).copy()
    ptt = np.broadcast_to(cfg["ptt_h"][:, None], (cfg["N"], ne)).copy()
    m.set_forcing(precip=prec, pet=pet, ptt=ptt)

    try:
        io_before = psutil.disk_io_counters()
    except Exception:
        io_before = None
    ru_before = resource.getrusage(resource.RUSAGE_SELF)
    t0 = time.perf_counter()

    m.init()
    for t in range(cfg["N"]):
        m.run_step(t, track=False)
    m.deinit()

    t1 = time.perf_counter()
    ru_after = resource.getrusage(resource.RUSAGE_SELF)
    try:
        io_after = psutil.disk_io_counters()
    except Exception:
        io_after = None

    nc_size = out_nc.stat().st_size if out_nc.exists() else 0
    if out_nc.exists():
        try:
            out_nc.unlink()
        except Exception:
            pass

    metrics = dict(
        solver=solver,
        ne=ne,
        trial=trial,
        n_steps=cfg["N"],
        wall_s=t1 - t0,
        user_cpu_s=ru_after.ru_utime - ru_before.ru_utime,
        sys_cpu_s=ru_after.ru_stime - ru_before.ru_stime,
        max_rss_kb_self=int(ru_after.ru_maxrss),
        omp_threads=int(omp_threads),
        courant_number=float(courant_number),
        nc_bytes=nc_size,
    )
    if io_before is not None and io_after is not None:
        metrics["disk_read_bytes"] = int(io_after.read_bytes - io_before.read_bytes)
        metrics["disk_write_bytes"] = int(io_after.write_bytes - io_before.write_bytes)
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ne-list", type=str, default=",".join(str(x) for x in NE_DEFAULT))
    parser.add_argument("--n-trials", type=int, default=N_TRIALS_DEFAULT)
    parser.add_argument("--solvers", type=str, default="explicit,implicit")
    parser.add_argument(
        "--omp-threads", type=int, default=_DEFAULT_OMP, help="OpenMP thread count (default: 8 = M1 Pro P-core count)"
    )
    parser.add_argument(
        "--courant-number",
        type=float,
        default=0.98,
        help="Courant safety factor for the explicit solver "
        "(default 0.98 for throughput benchmarking; note "
        "the comparison notebook uses 0.3 for accuracy)",
    )
    args = parser.parse_args()

    ensure_dirs()
    ne_list = [int(x) for x in args.ne_list.split(",") if x.strip()]
    solvers = [s.strip() for s in args.solvers.split(",") if s.strip()]

    cfg = derive_intensive_case()
    out_path = RESULTS_DIR / "gwswex_results.json"
    payload = {
        "omp_threads": args.omp_threads,
        "physical_cpus": _PHYS,
        "courant_number": args.courant_number,
        "soil_tag": __import__("benchmark_common").SOIL_TAG,
        "trials": [],
    }
    for solver in solvers:
        for ne in ne_list:
            for trial in range(args.n_trials):
                print(
                    f"[gwswex] solver={solver} ne={ne} trial={trial} "
                    f"omp_threads={args.omp_threads} cr={args.courant_number}",
                    flush=True,
                )
                t0 = time.perf_counter()
                try:
                    m = _run_one_trial(solver, ne, trial, args.omp_threads, args.courant_number, cfg)
                except Exception as e:
                    import traceback

                    traceback.print_exc()
                    m = {"error": f"{type(e).__name__}: {e}", "solver": solver, "ne": ne, "trial": trial}
                m["wall_outer_s"] = time.perf_counter() - t0
                payload["trials"].append(m)
                out_path.write_text(json.dumps(payload, indent=2))
                if "error" in m:
                    print(f"  ERROR: {m['error']}", flush=True)
                else:
                    print(
                        f"  wall={m['wall_s']:.3f}s "
                        f"user={m['user_cpu_s']:.3f}s "
                        f"rss={m['max_rss_kb_self']} kB "
                        f"nc_bytes={m['nc_bytes']}",
                        flush=True,
                    )
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
