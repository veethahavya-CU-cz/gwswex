# type: ignore
"""HYDRUS-1D benchmark: run n_e independent simulations of the intensive
column, distributing them across a ``multiprocessing.pool.ThreadPool``
of size `n_physical_cores`.

Each worker drives HYDRUS-1D end-to-end through ``phydrus``: it instantiates
a fresh ``ps.Model`` per simulation, registers materials / profile / forcing
/ root parameters, writes the input deck, calls ``ml.simulate()``, and then
calls ``ml.read_nod_inf()`` to parse the ASCII Nod_Inf.out output. This
matches the full path used by the comparison-intensive-sand-loam.ipynb
notebook exactly: the notebook also parses Nod_Inf.out (and T_LEVEL.OUT)
after simulate(), so the per-element wall-clock here is directly comparable
to a real phydrus user's experience.

Solver settings mirror the notebook: ``dtmax=1/48`` d (30-min cap) and
``maxit=500`` Picard iterations, matching ``add_time_info`` and
``add_waterflow`` in the notebook.

Uses ``multiprocessing.pool.ThreadPool`` so that workers share the same
process address space — matching the real-world use case where multiple
phydrus.Model instances share Python-side memory and file-system I/O state.
The HYDRUS-1D binary itself runs as a subprocess in each thread, so threads
run concurrently (the GIL is released during subprocess.wait). cfg is passed
directly in the task tuple without pickling.

Per (ne) cell: `n_trials` repetitions; each trial captures
  - total wall time (pool wall)
  - aggregated per-sim user/sys CPU time (sum of subprocess RUSAGE_CHILDREN)
  - peak resident set size across all HYDRUS processes (psutil sampling)
  - aggregated disk I/O delta
"""
from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
import pandas as pd
import psutil
from benchmark_common import (HYDRUS_EXE, N_TRIALS_DEFAULT, NE_DEFAULT,
                              RESULTS_DIR, TMP_DIR, derive_intensive_case,
                              ensure_dirs)

WORKSPACES = TMP_DIR / "hydrus_workers"


def _build_model_in(ws_dir: Path, cfg: dict, sim_id: int):
    """Build a fresh phydrus.Model in `ws_dir` from scratch using cfg.

    This is the per-simulation setup path used by both the one-shot
    template build (smoke test) and the per-worker runs. Calling it
    inside each worker is the intent: every benchmark sample includes
    the full Model() / add_* / write_input cost, exactly as a phydrus
    user experiences it (and exactly as comparison-intensive-loam.ipynb
    does).
    """
    import phydrus as ps

    Z_TOP = cfg["Z_TOP"]
    NL = cfg["NL"]
    Z_WT = cfg["Z_WT"]
    T_TOTAL = cfg["T_TOTAL_D"]

    Z_COL_CM = Z_TOP * 100.0
    Z_WT_CM = (Z_TOP - Z_WT) * 100.0
    bot_cm = [b * 100.0 for b in cfg["bot_layer_elevations"]]
    n_mats = len(cfg["materials"])

    from benchmark_common import SOIL_TAG
    ml = ps.Model(
        exe_name=str(HYDRUS_EXE),
        ws_name=str(ws_dir),
        name=f"perf-bench-{SOIL_TAG}-{sim_id}",
        description=f"intensive-{SOIL_TAG} column, perf benchmark",
        mass_units="-",
        time_unit="days",
        length_unit="cm",
        print_screen=False,
    )
    ml.add_time_info(tinit=0, tmax=T_TOTAL, print_times=True,
                     dt=0.001, dtmin=1e-6, dtmax=1 / 48, dtprint=1 / 24)
    ml.add_waterflow(top_bc=2, bot_bc=1, rbot=0.0, rroot=0.0,
                     maxit=500, tolh=1.0, tolth=1e-3, ha=1e-6, hb=1e4)

    mat_df = ml.get_empty_material_df(n=n_mats)
    for m in cfg["materials"]:
        v = m["vanG"]
        K_cmd = m["K_sat"] * 100.0 * 24.0  # cfg stores m h-1 (intensive); HYDRUS wants cm d-1
        mat_df.loc[m["id"]] = [v["theta_r"], v["theta_s"],
                                v["alpha"] / 100.0, v["n"],
                                K_cmd, m["lam"]]
    ml.add_material(mat_df)

    profile = ps.create_profile(top=0.0, bot=-Z_COL_CM, dx=1.0, mat=1)
    z_nodes = profile["x"].to_numpy(float)
    layer_top_depth = [0.0] + [(Z_TOP * 100.0 - bot_cm[i]) for i in range(NL)]
    mat_per_node = []
    for zn in z_nodes:
        d_from_surface = -zn
        chosen = cfg["layer_mat_ids"][-1]
        for i in range(NL):
            if d_from_surface <= layer_top_depth[i + 1] + 1e-9:
                chosen = cfg["layer_mat_ids"][i]
                break
        mat_per_node.append(int(chosen))
    profile["Mat"] = mat_per_node
    profile["h"] = -Z_WT_CM - profile["x"].to_numpy(float)
    z = profile["x"].to_numpy(float)
    rd_cm = cfg["ROOT_D0"] * 100.0
    beta = np.where(z >= -rd_cm, np.maximum(1.0 + z / rd_cm, 0.0), 0.0)
    dz_cm = abs(float(profile["x"].iloc[1] - profile["x"].iloc[0]))
    integ = float(np.trapezoid(beta, dx=dz_cm)) if hasattr(np, "trapezoid") \
        else float(np.trapz(beta, dx=dz_cm))
    if integ > 0:
        beta = beta / integ
    profile["Beta"] = beta
    ml.add_profile(profile)
    ml.add_obs_nodes([-25.0, -50.0, -75.0, -100.0, -120.0, -140.0])

    atm = pd.DataFrame({
        "tAtm": cfg["day_idx"],
        "Prec": cfg["prec_d"],
        "rSoil": cfg["pet_d"],
        "rRoot": cfg["ptt_d"],
        "hCritA": np.full_like(cfg["day_idx"], 1e5),
        "rB": np.zeros_like(cfg["day_idx"]),
        "hB": np.zeros_like(cfg["day_idx"]),
        "ht": np.zeros_like(cfg["day_idx"]),
    })
    ml.add_atmospheric_bc(atm)
    ml.add_root_uptake(model=0, poptm=[-25] * n_mats, p0=-10, p2h=-200,
                       p2l=-800, p3=-8000)
    ml.add_root_growth(irootin=2, irfak=1, trmin=0, trmed=0, trmax=T_TOTAL,
                       xrmin=cfg["ROOT_D0"] * 100.0, xrmed=0,
                       xrmax=cfg["ROOT_D1"] * 100.0, trperiod=365)
    ml.write_input()
    return ml


def _smoke_build(template_dir: Path, cfg: dict) -> None:
    """One-shot smoke test: build a model in `template_dir` so the user can
    inspect a representative HYDRUS workspace before the parallel sweep
    starts. Not used by workers."""
    if template_dir.exists():
        shutil.rmtree(template_dir)
    template_dir.mkdir(parents=True)
    _build_model_in(template_dir, cfg, sim_id=-1)


# ---- worker --------------------------------------------------------------

def _hydrus_worker(args: tuple) -> dict:
    """One full phydrus-driven HYDRUS-1D simulation in its own scratch
    directory.

    Each call:
      1. mkdir's a fresh per-sim workspace,
      2. constructs a fresh ``phydrus.Model`` (materials, profile,
         atmospheric BC, root parameters),
      3. ``write_input()`` -> ``ml.simulate()`` -> ``ml.read_nod_inf()``.

    The wall-clock returned spans steps 1-3 and so includes the full
    setup + simulate + ASCII-parse cost a real phydrus user incurs per
    run, matching what the comparison notebook measures.

    cfg is passed directly in the tuple — no pickling overhead because
    workers are threads sharing the same process address space.
    """
    worker_id, sim_id, cfg = args
    ws = WORKSPACES / f"w{worker_id}_s{sim_id}"
    if ws.exists():
        shutil.rmtree(ws)
    ws.mkdir(parents=True)

    ru_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.perf_counter()
    try:
        ml = _build_model_in(ws, cfg, sim_id=sim_id)
        result = ml.simulate()
        rc = result.returncode
        if rc == 0:
            ml.read_nod_inf()  # match notebook: parse ASCII Nod_Inf.out inside timed region
    except Exception as exc:
        rc = 99
        result = None
        err_exc = repr(exc)
    else:
        err_exc = ""
    t1 = time.perf_counter()
    ru_after = resource.getrusage(resource.RUSAGE_CHILDREN)

    if rc != 0 or not (ws / "Nod_Inf.out").exists():
        err = (ws / "Error.msg").read_text() if (ws / "Error.msg").exists() else ""
        return dict(
            rc=(rc if rc != 0 else 2),
            wall_s=t1 - t0,
            user_cpu_s=ru_after.ru_utime - ru_before.ru_utime,
            sys_cpu_s=ru_after.ru_stime - ru_before.ru_stime,
            max_rss_kb=int(ru_after.ru_maxrss),
            worker_id=worker_id, sim_id=sim_id,
            error_msg=(err + "\n" + err_exc)[-2000:],
        )

    # Best-effort cleanup of large outputs to keep tmp/ small
    try:
        shutil.rmtree(ws)
    except Exception:
        pass

    return dict(
        rc=rc,
        wall_s=t1 - t0,
        user_cpu_s=ru_after.ru_utime - ru_before.ru_utime,
        sys_cpu_s=ru_after.ru_stime - ru_before.ru_stime,
        max_rss_kb=int(ru_after.ru_maxrss),
        worker_id=worker_id, sim_id=sim_id,
    )


def _run_trial(cfg: dict, ne: int, n_phys: int, trial: int) -> dict:
    sim_args = [(i % n_phys, i, cfg) for i in range(ne)]
    try:
        io_before = psutil.disk_io_counters()
    except Exception:
        io_before = None
    t0 = time.perf_counter()
    with ThreadPool(processes=n_phys) as pool:
        per_sim = pool.map(_hydrus_worker, sim_args)
    wall = time.perf_counter() - t0
    try:
        io_after = psutil.disk_io_counters()
    except Exception:
        io_after = None

    n_ok = sum(1 for r in per_sim if r["rc"] == 0)
    if n_ok != ne:
        sample = next((r for r in per_sim if r["rc"] != 0), None)
        return {"error": f"only {n_ok}/{ne} HYDRUS sims succeeded",
                "sample_failure": sample}

    walls = [r["wall_s"] for r in per_sim]
    users = [r["user_cpu_s"] for r in per_sim]
    syss = [r["sys_cpu_s"] for r in per_sim]
    rss = [r["max_rss_kb"] for r in per_sim]

    metrics = dict(
        ne=ne, trial=trial, n_phys_workers=n_phys, n_sims=ne,
        wall_s=wall,
        per_sim_wall_mean_s=float(np.mean(walls)),
        per_sim_wall_min_s=float(np.min(walls)),
        per_sim_wall_max_s=float(np.max(walls)),
        per_sim_user_cpu_mean_s=float(np.mean(users)),
        per_sim_sys_cpu_mean_s=float(np.mean(syss)),
        total_user_cpu_s=float(np.sum(users)),
        total_sys_cpu_s=float(np.sum(syss)),
        per_sim_max_rss_kb_max=int(np.max(rss)),
        per_sim_max_rss_kb_mean=int(np.mean(rss)),
    )
    if io_before is not None and io_after is not None:
        metrics["disk_read_bytes"] = int(io_after.read_bytes - io_before.read_bytes)
        metrics["disk_write_bytes"] = int(io_after.write_bytes - io_before.write_bytes)
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    # Default: skip 10000 — the user can pass it explicitly when willing to
    # wait several hours. The script is functionally identical at that scale.
    parser.add_argument("--ne-list", type=str,
                        default=",".join(str(x) for x in NE_DEFAULT if x <= 1000))
    parser.add_argument("--n-trials", type=int, default=N_TRIALS_DEFAULT)
    parser.add_argument("--n-phys", type=int,
                        default=min(8, psutil.cpu_count(logical=False) or 1),
                        help="HYDRUS pool size (default: min(8, physical cores); "
                             "matches the GWSWEX OMP cap so both codes share "
                             "the same set of P-cores on Apple silicon)")
    args = parser.parse_args()

    if not HYDRUS_EXE.exists():
        sys.exit(f"HYDRUS-1D binary not found at {HYDRUS_EXE}")

    ensure_dirs()
    WORKSPACES.mkdir(parents=True, exist_ok=True)

    cfg = derive_intensive_case()
    template_dir = TMP_DIR / "hydrus_template"
    print(f"[hydrus] smoke-building one phydrus workspace at {template_dir} "
          f"(workers will rebuild from scratch per simulation)", flush=True)
    _smoke_build(template_dir, cfg)

    ne_list = [int(x) for x in args.ne_list.split(",") if x.strip()]
    out_path = RESULTS_DIR / "hydrus_results.json"
    payload = {"n_phys_workers": args.n_phys,
               "setup_path": "phydrus.Model rebuilt per simulation",
               "trials": []}

    for ne in ne_list:
        for trial in range(args.n_trials):
            print(f"[hydrus] ne={ne} trial={trial} ...", flush=True)
            m = _run_trial(cfg, ne, args.n_phys, trial)
            payload["trials"].append(m)
            out_path.write_text(json.dumps(payload, indent=2))
            if "error" in m:
                print(f"  ERROR: {m['error']}", flush=True)
            else:
                print(f"  wall={m['wall_s']:.2f}s "
                      f"per_sim_wall_mean={m['per_sim_wall_mean_s']:.3f}s "
                      f"rss_max={m['per_sim_max_rss_kb_max']} kB", flush=True)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
