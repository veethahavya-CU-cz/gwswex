# type: ignore
"""Aggregate per-trial metrics into median tables and log-log scaling figures.

Reads `results/{system_info,gwswex_results,hydrus_results}.json`. Writes:
  - `results/summary_table.csv`   : long-form per-(model, ne) median + min/max
  - `results/summary_table.md`    : same content, GitHub-flavoured markdown
  - `results/scaling.pdf`         : wall-clock vs n_e on log-log axes
  - `results/scaling_breakdown.pdf`: wall vs CPU-time stack per (model, ne)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from benchmark_common import RESULTS_DIR

# `getrusage` returns ru_maxrss in kilobytes on Linux, but in bytes on
# macOS / *BSD. Convert to MiB accordingly.
if sys.platform == "darwin":
    _RSS_TO_MIB = 1.0 / (1024 * 1024)
else:
    _RSS_TO_MIB = 1.0 / 1024

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "pdf.fonttype": 42,
})


def _summarise(trials: list[dict], group_keys: list[str], value_keys: list[str]) -> dict:
    """Group trials by the tuple (group_keys...) and return median + min + max
    for each value_key."""
    groups: dict = {}
    for t in trials:
        if "error" in t:
            continue
        k = tuple(t[g] for g in group_keys)
        groups.setdefault(k, []).append(t)
    out = {}
    for k, ts in groups.items():
        d = {}
        for v in value_keys:
            vals = [tt[v] for tt in ts if v in tt]
            if vals:
                d[v + "_med"] = float(np.median(vals))
                d[v + "_min"] = float(np.min(vals))
                d[v + "_max"] = float(np.max(vals))
        out[k] = d
    return out


def main() -> None:
    gw_default_path = RESULTS_DIR / "gwswex_results.json"
    if not gw_default_path.exists():
        gw_default_path = RESULTS_DIR / "gwswex_results_default.json"
    h_path = RESULTS_DIR / "hydrus_results.json"
    sys_path = RESULTS_DIR / "system_info.json"

    def _load(p):
        return json.loads(p.read_text()) if p.exists() else {"trials": []}

    gw_def = _load(gw_default_path)
    hy = _load(h_path)
    si = json.loads(sys_path.read_text()) if sys_path.exists() else {}

    soil_tag = gw_def.get("soil_tag", "unknown")
    all_trials = list(gw_def.get("trials", []))

    gw_keys = ["wall_s", "user_cpu_s", "sys_cpu_s", "max_rss_kb_self",
               "nc_bytes", "disk_write_bytes"]
    hy_keys = ["wall_s", "total_user_cpu_s", "total_sys_cpu_s",
               "per_sim_wall_mean_s", "per_sim_max_rss_kb_max",
               "disk_write_bytes"]

    gw_sum = _summarise(all_trials, ["solver", "ne"], gw_keys)
    hy_sum = _summarise(hy["trials"], ["ne"], hy_keys)

    # ---- table ----
    rows = []
    for (solver, ne), v in sorted(gw_sum.items()):
        rows.append({
            "model": f"GWSWEX-{solver}", "ne": ne,
            "wall_s_med": v.get("wall_s_med"),
            "wall_s_min": v.get("wall_s_min"),
            "wall_s_max": v.get("wall_s_max"),
            "user_cpu_s_med": v.get("user_cpu_s_med"),
            "sys_cpu_s_med": v.get("sys_cpu_s_med"),
            "max_rss_mib_med": v.get("max_rss_kb_self_med", 0) * _RSS_TO_MIB,
            "disk_write_mib_med": v.get("disk_write_bytes_med", 0) / 2**20,
            "nc_mib_med": v.get("nc_bytes_med", 0) / 2**20,
        })
    for (ne,), v in sorted(hy_sum.items()):
        rows.append({
            "model": "HYDRUS-1D", "ne": ne,
            "wall_s_med": v.get("wall_s_med"),
            "wall_s_min": v.get("wall_s_min"),
            "wall_s_max": v.get("wall_s_max"),
            "user_cpu_s_med": v.get("total_user_cpu_s_med"),
            "sys_cpu_s_med": v.get("total_sys_cpu_s_med"),
            "max_rss_mib_med": v.get("per_sim_max_rss_kb_max_med", 0) * _RSS_TO_MIB,
            "disk_write_mib_med": v.get("disk_write_bytes_med", 0) / 2**20,
            "nc_mib_med": None,
        })

    csv_path = RESULTS_DIR / "summary_table.csv"
    cols = list(rows[0].keys()) if rows else []
    with csv_path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join("" if r[c] is None else f"{r[c]}" for c in cols) + "\n")

    md_path = RESULTS_DIR / "summary_table.md"
    with md_path.open("w") as f:
        f.write(f"# Computational performance \u2014 intensive-{soil_tag}\n\n")
        if si:
            cpu = si.get("cpu", {})
            mem = si.get("memory", {})
            os_i = si.get("os", {})
            f.write(f"- CPU: {cpu.get('model', '?')} "
                    f"({cpu.get('physical_cores')}p / {cpu.get('logical_cores')}l)\n")
            f.write(f"- RAM: {mem.get('total_gib', '?')} GiB\n")
            f.write(f"- OS: {os_i.get('system', '?')} {os_i.get('release', '')} "
                    f"({os_i.get('arch', '')})\n")
            f.write(f"- OMP threads (GWSWEX): {gw_def.get('omp_threads')}\n")
            f.write(f"- HYDRUS pool size (workers): {hy.get('n_phys_workers')}\n")
            f.write(f"- HYDRUS setup path: {hy.get('setup_path', 'phydrus.Model rebuilt per simulation')}\n\n")
        f.write("| Model | n_e | wall [s] (med) | wall [s] (min..max) | user CPU [s] | sys CPU [s] | max RSS [MiB] | disk write [MiB] |\n")
        f.write("|---|---:|---:|---|---:|---:|---:|---:|\n")
        for r in rows:
            wm = r["wall_s_med"]
            wlo = r["wall_s_min"]
            whi = r["wall_s_max"]
            f.write(f"| {r['model']} | {r['ne']} | "
                    f"{wm:.3f} | {wlo:.3f}..{whi:.3f} | "
                    f"{r['user_cpu_s_med']:.3f} | {r['sys_cpu_s_med']:.3f} | "
                    f"{r['max_rss_mib_med']:.1f} | "
                    f"{r['disk_write_mib_med']:.2f} |\n")

    # ---- scaling figure ----
    fig, ax = plt.subplots(figsize=(5.4, 4.0), constrained_layout=True)
    series = {}
    for (solver, ne), v in gw_sum.items():
        series.setdefault(f"GWSWEX-{solver}", []).append((ne, v["wall_s_med"], v["wall_s_min"], v["wall_s_max"]))
    for (ne,), v in hy_sum.items():
        series.setdefault("HYDRUS-1D", []).append((ne, v["wall_s_med"], v["wall_s_min"], v["wall_s_max"]))
    colours = {"GWSWEX-explicit": "#d62728",
               "GWSWEX-implicit": "#1f77b4",
               "HYDRUS-1D": "#2ca02c"}
    markers = {"GWSWEX-explicit": "s",
               "GWSWEX-implicit": "o",
               "HYDRUS-1D": "^"}
    for name, pts in series.items():
        pts.sort()
        ne_arr = np.array([p[0] for p in pts])
        med = np.array([p[1] for p in pts])
        lo = np.array([p[2] for p in pts])
        hi = np.array([p[3] for p in pts])
        ax.fill_between(ne_arr, lo, hi, color=colours.get(name, "grey"), alpha=0.15)
        ax.plot(ne_arr, med, marker=markers.get(name, "x"),
                color=colours.get(name, "grey"), label=name, lw=1.4, ms=5)
        # log-log scaling exponent over asymptotic regime (ne >= 10 only)
        mask = ne_arr >= 10
        if mask.sum() >= 2 and np.all(med[mask] > 0):
            beta, logA = np.polyfit(np.log(ne_arr[mask]), np.log(med[mask]), 1)
            ax.plot([], [], ' ', label=fr"$\beta_{{{name}}}={beta:.2f}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"number of elements $n_e$")
    ax.set_ylabel("wall-clock time [s]")
    ax.set_title(f"Intensive-{soil_tag}, 32 d hourly, deferred NetCDF (track=False)")
    ax.legend(fontsize=8, loc="upper left")
    fig.savefig(RESULTS_DIR / "scaling.pdf")

    _ms_figs = Path(__file__).parents[3] / "docs" / "manuscript" / "figures"
    _ms_figs.mkdir(parents=True, exist_ok=True)
    fig.savefig(_ms_figs / "fig_perf_scaling.pdf")
    plt.close(fig)

    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"wrote {RESULTS_DIR / 'scaling.pdf'}")
    print(f"wrote {_ms_figs / 'fig_perf_scaling.pdf'}")


if __name__ == "__main__":
    main()
