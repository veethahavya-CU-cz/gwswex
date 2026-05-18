"""Snapshot the hardware and software environment for the benchmark.

Writes `results/system_info.json`. Captures CPU model, physical/logical core
count, total RAM, OS, architecture, Python version, and the versions of every
library that materially affects benchmark wall time.
"""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import psutil

from benchmark_common import HYDRUS_EXE, RESULTS_DIR, ensure_dirs


def _safe(cmd: list[str]) -> str:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=10)
        return out.strip().splitlines()[0]
    except Exception:
        return ""


def _cpu_model() -> str:
    if sys.platform == "darwin":
        return _safe(["sysctl", "-n", "machdep.cpu.brand_string"])
    if sys.platform.startswith("linux"):
        try:
            for line in open("/proc/cpuinfo"):
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except Exception:
            pass
    return platform.processor() or platform.machine()


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, timeout=5
        ).strip()
    except Exception:
        return "unknown"


def _pkg_version(name: str) -> str:
    try:
        from importlib.metadata import version as _v
        return _v(name)
    except Exception:
        return "unknown"


def _hydrus_info() -> dict:
    info: dict = {"path": str(HYDRUS_EXE), "exists": HYDRUS_EXE.exists()}
    if HYDRUS_EXE.exists():
        try:
            info["size_bytes"] = HYDRUS_EXE.stat().st_size
        except Exception:
            pass
        # HYDRUS-1D prints a banner with its version when invoked without args
        try:
            out = subprocess.run(
                [str(HYDRUS_EXE)], capture_output=True, text=True, timeout=2,
                input="\n",
            )
            banner = (out.stdout + out.stderr).strip().splitlines()[:6]
            info["banner"] = banner
        except Exception:
            pass
    return info


def collect() -> dict:
    info = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": platform.node(),
        "os": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "arch": platform.architecture()[0],
        },
        "cpu": {
            "model": _cpu_model(),
            "physical_cores": psutil.cpu_count(logical=False),
            "logical_cores": psutil.cpu_count(logical=True),
            "max_freq_mhz": getattr(psutil.cpu_freq(), "max", None),
        },
        "memory": {
            "total_bytes": psutil.virtual_memory().total,
            "total_gib": round(psutil.virtual_memory().total / 2**30, 2),
        },
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "env": {
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", ""),
            "OMP_PROC_BIND": os.environ.get("OMP_PROC_BIND", ""),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS", ""),
        },
        "packages": {
            name: _pkg_version(name) for name in
            ["numpy", "scipy", "pandas", "netCDF4", "matplotlib",
             "psutil", "phydrus", "gwswex"]
        },
        "gwswex_commit": _git_commit(),
        "hydrus": _hydrus_info(),
        "tools": {
            "gfortran": _safe(["gfortran", "--version"]),
            "ifort": _safe(["ifort", "--version"]),
            "ninja": _safe(["ninja", "--version"]),
            "meson": _safe(["meson", "--version"]),
        },
        "disk_free_bytes": shutil.disk_usage(str(RESULTS_DIR.parent)).free,
    }
    return info


def main() -> None:
    ensure_dirs()
    info = collect()
    out = RESULTS_DIR / "system_info.json"
    out.write_text(json.dumps(info, indent=2, default=str))
    print(f"wrote {out}")
    print(json.dumps(info, indent=2, default=str))


if __name__ == "__main__":
    main()
