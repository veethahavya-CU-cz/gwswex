"""
GWSWEX: Groundwater-Vadose-Surface-Water Exchange Model

A field-scale vadose zone model with adaptive sub-stepping,
mass-conserving numerics, and coupling to external SW and GW models.
"""

from .version import __version__


def _gwswex_build_prompt(err: Exception) -> bool:
    """Interactive fallback when the compiled Fortran extension is missing.

    Offers two install routes and runs them via subprocess.  Returns True if
    the user chose an option and pip succeeded (caller should then ask the
    user to restart Python); False if the user declined or any step failed.
    """
    import shutil
    import subprocess
    import sys
    import sysconfig
    from pathlib import Path

    _pkg = Path(__file__).parent  # gwswex/
    _repo = _pkg.parent  # must be the repo root for builds to work

    if not (_repo / "meson.build").exists():
        return False  # not a source checkout; can't help

    _py = sys.executable
    _so = f"f_gwswex{sysconfig.get_config_var('EXT_SUFFIX') or '.so'}"
    _searched = [str(Path(p) / "gwswex" / _so) for p in sys.path if p]

    print(f"\n[gwswex] ImportError: {err}", file=sys.stderr)
    print(f"  Python executable : {_py}", file=sys.stderr)
    print(f"  Extension name    : {_so}", file=sys.stderr)
    print(f"  Locations checked :", file=sys.stderr)
    for _s in _searched:
        print(f"    {_s}", file=sys.stderr)
    print(f"  Repository root   : {_repo}", file=sys.stderr)
    print(file=sys.stderr)
    print("Compile from source?  (repo root above must be the GWSWEX repo root)", file=sys.stderr)
    print("  1  Non-editable — compile and copy the extension to site-packages.", file=sys.stderr)
    print("     Once installed the package is independent of this directory.", file=sys.stderr)
    print("  2  Editable — compile here; package refers back to this repo.", file=sys.stderr)
    print("     Changes to Python source are immediately visible; Fortran changes", file=sys.stderr)
    print("     require a recompile.  Suitable for development.", file=sys.stderr)
    print("  q  Quit — raise ImportError without building.", file=sys.stderr)

    try:
        _choice = input("\nChoice [1/2/q]: ").strip().lower()
    except (EOFError, OSError):
        print("[gwswex] Non-interactive environment — cannot prompt.", file=sys.stderr)
        return False

    if _choice not in ("1", "2"):
        return False

    _meson = shutil.which("meson")
    if _meson is None:
        print(
            f"[gwswex] 'meson' not found on PATH.  Install it first:\n" f"  {_py} -m pip install meson",
            file=sys.stderr,
        )
        return False

    _build = _repo / "build"
    if not _build.exists():
        print("[gwswex] Running: meson setup build ...", file=sys.stderr)
        if subprocess.run([_meson, "setup", "build"], cwd=_repo).returncode != 0:
            print("[gwswex] meson setup failed.", file=sys.stderr)
            return False

    print("[gwswex] Running: meson compile -C build ...", file=sys.stderr)
    if subprocess.run([_meson, "compile", "-C", "build"], cwd=_repo).returncode != 0:
        print("[gwswex] meson compile failed.", file=sys.stderr)
        return False

    _cmd = [_py, "-m", "pip", "install"]
    if _choice == "2":
        _cmd += ["-e"]
    _cmd += [str(_repo), "--no-build-isolation", "-Cbuilddir=build"]

    print(f"[gwswex] Running: {' '.join(_cmd)} ...", file=sys.stderr)
    if subprocess.run(_cmd, cwd=_repo).returncode != 0:
        print("[gwswex] pip install failed.", file=sys.stderr)
        return False

    return True


# Eagerly probe for the compiled Fortran extension so any ImportError surfaces
# at `import gwswex` time rather than later at model.init() time.
try:
    from . import f_gwswex as _f_gwswex_probe  # type: ignore[attr-defined]

    del _f_gwswex_probe
except ImportError as _ext_err:
    if _gwswex_build_prompt(_ext_err):
        raise ImportError(
            "gwswex extension built and installed — please restart Python and re-import gwswex."
        ) from None
    raise
finally:
    del _gwswex_build_prompt


from .config import (
    ETStressParams,
    InitialConditions,
    LateralFluxes,
    Material,
    ModelParams,
    RootGrowthModel,
    RootParams,
    SolverConfig,
    SpatialDomain,
    TemporalDomain,
    VanGenuchtenParams,
    Vegetation,
)
from .model import GWSWEXmodel

__all__ = [
    "GWSWEXmodel",
    "Material",
    "VanGenuchtenParams",
    "Vegetation",
    "ETStressParams",
    "ModelParams",
    "RootParams",
    "RootGrowthModel",
    "SpatialDomain",
    "TemporalDomain",
    "SolverConfig",
    "InitialConditions",
    "LateralFluxes",
]
