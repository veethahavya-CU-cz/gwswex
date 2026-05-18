#!/usr/bin/env python3
"""Remove all Meson build artefacts and compiled extension files.

Deletes the following directories from the project root (if present):
  build/       — release / standard Meson build directory
  build-debug/ — debug Meson build directory
  *.egg-info/  — editable-install metadata created by pip / meson-python
  .mesonpy-*/  — meson-python temporary directories

Safe to run at any time; missing directories are silently skipped.
"""
import shutil
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    for rel in ("build", "build-debug", "dist", "*.egg-info", ".mesonpy-*"):
        if "*" in rel:
            for path in root.glob(rel):
                if path.is_dir():
                    shutil.rmtree(path, ignore_errors=True)
            continue

        path = root / rel
        if path.exists() and path.is_dir():
            shutil.rmtree(path, ignore_errors=True)

    print("Cleaned build artifacts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
