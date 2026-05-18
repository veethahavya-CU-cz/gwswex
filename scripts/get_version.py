#!/usr/bin/env python3
"""Extract version from pyproject.toml.

This script is called by meson.build to establish pyproject.toml as the
single source of truth (SOT) for the package version.
"""

import sys
import tomllib


def get_version():
    """Read version from [project] section of pyproject.toml."""
    with open("pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    return data["project"]["version"]


if __name__ == "__main__":
    print(get_version(), end="")
