"""Version information for GWSWEX."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gwswex")
except PackageNotFoundError:
    __version__ = "unknown"
