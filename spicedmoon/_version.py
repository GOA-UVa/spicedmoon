"""
spicedmoon version handler module

Handles versioning for the `spicedmoon` package.
"""
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("spicedmoon")
except PackageNotFoundError:
    __version__ = "0.0.0"
