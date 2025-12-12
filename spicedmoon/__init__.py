"""
Calculation of observer-lunar geometries using NASA's SPICE toolbox.
"""
from typing import TYPE_CHECKING as _TYPE_CHECKING

from .custombody.preexisting import get_moon_datas_from_extra_kernels
from .geometry import (
    get_moon_datas_xyzs, get_moon_datas_llhs
)
from .custombody.geotic import get_moon_datas
from .custombody.selenic import get_moon_datas_from_moon
from .heliac import get_sun_moon_datas
from .types import MoonData, MoonSunData


def __getattr__(name):
    if name == "spicedmoon":
        from . import deprecated as _spicedmoon
        return _spicedmoon
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if _TYPE_CHECKING:
    from . import deprecated
