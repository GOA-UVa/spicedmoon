"""
Implementation of the now deprecated `spicedmoon.spicedmoon` submodule.
It will be eventually removed, but for now it preserves 1.0.13 functionalities that
should remain temporarily.
"""
from typing import List, Tuple
import warnings

from .custombody.preexisting import get_moon_datas_from_extra_kernels
from .geometry import (
    get_moon_datas_xyzs,
)
from .custombody.geotic import get_moon_datas
from .custombody.selenic import get_moon_datas_from_moon
from .heliac import get_sun_moon_datas
from .types import MoonData, MoonSunData
from .basics import furnsh_safer as _furnsh_safer


warnings.warn(
    "`spicedmoon.spicedmoon` is deprecated and will be removed soon. Use directly `spicedmoon`.",
    FutureWarning,
    stacklevel=3,
)


def get_moon_datas_xyzs_no_zenith_azimuth(
    xyzs: List[Tuple[float, float, float]],
    dts: List[str],
    kernels_path: str,
    source_frame: str = "J2000",
    target_frame: str = "MOON_ME",
):
    warnings.warn(
        "`spicedmoon.spicedmoon` is deprecated. Use directly `spicedmoon.geometry.get_moon_datas_xyzs`.",
        FutureWarning,
        stacklevel=2,
    )
    return get_moon_datas_xyzs(
        xyzs, dts, kernels_path, source_frame, target_frame, "ITRF93", False
    )
