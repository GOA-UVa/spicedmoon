"""Spiced Moon

Calculation of lunar data using NASA's SPICE toolbox.

It exports the following functions:

    * get_moon_datas - Calculates needed MoonData from SPICE toolbox
    * get_moon_datas_from_extra_kernels - Calculates needed MoonData from SPICE toolbox
        and using data from extra kernels for the observer body
    * get_sun_moon_datas - Calculates solar selenographic coordinates.
    * get_moon_datas_from_moon - Calculates needed MoonData from SPICE toolbox from selenographic coordinates
"""
from typing import List, Tuple

from .custombody.preexisting import get_moon_datas_from_extra_kernels
from .geometry import (
    get_moon_datas_xyzs,
)
from .custombody.geotic import get_moon_datas
from .custombody.selenic import get_moon_datas_from_moon
from .heliac import get_sun_moon_datas
from .types import MoonData, MoonSunData


def get_moon_datas_xyzs_no_zenith_azimuth(
    xyzs: List[Tuple[float, float, float]],
    dts: List[str],
    kernels_path: str,
    source_frame: str = "J2000",
    target_frame: str = "MOON_ME",
):
    get_moon_datas_xyzs(
        xyzs, dts, kernels_path, source_frame, target_frame, "ITRF93", False
    )
