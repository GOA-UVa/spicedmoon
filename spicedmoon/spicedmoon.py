"""Spiced Moon

Calculation of lunar data using NASA's SPICE toolbox.

It exports the following functions:

    * get_moon_datas - Calculates needed MoonData from SPICE toolbox
    * get_moon_datas_from_extra_kernels - Calculates needed MoonData from SPICE toolbox
        and using data from extra kernels for the observer body
    * get_sun_moon_datas - Calculates solar selenographic coordinates.
    * get_moon_datas_from_moon - Calculates needed MoonData from SPICE toolbox from selenographic coordinates
"""
from .geoselenic import (
    get_moon_datas_from_extra_kernels,
    get_moon_datas_xyzs_no_zenith_azimuth,
)
from .geotic import get_moon_datas
from .selenic import get_moon_datas_from_moon
from .heliac import get_sun_moon_datas
from .types import MoonData, MoonSunData
