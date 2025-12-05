import os
from typing import Union, List
from datetime import datetime

import numpy as np
import spiceypy as spice

from .types import MoonSunData
from .constants import BASIC_KERNELS, MOON_KERNELS
from .basics import dt_to_str, furnsh_safer, get_radii_moon


def get_sun_moon_data(
    utc_time: str,
    ignore_bodvrd: bool = True,
):
    et_date = spice.str2et(utc_time)
    m_eq_rad, m_pol_rad = get_radii_moon(ignore_bodvrd)
    flattening = (m_eq_rad - m_pol_rad) / m_eq_rad
    # Calculate selenographic longitude of sun
    sun_spoint, _, _ = spice.subslr(
        "INTERCEPT/ELLIPSOID", "MOON", et_date, "MOON_ME", "NONE", "SUN"
    )
    lon_sun_rad, lat_sun_rad, _ = spice.recpgr("MOON", sun_spoint, m_eq_rad, flattening)

    # Calculate the distance between sun and moon (AU)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", "SUN")
    dist_sun_moon_km = np.sqrt(state[0] ** 2 + state[1] ** 2 + state[2] ** 2)
    dist_sun_moon_au = spice.convrt(dist_sun_moon_km, "KM", "AU")

    limit_lat_rad = np.pi / 2
    limit_lon_rad = np.pi
    if lat_sun_rad > limit_lat_rad:
        lat_sun_rad = limit_lat_rad + (limit_lat_rad - lat_sun_rad)
        lon_sun_rad -= limit_lon_rad
    elif lat_sun_rad < -limit_lat_rad:
        lat_sun_rad = -limit_lat_rad - (limit_lat_rad + lat_sun_rad)
        lon_sun_rad += limit_lon_rad

    while lon_sun_rad > limit_lon_rad:
        lon_sun_rad -= limit_lon_rad * 2
    while lon_sun_rad < -limit_lon_rad:
        lon_sun_rad += limit_lon_rad * 2

    return MoonSunData(lon_sun_rad, lat_sun_rad, dist_sun_moon_km, dist_sun_moon_au)


def get_sun_moon_datas(
    times: Union[List[str], List[datetime]],
    kernels_path: str,
    ignore_bodvrd: bool = True,
) -> List[MoonSunData]:
    """
    Obtain solar selenographic coordinates of at multiple times.

    times : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    ignore_bodvrd : bool
        Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
        1738.1 and 1736
    """
    utc_times = dt_to_str(times)
    if len(utc_times) == 0:
        return []
    kernels = BASIC_KERNELS + MOON_KERNELS
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)

    msds = []
    for utc_time in utc_times:
        msd = get_sun_moon_data(utc_time, ignore_bodvrd)
        msds.append(msd)
    spice.kclear()
    return msds
