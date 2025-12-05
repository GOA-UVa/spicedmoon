import os
from typing import List

import numpy as np
import spiceypy as spice

from .constants import (
    CUSTOM_KERNEL_NAME,
    DEFAULT_OBSERVER_FRAME,
    DEFAULT_OBSERVER_NAME,
    DEFAULT_OBSERVER_ZENITH_NAME,
    BASIC_KERNELS,
    MOON_KERNELS,
)
from .types import MoonData
from .basics import get_radii_moon, furnsh_safer
from .heliac import get_sun_moon_data


def get_moon_data(
    utc_time: str,
    observer_name: str = DEFAULT_OBSERVER_NAME,
    observer_frame: str = DEFAULT_OBSERVER_FRAME,
    observer_zenith_name: str = DEFAULT_OBSERVER_ZENITH_NAME,
    correct_zenith_azimuth: bool = False,
    longitude: float = 0,
    colat: float = 0,
    ignore_bodvrd: bool = True,
) -> MoonData:
    """Calculation of the moon data for the given utc_time for the loaded observer

    Parameters
    ----------
    utc_time : str
        Time at which the lunar data will be calculated, in a valid UTC DateTime format
    observer_name : str
        Name of the body of the observer that should be loaded from the extra kernels.
        By default is "Observer", in which case it shouldn't be loaded from the extra
        kernels but from the custom kernel.
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
    observer_zenith_name : str
        The observer used for the zenith and azimuth calculation. By default it's "EARTH".
    correct_zenith_azimuth : bool
        In case that it's calculated without using the extra kernels, the coordinates should be
        corrected rotating them into the correct location.
    longitude : float
        Geographic longitude of the observer point. Used if it's calculated without using the
        extra kernels.
    colat : float
        Geographic colatitude of the observer point. Used if it's calculated without using the
        extra kernels.
    ignore_bodvrd : bool
        Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
        1738.1 and 1736
    Returns
    -------
    MoonData
        Moon data obtained from SPICE toolbox
    """
    et_date = spice.str2et(utc_time)

    m_eq_rad, m_pol_rad = get_radii_moon(ignore_bodvrd)
    flattening = (m_eq_rad - m_pol_rad) / m_eq_rad

    # Calculate moon zenith and azimuth
    state_zenith, _ = spice.spkezr(
        "MOON", et_date, observer_frame, "NONE", observer_zenith_name
    )
    rectan_zenith = np.split(state_zenith, 2)[0]
    if correct_zenith_azimuth:
        lon_rad = (longitude + 180) * spice.rpd()
        colat_rad = colat * spice.rpd()
        bf2tp = spice.eul2m(-lon_rad, -colat_rad, 0, 3, 2, 3)
        rectan_zenith = spice.mtxv(bf2tp, rectan_zenith)

    _, longi, lati = spice.reclat(rectan_zenith)

    zenith = 90.0 - lati * spice.dpr()
    azimuth = 180 - longi * spice.dpr()

    # Calculate moon phase angle
    spoint, _, _ = spice.subpnt(
        "INTERCEPT/ELLIPSOID", "MOON", et_date, "MOON_ME", "NONE", observer_name
    )
    phase = spice.phaseq(et_date, "MOON", "SUN", observer_name, "NONE")
    phase = phase * spice.dpr()

    et_date_2 = et_date + 1
    phase2 = spice.phaseq(et_date_2, "MOON", "SUN", observer_name, "NONE")
    phase2 = phase2 * spice.dpr()
    if phase2 < phase:
        phase = -phase

    # Calculate selenographic coordinates of the observer
    lon_obs, lat_obs, _ = spice.recpgr("MOON", spoint, m_eq_rad, flattening)
    lon_obs = lon_obs * spice.dpr()
    lat_obs = lat_obs * spice.dpr()

    # Calculate the distance between observer and moon (KM)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", observer_name)
    dist_obs_moon = np.sqrt(state[0] ** 2 + state[1] ** 2 + state[2] ** 2)

    smd = get_sun_moon_data(utc_time, ignore_bodvrd)
    lon_sun_rad = smd.lon_sun_rad
    lat_sun_rad = smd.lat_sun_rad
    dist_sun_moon_km = smd.dist_sun_moon_km
    dist_sun_moon_au = smd.dist_sun_moon_au

    limit_lat = 90
    limit_lon = 180
    if lat_obs > limit_lat:
        lat_obs = limit_lat + (limit_lat - lat_obs)
        lon_obs -= limit_lon
    elif lat_obs < -limit_lat:
        lat_obs = -limit_lat - (limit_lat + lat_obs)
        lon_obs += limit_lon

    while lon_obs > limit_lon:
        lon_obs -= limit_lon * 2
    while lon_obs < -limit_lon:
        lon_obs += limit_lon * 2

    moon_data = MoonData(
        dist_sun_moon_au,
        dist_sun_moon_km,
        dist_obs_moon,
        lon_sun_rad,
        lat_sun_rad,
        lat_obs,
        lon_obs,
        phase,
        azimuth,
        zenith,
    )
    return moon_data


def get_moon_datas_id(
    utc_times: List[str],
    kernels_path: str,
    observer_id: int,
    observer_frame: str,
    custom_kernel_dir: str,
    correct_zenith_azimuth: bool = False,
    latitude: float = 0,
    longitude: float = 0,
    earth_as_zenith_observer: bool = False,
    ignore_bodvrd: bool = True,
) -> List[MoonData]:
    """Calculation of needed MoonDatas from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    utc_times : list of str
        Times at which the lunar data will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    observer_id : int
        Observer's body ID
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
    custom_kernel_dir: str
        Path where the writable kernel custom.bsp will be stored.
    correct_zenith_azimuth : bool
        In case that it's calculated without using the extra kernels, the coordinates should be
        corrected rotating them into the correct location.
    latitude : float
        Geographic latitude of the observer point.
    longitude : float
        Geographic longitude of the observer point.
    earth_as_zenith_observer : bool
        If True the Earth will be used as the observer for the zenith and azimuth calculation.
        Otherwise it will be the actual observer. By default is False.
    ignore_bodvrd : bool
        Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
        1738.1 and 1736
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    kernels = BASIC_KERNELS + MOON_KERNELS

    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)
    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    furnsh_safer(custom_kernel_path)

    observer_name = DEFAULT_OBSERVER_NAME
    spice.boddef(observer_name, observer_id)
    if earth_as_zenith_observer:
        zenith_observer = "EARTH"
    else:
        zenith_observer = observer_name
    moon_datas = []
    colat = 90 - (latitude % 90)
    lon = longitude % 180
    for utc_time in utc_times:
        new_md = get_moon_data(
            utc_time,
            observer_name,
            observer_frame,
            zenith_observer,
            correct_zenith_azimuth,
            lon,
            colat,
            ignore_bodvrd,
        )
        moon_datas.append(new_md)

    spice.kclear()

    return moon_datas
