import os
from typing import List, Union, Tuple
from datetime import datetime

import numpy as np
import spiceypy as spice

from .basics import furnsh_safer, dt_to_str, get_colat_deg, get_radii_moon
from .core import get_moon_data_body_ellipsoid, get_zn_az
from .types import MoonData
from .constants import BASIC_KERNELS, MOON_KERNELS
from .coordinates import to_planetographic_multiple


def get_moon_datas_from_extra_kernels(
    times: Union[List[str], List[datetime]],
    kernels_path: str,
    extra_kernels: List[str],
    extra_kernels_path: str,
    observer_name: str,
    observer_frame: str,
    earth_as_zenith_observer: bool = False,
    ignore_bodvrd: bool = True,
) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    times : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    extra_kernels : list of str
        Custom kernels from which the observer body will be loaded, instead of calculating it.
    extra_kernels_path : str
        Folder where the extra kernels are located.
    observer_name : str
        Name of the body of the observer that will be loaded from the extra kernels.
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
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
    base_kernels = BASIC_KERNELS + MOON_KERNELS
    for kernel in base_kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)
    for kernel in extra_kernels:
        k_path = os.path.join(extra_kernels_path, kernel)
        furnsh_safer(k_path)

    if earth_as_zenith_observer:
        zenith_observer = "EARTH"
    else:
        zenith_observer = observer_name
    moon_datas = []
    utc_times = dt_to_str(times)
    for utc_time in utc_times:
        moon_datas.append(
            get_moon_data_body_ellipsoid(
                utc_time,
                observer_name,
                observer_frame,
                zenith_observer,
                ignore_bodvrd=ignore_bodvrd,
            )
        )
    spice.kclear()
    return moon_datas


def _get_moon_datas_xyzs(
    xyz: Tuple[float, float, float],
    et: float,
    source_frame: str,
    target_frame: str,
    angular_frame: str,
    intercept_ellipsoid: bool,
) -> MoonData:
    sun_pos_moonref, lightime = spice.spkpos("SUN", et, target_frame, "NONE", "MOON")
    sun_pos_satref, lighttime = spice.spkpos("SUN", et, source_frame, "NONE", "EARTH")
    obs_body = "EARTH"
    if "MOON" in source_frame and "MOON" in target_frame:
        obs_body = "MOON"
    moon_pos_satref, lightime = spice.spkpos("MOON", et, source_frame, "NONE", obs_body)
    rotation = spice.pxform(source_frame, target_frame, et)
    # set moon center as zero point
    sat_pos_translate = np.zeros(3)
    sat_pos_translate[0] = xyz[0] - moon_pos_satref[0]
    sat_pos_translate[1] = xyz[1] - moon_pos_satref[1]
    sat_pos_translate[2] = xyz[2] - moon_pos_satref[2]
    sat_pos_moonref = spice.mxv(rotation, sat_pos_translate)
    # selenographic coordinates
    # sun
    sel_lon_sun = np.arctan2(sun_pos_moonref[1], sun_pos_moonref[0])
    sel_lat_sun = np.arctan2(
        sun_pos_moonref[2],
        np.sqrt(
            sun_pos_moonref[0] * sun_pos_moonref[0]
            + sun_pos_moonref[1] * sun_pos_moonref[1]
        ),
    )
    distance_sun_moon = np.sqrt(
        sun_pos_moonref[0] * sun_pos_moonref[0]
        + sun_pos_moonref[1] * sun_pos_moonref[1]
        + sun_pos_moonref[2] * sun_pos_moonref[2]
    )
    dist_sun_moon_au = spice.convrt(distance_sun_moon, "KM", "AU")
    # sat
    if intercept_ellipsoid:
        x, y, z = sat_pos_moonref
        m_eq_rad, m_pol_rad = get_radii_moon(ignore_bodvrd=True)
        flattening = (m_eq_rad - m_pol_rad) / m_eq_rad
        # Intersection ray center-observer with the ellipsoid
        k = 1.0 / np.sqrt(
            (x * x + y * y) / (m_eq_rad**2) + (z * z) / (m_pol_rad**2)
        )
        spoint = np.array([k * x, k * y, k * z])
        sel_lon_sat, sel_lat_sat, _ = spice.recpgr("MOON", spoint, m_eq_rad, flattening)
        sel_lon_sat, sel_lat_sat = np.array([sel_lon_sat, sel_lat_sat]) * 180.0 / np.pi
    else:
        sel_lon_sat = np.arctan2(sat_pos_moonref[1], sat_pos_moonref[0]) * 180.0 / np.pi
        sel_lat_sat = (
            np.arctan2(
                sat_pos_moonref[2],
                np.sqrt(
                    sat_pos_moonref[0] * sat_pos_moonref[0]
                    + sat_pos_moonref[1] * sat_pos_moonref[1]
                ),
            )
            * 180.0
            / np.pi
        )
    distance_sat_moon = np.sqrt(
        sat_pos_moonref[0] * sat_pos_moonref[0]
        + sat_pos_moonref[1] * sat_pos_moonref[1]
        + sat_pos_moonref[2] * sat_pos_moonref[2]
    )
    # zn az
    ang_rotation = spice.pxform(source_frame, angular_frame, et)
    sat_pos_angref = spice.mxv(ang_rotation, sat_pos_translate)
    plt = to_planetographic_multiple(
        [xyz * 1000],
        obs_body,
        [spice.et2utc(et, "ISOC", 0)],
        source_frame,
        angular_frame,
    )
    colat = get_colat_deg(plt[0][0])
    lon = plt[0][1] % 180
    zn, az = get_zn_az(-sat_pos_angref, True, lon, colat)
    # phase
    phase = (180.0 / np.pi) * np.arccos(
        (
            sun_pos_moonref[0] * sat_pos_moonref[0]
            + sun_pos_moonref[1] * sat_pos_moonref[1]
            + sun_pos_moonref[2] * sat_pos_moonref[2]
        )
        / (distance_sat_moon * distance_sun_moon)
    )
    return MoonData(
        dist_sun_moon_au,
        distance_sun_moon,
        distance_sat_moon,
        sel_lon_sun,
        sel_lat_sun,
        sel_lat_sat,
        sel_lon_sat,
        phase,
        az,
        zn,
    )


def get_moon_datas_xyzs(
    xyzs: List[Tuple[float, float, float]],
    dts: List[str],
    kernels_path: str,
    source_frame: str = "J2000",
    target_frame: str = "MOON_ME",
    angular_frame: str = "ITRF93",
    intercept_ellipsoid: bool = True,
) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox, without using intermediate custom kernels.

    xyzs: list of tuple of 3 floats
        Observer rectangular positions
    dts : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    source_frame : str
        Name of the EARTH or MOON frame to transform the coordinates from.
    target_frame : str
        Name of the MOON frame which the location point will be referencing.
    intercept_ellipsoid: bool
        Controls how the observer selenographic latitude and longitude are defined.
        If True, they correspond to the sub-observer point on the lunar surface,
        computed by intersecting the observer direction with the Moon reference
        ellipsoid (SPICE "INTERCEPT/ELLIPSOID" behavior).
        If False, they correspond to the angular direction of the observer as seen
        from the Moon center, without intersecting the lunar surface.
    Returns
    -------
    list of MoonData
        List of the calculated MoonDatas
    """
    kernels = BASIC_KERNELS + MOON_KERNELS
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)
    mds = []
    for xyz, dt in zip(xyzs, dts):
        et = spice.str2et(dt)
        md = _get_moon_datas_xyzs(
            xyz, et, source_frame, target_frame, angular_frame, intercept_ellipsoid
        )
        et_2 = et + 1
        md2 = _get_moon_datas_xyzs(
            xyz, et_2, source_frame, target_frame, angular_frame, intercept_ellipsoid
        )
        if md2.mpa_deg < md.mpa_deg:
            md.mpa_deg = -md.mpa_deg
        mds.append(md)
    spice.kclear()
    return mds
