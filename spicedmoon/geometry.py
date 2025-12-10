"""
Compute main lunar geometries
"""
import os
from typing import List, Tuple

import numpy as np
import spiceypy as spice

from .angular import get_zn_az

from .basics import furnsh_safer, get_radii_moon
from .types import MoonData
from .constants import BASIC_KERNELS, MOON_KERNELS
from .coordinates import to_planetographic_multiple, to_rectangular_multiple


def _get_sel_lon_lat_intercept(pos_moonref: np.ndarray):
    m_eq_rad, m_pol_rad = get_radii_moon(ignore_bodvrd=True)
    flattening = (m_eq_rad - m_pol_rad) / m_eq_rad
    x, y, z = pos_moonref
    # Intersection ray center-body with the moon ellipsoid
    k = 1.0 / np.sqrt((x * x + y * y) / (m_eq_rad**2) + (z * z) / (m_pol_rad**2))
    spoint = np.array([k * x, k * y, k * z])
    sel_lon, sel_lat, _ = spice.recpgr("MOON", spoint, m_eq_rad, flattening)
    return sel_lon, sel_lat


def _get_sel_lon_lat_simple(pos_moonref: np.ndarray):
    sel_lon = np.arctan2(pos_moonref[1], pos_moonref[0])
    sel_lat = np.arctan2(
        pos_moonref[2],
        np.sqrt(pos_moonref[0] * pos_moonref[0] + pos_moonref[1] * pos_moonref[1]),
    )
    return sel_lon, sel_lat


def _get_distance_moon(pos_moonref: np.ndarray):
    distance = np.sqrt(
        pos_moonref[0] * pos_moonref[0]
        + pos_moonref[1] * pos_moonref[1]
        + pos_moonref[2] * pos_moonref[2]
    )
    return distance


def _get_moon_data_xyzs(
    xyz: Tuple[float, float, float],
    et: float,
    source_frame: str,
    target_frame: str,
    angular_frame: str,
    intercept_ellipsoid: bool,
) -> MoonData:
    sun_pos_moonref, lightime = spice.spkpos("SUN", et, target_frame, "NONE", "MOON")
    # sun_pos_satref, lighttime = spice.spkpos("SUN", et, source_frame, "NONE", "EARTH")
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
    if intercept_ellipsoid:
        sel_lon_sun, sel_lat_sun = _get_sel_lon_lat_intercept(sun_pos_moonref)
        sel_lon_sat, sel_lat_sat = _get_sel_lon_lat_intercept(sat_pos_moonref)
    else:
        sel_lon_sun, sel_lat_sun = _get_sel_lon_lat_simple(sun_pos_moonref)
        sel_lon_sat, sel_lat_sat = _get_sel_lon_lat_simple(sat_pos_moonref)
    sel_lon_sat, sel_lat_sat = np.array([sel_lon_sat, sel_lat_sat]) * 180.0 / np.pi
    # distances
    distance_sun_moon = _get_distance_moon(sun_pos_moonref)
    dist_sun_moon_au = spice.convrt(distance_sun_moon, "KM", "AU")
    distance_sat_moon = _get_distance_moon(sat_pos_moonref)
    # zn az
    ang_rotation = spice.pxform(source_frame, angular_frame, et)
    sat_pos_angref = spice.mxv(ang_rotation, sat_pos_translate)
    plt = to_planetographic_multiple(
        [xyz],
        obs_body,
        [et],
        source_frame,
        angular_frame,
    )
    zn, az = get_zn_az(
        -sat_pos_angref, in_sez=False, latitude=plt[0][0], longitude=plt[0][1]
    )
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


def _get_moon_datas_xyzs(
    xyzs: List[Tuple[float, float, float]],
    dts: List[str],
    source_frame: str,
    target_frame: str,
    angular_frame: str,
    intercept_ellipsoid: bool,
) -> List[MoonData]:
    mds = []
    for xyz, dt in zip(xyzs, dts):
        et = spice.str2et(dt)
        md = _get_moon_data_xyzs(
            xyz, et, source_frame, target_frame, angular_frame, intercept_ellipsoid
        )
        et_2 = et + 1
        md2 = _get_moon_data_xyzs(
            xyz, et_2, source_frame, target_frame, angular_frame, intercept_ellipsoid
        )
        if md2.mpa_deg < md.mpa_deg:
            md.mpa_deg = -md.mpa_deg
        mds.append(md)
    return mds


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
        Observer rectangular positions in km.
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
    angular : str
        Name of the EARTH frame which the calculated zenith and azimuth point will be referencing.
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
    mds = _get_moon_datas_xyzs(
        xyzs, dts, source_frame, target_frame, angular_frame, intercept_ellipsoid
    )
    spice.kclear()
    return mds


def get_moon_datas_llhs(
    llhs: List[Tuple[float, float, float]],
    dts: List[str],
    kernels_path: str,
    body: str = "EARTH",
    source_planetographic_frame: str = "ITRF93",
    source_rectangular_frame: str = "J2000",
    target_frame: str = "MOON_ME",
    angular_frame: str = "ITRF93",
    intercept_ellipsoid: bool = True,
) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox, without using intermediate custom kernels.
    Accepts planetographic coordinates, that will be internally transformed into rectangular ones.

    llhs: list of tuple of 3 floats
        Observer geometrical positions in decimal degrees and km.
    dts : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    body: str
        Body the planetographic coordinates reference.
    source_planetographic_frame : str
        Name of the EARTH or MOON frame the planetographic coordinates are.
    source_rectangular_frame: str
        Name of the frame the rectangular coordinates that reference the body will be.
    target_frame : str
        Name of the MOON frame which the location point will be referencing.
    angular : str
        Name of the EARTH frame which the calculated zenith and azimuth point will be referencing.
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
    ets = spice.str2et(dts)
    xyzs = [
        to_rectangular_multiple(
            [llh], body, [et], source_planetographic_frame, source_rectangular_frame
        )[0]
        for llh, et in zip(llhs, ets)
    ]
    mds = _get_moon_datas_xyzs(
        xyzs,
        dts,
        source_rectangular_frame,
        target_frame,
        angular_frame,
        intercept_ellipsoid,
    )
    spice.kclear()
    return mds
