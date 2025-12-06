from typing import List, Tuple

import spiceypy as spice
import numpy as np


def to_rectangular_same_frame(
    latlonheights: List[Tuple[float, float, float]],
    body: str,
):  # h in meters
    _, radios = spice.bodvrd(body, "RADII", 3)
    eq_rad = radios[0]  # Equatorial Radius
    pol_rad = radios[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    poss_iaus = []
    for llh in latlonheights:
        pos_iau = spice.georec(
            spice.rpd() * llh[1],
            spice.rpd() * llh[0],
            llh[2] / 1000,
            eq_rad,
            flattening,
        )
        poss_iaus.append(pos_iau)
    poss_iaus = list(map(lambda n: n * 1000, poss_iaus))
    return poss_iaus  # in meters


def to_planetographic_same_frame(
    xyz_list: List[Tuple[float]],
    body: str,
):
    _, radii = spice.bodvrd(body, "RADII", 3)
    eq_rad = radii[0]  # Equatorial Radius
    pol_rad = radii[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    llh_list = []  # alt km
    for xyz in xyz_list:
        pos_iau = np.array(list(map(lambda n: n / 1000, xyz)))
        llh = spice.recgeo(pos_iau, eq_rad, flattening)
        llh_list.append(llh)
    for i, llh in enumerate(llh_list):
        lat = llh[1] * spice.dpr()
        lon = llh[0] * spice.dpr()
        alt = llh[2] * 1000
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360
        llh_list[i] = (lat, lon, alt)
    # alt in meters
    return llh_list


def _change_frames(
    coords: np.ndarray, source_frame: str, target_frame: str, et: float
) -> np.ndarray:
    if "MOON" not in target_frame:
        trans_matrix = spice.pxform(source_frame, target_frame, et)
        return spice.mxv(trans_matrix, coords)
    moon_pos_satref, _ = spice.spkpos("MOON", et, source_frame, "NONE", "EARTH")
    rotation = spice.pxform(source_frame, target_frame, et)
    # set moon center as zero point
    sat_pos_translate = np.zeros(3)
    sat_pos_translate[0] = coords[0] - moon_pos_satref[0]
    sat_pos_translate[1] = coords[1] - moon_pos_satref[1]
    sat_pos_translate[2] = coords[2] - moon_pos_satref[2]
    return spice.mxv(rotation, sat_pos_translate)


def to_rectangular_multiple(
    latlonheights: List[Tuple[float, float, float]],
    body: str,
    dts: List[str],
    source_frame: str = "IAU_EARTH",
    target_frame: str = "J2000",
):  # h in meters
    _, radios = spice.bodvrd(body, "RADII", 3)
    eq_rad = radios[0]  # Equatorial Radius
    pol_rad = radios[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    poss_iaus = []
    ets = spice.str2et(dts)
    for llh, et in zip(latlonheights, ets):
        pos_iau = spice.georec(
            spice.rpd() * llh[1],
            spice.rpd() * llh[0],
            llh[2] / 1000,
            eq_rad,
            flattening,
        )
        poss_iaus.append(_change_frames(pos_iau, source_frame, target_frame, et))
    poss_iaus = list(map(lambda n: n * 1000, poss_iaus))
    return poss_iaus  # in meters


def to_planetographic_multiple(
    xyz_list: List[Tuple[float]],
    body: str,
    dts: List[str],
    source_frame: str = "J2000",
    target_frame: str = "IAU_EARTH",
):  # in meters
    _, radii = spice.bodvrd(body, "RADII", 3)
    eq_rad = radii[0]  # Equatorial Radius
    pol_rad = radii[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    llh_list = []  # alt km
    ets = spice.str2et(dts)
    for xyz, et in zip(xyz_list, ets):
        pos_iau = np.array(list(map(lambda n: n / 1000, xyz)))
        pos_iau_proc = _change_frames(pos_iau, source_frame, target_frame, et)
        llh = spice.recgeo(pos_iau_proc, eq_rad, flattening)
        llh_list.append(llh)
    for i, llh in enumerate(llh_list):
        lat = llh[1] * spice.dpr()
        lon = llh[0] * spice.dpr()
        alt = llh[2] * 1000
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360
        llh_list[i] = (lat, lon, alt)
    return llh_list
