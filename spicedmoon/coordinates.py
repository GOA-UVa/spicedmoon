"""
Switch between planetographic and rectangular coordinates.
"""
from typing import List, Tuple

import spiceypy as spice
import numpy as np


def to_rectangular_same_frame(
    latlonheights: List[Tuple[float, float, float]],
    body: str,
) -> List[List[float]]:
    """Convert planetographic coordinates to rectangular, using the same reference frame.
    
    Parameters
    ----------
    latlonheights: list of tuples of float
        Each tuple has 3 float values:
        - lat: Latitude in decimal degrees
        - lon: Longitude in decimal degrees
        - hhh: Height (distance to surface) in kilometers
    body: str
        Name of the standard body the coordinates system references.

    Returns
    -------
    pos: list of list of float
        Each inner list has 3 values: The rectangular coordinates in kilometers.
    """
    _, radios = spice.bodvrd(body, "RADII", 3)
    eq_rad = radios[0]  # Equatorial Radius
    pol_rad = radios[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    poss_iaus = []
    for llh in latlonheights:
        pos_iau = spice.georec(
            spice.rpd() * llh[1],
            spice.rpd() * llh[0],
            llh[2],
            eq_rad,
            flattening,
        )
        poss_iaus.append(pos_iau)
    poss_iaus = list(map(lambda n: n, poss_iaus))
    return poss_iaus


def to_planetographic_same_frame(
    xyz_list: List[Tuple[float]],
    body: str,
) -> List[List[float]]:
    """Convert rectangular coordinates to planetographic, using the same reference frame.
    
    Parameters
    ----------
    xyz_list: list of tuples of float
        Each tuple has 3 float values, corresponding to the x, y and z coordinates in kilometers.
    body: str
        Name of the standard body the coordinates system references.

    Returns
    -------
    llh: list of list of float
        Each inner list has 3 float values:
        - lat: Latitude in decimal degrees
        - lon: Longitude in decimal degrees
        - hhh: Height (distance to surface) in kilometers
    """
    _, radii = spice.bodvrd(body, "RADII", 3)
    eq_rad = radii[0]  # Equatorial Radius
    pol_rad = radii[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    llh_list = []
    for xyz in xyz_list:
        pos_iau = np.array(list(xyz))
        llh = spice.recgeo(pos_iau, eq_rad, flattening)
        llh_list.append(llh)
    for i, llh in enumerate(llh_list):
        lat = llh[1] * spice.dpr()
        lon = llh[0] * spice.dpr()
        alt = llh[2]
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360
        llh_list[i] = (lat, lon, alt)
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
    ets: List[float],
    source_frame: str = "IAU_EARTH",
    target_frame: str = "J2000",
):
    """Convert planetographic coordinates to rectangular, also changing the reference frame.
    
    Parameters
    ----------
    latlonheights: list of tuples of float
        Each tuple has 3 float values:
        - lat: Latitude in decimal degrees
        - lon: Longitude in decimal degrees
        - hhh: Height (distance to surface) in kilometers
    body: str
        Name of the standard body the coordinates system references.
    ets: list of float
        Timestamps at which the positions coordinates are being changed.
    source_frame: str
        Original reference frame `xyz_list` is on.
    target_frame: str
        Target reference frame output will be on.
    Returns
    -------
    pos: list of list of float
        Each inner list has 3 values: The rectangular coordinates in kilometers.
    """
    _, radios = spice.bodvrd(body, "RADII", 3)
    eq_rad = radios[0]  # Equatorial Radius
    pol_rad = radios[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    poss_iaus = []
    for llh, et in zip(latlonheights, ets):
        pos_iau = spice.georec(
            spice.rpd() * llh[1],
            spice.rpd() * llh[0],
            llh[2],
            eq_rad,
            flattening,
        )
        poss_iaus.append(_change_frames(pos_iau, source_frame, target_frame, et))
    poss_iaus = list(poss_iaus)
    return poss_iaus


def to_planetographic_multiple(
    xyz_list: List[Tuple[float]],
    body: str,
    ets: List[float],
    source_frame: str = "J2000",
    target_frame: str = "IAU_EARTH",
) -> List[List[float]]:
    """Convert rectangular coordinates to planetographic, also changing the reference frame.
    
    Parameters
    ----------
    xyz_list: list of tuples of float
        Each tuple has 3 float values, corresponding to the x, y and z coordinates in kilometers.
    body: str
        Name of the standard body the coordinates system references.
    ets: list of float
        Timestamps at which the positions coordinates are being changed.
    source_frame: str
        Original reference frame `xyz_list` is on.
    target_frame: str
        Target reference frame output will be on.
    Returns
    -------
    llh: list of list of float
        Each inner list has 3 float values:
        - lat: Latitude in decimal degrees
        - lon: Longitude in decimal degrees
        - hhh: Height (distance to surface) in kilometers
    """
    _, radii = spice.bodvrd(body, "RADII", 3)
    eq_rad = radii[0]  # Equatorial Radius
    pol_rad = radii[2]  # Polar radius
    flattening = (eq_rad - pol_rad) / eq_rad
    llh_list = []
    for xyz, et in zip(xyz_list, ets):
        pos_iau = np.array(list(xyz))
        pos_iau_proc = _change_frames(pos_iau, source_frame, target_frame, et)
        llh = spice.recgeo(pos_iau_proc, eq_rad, flattening)
        llh_list.append(llh)
    for i, llh in enumerate(llh_list):
        lat = llh[1] * spice.dpr()
        lon = llh[0] * spice.dpr()
        alt = llh[2]
        while lon < -180:
            lon += 360
        while lon > 180:
            lon -= 360
        llh_list[i] = (lat, lon, alt)
    return llh_list
