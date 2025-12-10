"""
Compute angular coordinates like zenith and azimuth, and related calculations.
"""
from typing import Tuple, overload

import numpy as np
import spiceypy as spice



@overload
def get_zn_az(
    state_pos_zenith: np.ndarray,
    *,
    in_sez: bool = True,
) -> Tuple[float, float]:
    ...

@overload
def get_zn_az(
    state_pos_zenith: np.ndarray,
    *,
    in_sez: bool = False,
    latitude: float = None,
    longitude: float = None,
) -> Tuple[float, float]:
    ...

def get_zn_az(
    state_pos_zenith: np.ndarray,
    *,
    in_sez: bool = False,
    latitude: float = None,
    longitude: float = None,
) -> Tuple[float, float]:
    """
    Calculate the zenith and azimuth for a position of a target body, relative to an observing body

    Parameters
    ----------
    state_pos_zenith: np.ndarray
        3D position of the target relative to the observer (first 3 elements of
        a SPICE state vector).
        If `in_sez` is False, this vector is assumed to be expressed in a
        body-fixed frame of the observing body (e.g. IAU_EARTH), with axes
        aligned to the planet's rotation and prime meridian.
        If `in_sez` is True, this vector is assumed to already be expressed
        in the local topocentric SEZ frame:
          X = S (South),
          Y = E (East),
          Z = Z (local zenith).
    in_sez : bool
        - False (default): `state_pos_zenith` is in a body-fixed frame and will
          be rotated into the local SEZ frame using `latitude` and `longitude`.
        - True: `state_pos_zenith` is already in SEZ and no rotation is applied.
    latitude : float
        Geographic latitude of the observer, north positive. Required if `in_sez` is False.
    longitude : float
        Geographic longitude of the observer, east positive. Required if `in_sez` is False.

    Returns
    -------
    zenith: float
        Zenith distance of the target in decimal degrees.
    azimuth: float
        Azimuth of the target in decimal degrees, measured from North towards East:
        0 deg = North, 90 deg = East, 180 deg = South, 270 deg = West.
    """
    if not in_sez:
        if longitude is None or latitude is None:
            raise ValueError(
                "latitude and longitude must be provided when `in_sez` is False"
            )
        colat = get_colat_deg(latitude)
        lon_rad = ((longitude % 180) + 180) * spice.rpd()
        colat_rad = colat * spice.rpd()
        bf2tp = spice.eul2m(-lon_rad, -colat_rad, 0, 3, 2, 3)
        state_pos_zenith = spice.mtxv(bf2tp, state_pos_zenith)
    _, longi, lati = spice.reclat(state_pos_zenith)
    zenith = 90.0 - lati * spice.dpr()
    azimuth = 180 - longi * spice.dpr()
    return zenith, azimuth


def get_colat_deg(lat_deg: float) -> float:
    """
    Convert the latitude into colatitude.

    Parameters
    ----------
    lat: float
        Latitude in decimal degrees

    Returns
    -------
    colat: float
        Colatitude associated to `lat`.
    """
    return 90 - (lat_deg % 90)
