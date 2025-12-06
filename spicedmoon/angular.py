import numpy as np
import spiceypy as spice


def get_zn_az(
    state_pos_zenith: np.ndarray,
    correct_rotating: bool = False,
    longitude: float = None,
    colat: float = None,
):
    """
    Calculate the zenith and azimuth for a position of a target body, relative to an observing body

    Parameters
    ----------
    state_pos_zenith: np.ndarray
        The position (3 first elements of state) of a target body relative to an observing body
    correct_rotating : bool
        Correct the coordinates rotating them into the local SEZ orientation.
    longitude : float
        Geographic longitude of the observer point. Needed only if correcting coordinates.
    colat : float
        Geographic colatitude of the observer point. Needed only if correcting coordinates.

    Returns
    -------
    zenith: float
        Zenith of the target body in decimal degrees.
    azimuth: float
        Azimuth of the target body in decimal degrees.
    """
    if correct_rotating:
        if longitude is None or colat is None:
            raise ValueError(
                "longitude and colat must be provided when correct_rotating=True"
            )
        lon_rad = (longitude + 180) * spice.rpd()
        colat_rad = colat * spice.rpd()
        bf2tp = spice.eul2m(-lon_rad, -colat_rad, 0, 3, 2, 3)
        state_pos_zenith = spice.mtxv(bf2tp, state_pos_zenith)
    _, longi, lati = spice.reclat(state_pos_zenith)
    zenith = 90.0 - lati * spice.dpr()
    azimuth = 180 - longi * spice.dpr()
    return zenith, azimuth


def get_colat_deg(lat: float) -> float:
    return 90 - (lat % 90)
