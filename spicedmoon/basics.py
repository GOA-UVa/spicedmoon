import os
import time
from typing import List, Union, Tuple
from datetime import datetime, timezone

import numpy as np
import spiceypy as spice

from .constants import CUSTOM_KERNEL_NAME


def furnsh_safer(k_path: str):
    """
    Performs SPICE's `furnsh_c`, but in case that it fails it tries again after a small time
    interval.

    Furnsh very rarely crashes, but it can be solved trying again.

    Parameters
    ----------
    k_path : str
        Path of the kernel to load.
    """
    try:
        spice.furnsh(k_path)
    except:
        time.sleep(2)
        spice.furnsh(k_path)


def dt_to_str(dts: Union[List[datetime], List[str]]) -> List[str]:
    """Convert a list of datetimes into a list of string dates in a valid format.

    Parameters
    ----------
    dts: list of datetimes | list of str
        List of datetimes that will be converted to utc_times. They must be timezone aware.
        A list of already converted strings can be given instead, and it will be returned without
        change.

    Returns
    -------
    utc_times: list of str
        List of the datetimes in a valid string format for SPICE.
    """
    utc_times = []
    for dt in dts:
        if isinstance(dt, datetime):
            dt_utc = dt.astimezone(timezone.utc)
            utc_times.append(dt_utc.strftime("%Y-%m-%d %H:%M:%S"))
        else:
            utc_times.append(dt)
    return utc_times


def remove_custom_kernel_file(kernels_path: str) -> None:
    """Remove the custom SPK kernel file if it exists

    Parameters
    ----------
    kernels_path : str
        Path where the SPICE kernels are stored
    """
    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    if os.path.exists(custom_kernel_path):
        os.remove(custom_kernel_path)


def get_radii_moon(ignore_bodvrd: bool = True) -> Tuple[float, float]:
    """
    Obtain moon radii.

    Parameters
    ----------
    ignore_bodvrd: bool
        If True, assign default values instead of obtaining them through spice's bodvrd
        which yields less accurate lunar radii. True by default.

    Returns
    -------
    eq_rad: float
        Equatorial radius of the Moon.
    pol_rad: float
        Polar radius of the Moon.
    """
    eq_rad, pol_rad = 1738.1, 1736
    if not ignore_bodvrd:
        # The ones obtained with bodvrd are not correct
        _, radii_moon = spice.bodvrd("MOON", "RADII", 3)
        eq_rad, pol_rad = radii_moon[0], radii_moon[2]
    return eq_rad, pol_rad


def _calculate_states(
    ets: np.ndarray,
    pos_iau: np.ndarray,
    delta_t: float,
    source_frame: str,
    target_frame: str,
) -> np.ndarray:
    """
    Returns a ndarray containing the states of a point referencing the target frame.

    The states array is a time-ordered array of geometric states (x, y, z, dx/dt, dy/dt, dz/dt,
    in kilometers and kilometers per second) of body relative to center, specified relative
    to frame. Useful for spice function "spkw09_c", for example.

    Parameters
    ----------
    ets : np.ndarray
        Array of TDB seconds from J2000 (et dates) of which the data will be taken
    pos_iau : np.ndarray
        Rectangular coordinates of the point, referencing IAU frame.
    delta_t : float
        TDB seconds between states
    source_frame : str
        Name of the frame to transform from.
    target_frame : str
        Name of the frame which the location will be referencing.

    Returns
    -------
    ndarray of float
        ndarray containing the states calculated
    """
    num_coordinates = 3
    n_state_attributes = 6
    states = np.zeros((len(ets), n_state_attributes))
    for i, et_value in enumerate(ets):
        states[i, :num_coordinates] = np.dot(
            spice.pxform(source_frame, target_frame, et_value), pos_iau
        )

    for i in range(len(ets) - 1):
        states[i, num_coordinates:] = (
            states[i + 1, :num_coordinates] - states[i, :num_coordinates]
        ) / delta_t

    pos_np1 = np.dot(
        spice.pxform(source_frame, target_frame, ets[-1] + delta_t), pos_iau
    )
    states[-1, num_coordinates:] = (pos_np1 - states[-1, :num_coordinates]) / delta_t
    return states


class Location:
    """
    Data for the creation of an observer point on a body's surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on the body's surface
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """

    __slots__ = ["point_id", "states"]

    def __init__(
        self,
        point_id: int,
        body: str,
        lat: float,
        lon: float,
        altitude: float,
        eq_rad: float,
        pol_rad: float,
        ets: np.ndarray,
        delta_t: float,
        source_frame: str,
        target_frame: str,
    ):
        """
        Parameters
        ----------
        point_id: int
            ID code that will be associated with the point.
        body: str
            Name of the body the point is located at.
        lat : float
            Planetographic latitude of the observer point
        lon : float
            Planetographic longitude of the observer point
        altitude : float
            Altitude over the sea level in meters.
        eq_rad: float
            Body's equatorial radius
        pol_rad: float
            Body's polar radius
        ets : np.ndarray
            Array of TDB seconds from J2000 (et dates) of which the data will be taken.
        delta_t : float
            TDB seconds between states
        source_frame : str
            Name of the frame to transform from.
        target_frame : str
            Name of the frame which the location will be referencing.
        """
        self.point_id = point_id
        flattening = (eq_rad - pol_rad) / eq_rad
        alt_km = altitude / 1000
        pos_iau = spice.pgrrec(
            body, np.radians(lon), np.radians(lat), alt_km, eq_rad, flattening
        )
        self.states = _calculate_states(
            ets, pos_iau, delta_t, source_frame, target_frame
        )
