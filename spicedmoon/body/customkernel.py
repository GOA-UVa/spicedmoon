"""
Create and remove the custom kernel used in some `spicedmoon` calculations.
"""
import os
from typing import List

import numpy as np
import spiceypy as spice

from ..constants import CUSTOM_KERNEL_NAME


def _calculate_ets(
    utc_times: List[str], delta_t: int, polynomial_degree: int
) -> np.ndarray:
    """
    Parameters
    ----------
    utc_times: list of str
        Timestamps at which the data will be calculated, in UTC and a valid SPICE format.
    delta_t: int
        TDB seconds between states.
    polynomial_degree: int
        Degree of the lagrange polynomials that will be used to interpolate the states.
    """
    min_states_polynomial = polynomial_degree + 1
    # Min # states that are required to define a polynomial of that degree
    ets = np.array([])
    left_states = int(min_states_polynomial / 2)
    right_states = left_states + min_states_polynomial % 2
    for utc_time in utc_times:
        et0 = spice.str2et(utc_time)
        etprev = et0 - delta_t * left_states
        etf = et0 + delta_t * right_states
        ets_t = np.arange(etprev, etf, delta_t)
        for et_t in ets_t:
            if et_t not in ets:
                index = np.searchsorted(ets, et_t)
                ets = np.insert(ets, index, et_t)
    return ets


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
    polynomial_degree: int
        Degree of the lagrange polynomials that used to interpolate the states.
    ets: np.ndarray of float64
        Array of TDB seconds from J2000 (et dates) of which the data will be taken.
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """

    __slots__ = ["point_id", "polynomial_degree", "ets", "states"]

    def __init__(
        self,
        point_id: int,
        utc_times: List[str],
        body: str,
        lat: float,
        lon: float,
        altitude: float,
        eq_rad: float,
        pol_rad: float,
        source_frame: str,
        target_frame: str,
        polynomial_degree: int = 5,
    ):
        """
        Parameters
        ----------
        point_id: int
            ID code that will be associated with the point.
        utc_times: list of str
            Timestamps at which the data will be calculated, in UTC and a valid SPICE format
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
        source_frame : str
            Name of the frame to transform from.
        target_frame : str
            Name of the frame which the location will be referencing.
        polynomial_degree: int
            Degree of the lagrange polynomials that will be used to interpolate the states.
        """
        self.point_id = point_id
        delta_t = 1  # Arbitrary
        self.polynomial_degree = polynomial_degree
        self.ets = _calculate_ets(utc_times, delta_t, self.polynomial_degree)
        flattening = (eq_rad - pol_rad) / eq_rad
        alt_km = altitude / 1000
        pos_iau = spice.pgrrec(
            body, np.radians(lon), np.radians(lat), alt_km, eq_rad, flattening
        )
        self.states = _calculate_states(
            self.ets, pos_iau, delta_t, source_frame, target_frame
        )


def create_custom_point_kernel(
    obs: Location,
    center: int,
    custom_kernel_dir: str,
    target_frame: str,
) -> None:
    """Creates a SPK custom kernel file containing the data of a point on a body's surface

    Parameters
    ----------
    obs: Location
        Observer point on a body's surface.
    center: int
        NAIF code for center of motion of the body.
    custom_kernel_dir: str
        Path where the writable custom kernel custom.bsp will be stored.
    target_frame : str
        Name of the frame which the location point will be referencing.
    """
    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, "SPK_file", 0)
    spice.spkw09(
        handle,
        obs.point_id,
        center,
        target_frame,
        obs.ets[0],
        obs.ets[-1],
        "0",
        obs.polynomial_degree,
        len(obs.ets),
        obs.states.tolist(),
        obs.ets.tolist(),
    )
    spice.spkcls(handle)


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
