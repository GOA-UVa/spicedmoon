import os
from typing import List, Union
from datetime import datetime

import numpy as np
import spiceypy as spice

from .types import MoonData
from .basics import Location, furnsh_safer, dt_to_str, remove_custom_kernel_file
from .core import get_moon_datas_id
from .constants import BASIC_KERNELS, MOON_KERNELS, CUSTOM_KERNEL_NAME

EARTH_ID_CODE = 399


class _EarthLocation(Location):
    """
    Data for the creation of an observer point on Earth's surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Earth's surface
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """

    def __init__(
        self,
        point_id: int,
        lat: float,
        lon: float,
        altitude: float,
        ets: np.ndarray,
        delta_t: float,
        source_frame: str,
        target_frame: str,
    ):
        """
        Parameters
        ----------
        point_id : int
            ID code that will be associated with the point on Earth's surface
        lat : float
            Geographic latitude of the observer point
        lon : float
            Geographic longitude of the observer point
        altitude : float
            Altitude over the sea level in meters.
        ets : np.ndarray
            Array of TDB seconds from J2000 (et dates) of which the data will be taken
        delta_t : float
            TDB seconds between states
        source_frame : str
            Name of the frame to transform from.
        target_frame : str
            Name of the frame which the location will be referencing.
        """
        eq_rad = 6378.1366  # Earth equatorial radius
        pol_rad = 6356.7519  # Earth polar radius
        super().__init__(
            point_id,
            "EARTH",
            lat,
            lon,
            altitude,
            eq_rad,
            pol_rad,
            ets,
            delta_t,
            source_frame,
            target_frame,
        )


def _create_earth_point_kernel(
    utc_times: List[str],
    kernels_path: str,
    lat: int,
    lon: int,
    altitude: float,
    id_code: int,
    custom_kernel_dir: str,
    source_frame: str = "ITRF93",
    target_frame: str = "ITRF93",
) -> None:
    """Creates a SPK custom kernel file containing the data of a point on Earth's surface

    Parameters
    ----------
    utc_times : list of str
        Times at which the lunar data will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    id_code : int
        ID code that will be associated with the point on Earth's surface
    custom_kernel_dir: str
        Path where the writable custom kernel custom.bsp will be stored.
    source_frame : str
        Name of the frame to transform the coordinates from.
    target_frame : str
        Name of the frame which the location point will be referencing.
    """
    kernels = BASIC_KERNELS
    if "MOON" in source_frame or "MOON" in target_frame:
        kernels += MOON_KERNELS
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)

    polynomial_degree = 5
    # Degree of the lagrange polynomials that will be used to interpolate the states
    delta_t = 1  # TDB seconds between states. Arbitrary.
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

    obs = _EarthLocation(
        id_code, lat, lon, altitude, ets, delta_t, source_frame, target_frame
    )

    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, "SPK_file", 0)

    center = EARTH_ID_CODE
    spice.spkw09(
        handle,
        obs.point_id,
        center,
        target_frame,
        ets[0],
        ets[-1],
        "0",
        polynomial_degree,
        len(ets),
        obs.states.tolist(),
        ets.tolist(),
    )
    spice.spkcls(handle)
    spice.kclear()


def get_moon_datas(
    lat: float,
    lon: float,
    altitude: float,
    times: Union[List[str], List[datetime]],
    kernels_path: str,
    correct_zenith_azimuth: bool = True,
    observer_frame: str = "ITRF93",
    earth_as_zenith_observer: bool = False,
    custom_kernel_path: str = None,
    ignore_bodvrd: bool = True,
    source_frame: str = "ITRF93",
    target_frame: str = "ITRF93",
) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    lat : float
        Geographic latitude (in degrees) of the location.
    lon : float
        Geographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    times : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    correct_zenith_azimuth : bool
        In case that it's calculated without using the extra kernels, the coordinates should be
        corrected rotating them into the correct location.
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
    earth_as_zenith_observer : bool
        If True the Earth will be used as the observer for the zenith and azimuth calculation.
        Otherwise it will be the actual observer. By default is False.
    custom_kernel_path: str
        Path of the kernel custom.bsp that will be edited by the library, not only read.
        If none, it will be the same as kernels_path.
    ignore_bodvrd : bool
        Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
        1738.1 and 1736
    source_frame : str
        Name of the frame to transform the coordinates from.
    target_frame : str
        Name of the frame which the location point will be referencing.
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    if custom_kernel_path == None:
        custom_kernel_path = kernels_path
    id_code = 399100
    utc_times = dt_to_str(times)
    if len(utc_times) == 0:
        return []
    remove_custom_kernel_file(custom_kernel_path)
    _create_earth_point_kernel(
        utc_times,
        kernels_path,
        lat,
        lon,
        altitude,
        id_code,
        custom_kernel_path,
        source_frame,
        target_frame,
    )
    return get_moon_datas_id(
        utc_times,
        kernels_path,
        id_code,
        observer_frame,
        custom_kernel_path,
        correct_zenith_azimuth,
        lat,
        lon,
        earth_as_zenith_observer,
        ignore_bodvrd,
    )
