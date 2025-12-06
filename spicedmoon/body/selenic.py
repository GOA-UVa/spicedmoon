import os
from typing import List, Union
from datetime import datetime

import spiceypy as spice

from ..basics import (
    furnsh_safer,
    dt_to_str,
    get_radii_moon,
)
from .core import get_moon_datas_body_ellipsoid_id
from ..constants import MOON_ID_CODE, BASIC_KERNELS, MOON_KERNELS
from ..types import MoonData
from .customkernel import (
    Location,
    create_custom_point_kernel,
    remove_custom_kernel_file,
)


class _MoonLocation(Location):
    """
    Data for the creation of an observer point on Moon's surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Moon's surface
    polynomial_degree: int
        Degree of the lagrange polynomials that used to interpolate the states.
    ets: np.ndarray of float64
        Array of TDB seconds from J2000 (et dates) of which the data will be taken.
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """

    def __init__(
        self,
        point_id: int,
        utc_times: List[str],
        lat: float,
        lon: float,
        altitude: float,
        source_frame: str,
        target_frame: str,
        ignore_bodvrd: bool = True,
    ):
        """
        Parameters
        ----------
        point_id : int
            ID code that will be associated with the point on Moon's surface
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
        ignore_bodvrd : bool
            Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
            1738.1 and 1736
        """
        eq_rad, pol_rad = get_radii_moon(ignore_bodvrd)
        super().__init__(
            point_id,
            utc_times,
            "MOON",
            lat,
            lon,
            altitude,
            eq_rad,
            pol_rad,
            source_frame,
            target_frame,
        )


def _create_moon_point_kernel(
    utc_times: List[str],
    kernels_path: str,
    lat: int,
    lon: int,
    altitude: float,
    id_code: int,
    custom_kernel_dir: str,
    ignore_bodvrd: bool = True,
    source_frame: str = "MOON_ME",
    target_frame: str = "MOON_ME",
) -> None:
    """Creates a SPK custom kernel file containing the data of a point on Earth's surface

    Parameters
    ----------
    utc_times : list of str
        Times at which the lunar data will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    lat : float
        Selenographic latitude (in degrees) of the location.
    lon : float
        Selenographic longitude (in degrees) of the location.
    altitude : float
        Altitude over the sea level in meters.
    id_code : int
        ID code that will be associated with the point on Moon's surface
    custom_kernel_dir: str
        Path where the writable custom kernel custom.bsp will be stored.
    ignore_bodvrd : bool
        Ignore the SPICE function bodvrd for the calculation of the Moon's radii and use the values
        1738.1 and 1736
    source_frame : str
        Name of the frame to transform the coordinates from.
    target_frame : str
        Name of the frame which the location point will be referencing.
    """
    kernels = BASIC_KERNELS + MOON_KERNELS
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        furnsh_safer(k_path)
    obs = _MoonLocation(
        id_code,
        utc_times,
        lat,
        lon,
        altitude,
        source_frame,
        target_frame,
        ignore_bodvrd,
    )
    center = MOON_ID_CODE
    create_custom_point_kernel(obs, center, custom_kernel_dir, target_frame)
    spice.kclear()


def get_moon_datas_from_moon(
    lat: float,
    lon: float,
    altitude: float,
    times: Union[List[str], List[datetime]],
    kernels_path: str,
    correct_zenith_azimuth: bool = True,
    observer_frame: str = "MOON_ME",
    custom_kernel_path: str = None,
    ignore_bodvrd: bool = True,
    source_frame: str = "MOON_ME",
    target_frame: str = "MOON_ME",
) -> List[MoonData]:
    """Calculation of needed Moon data from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    lat : float
        Selenographic latitude (in degrees) of the location.
    lon : float
        Selenographic longitude (in degrees) of the location.
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
    id_code = MOON_ID_CODE * 1000 + 100
    utc_times = dt_to_str(times)
    if len(utc_times) == 0:
        return []
    remove_custom_kernel_file(custom_kernel_path)
    _create_moon_point_kernel(
        utc_times,
        kernels_path,
        lat,
        lon,
        altitude,
        id_code,
        custom_kernel_path,
        ignore_bodvrd,
        source_frame,
        target_frame,
    )
    return get_moon_datas_body_ellipsoid_id(
        utc_times,
        kernels_path,
        id_code,
        observer_frame,
        custom_kernel_path,
        correct_zenith_azimuth,
        lat,
        lon,
        False,
        ignore_bodvrd,
    )
