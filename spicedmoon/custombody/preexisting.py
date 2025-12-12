"""
Calculate lunar geometries for a  point already existing in an extra custom kernel file.
"""
import os
from datetime import datetime
from typing import List, Union

import spiceypy as spice

from ..basics import dt_to_str, furnsh_safer
from ..custombody.core import get_moon_data_body_ellipsoid
from ..constants import BASIC_KERNELS, MOON_KERNELS
from ..types import MoonData


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
