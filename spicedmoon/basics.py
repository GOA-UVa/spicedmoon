"""
Common basic functions that help and improve in SPICE usage.
"""
import os
import time
from typing import List, Union, Tuple
from datetime import datetime, timezone
import warnings

import spiceypy as spice

from .constants import MOON_EQ_RAD, MOON_POL_RAD


def furnsh_safer(k_path: str):
    """
    Performs SPICE's `furnsh_c`, but in case that it fails it tries again after a small time
    interval. Furnsh very rarely crashes, but it can be solved trying again.

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


def _is_dt_tz_aware(dt: datetime) -> bool:
    """Checks if a datetime is timezone aware or not

    Parameters
    ----------
    dt: datetime
        Datetime to check timezone-awareness

    Result
    ------
    is_aware: bool
        Timezone-awareness of the datetime
    """
    return dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None


def dt_to_str(dts: Union[List[datetime], List[str]]) -> List[str]:
    """Convert a list of datetimes into a list of string dates in a valid SPICE `str` format.

    Parameters
    ----------
    dts: list of datetimes | list of str
        List of datetimes that will be converted to `str`. They must be timezone aware.
        A list of already `str` can be given instead, and it will be returned without
        change.

    Returns
    -------
    utc_times: list of str
        List of the timestamps in a valid `str` format for SPICE.
    """
    utc_times = []
    for dt in dts:
        if isinstance(dt, datetime):
            if not _is_dt_tz_aware(dt):
                warnings.warn("Using timezone-naive datetime object", RuntimeWarning)
            dt_utc = dt.astimezone(timezone.utc)
            utc_times.append(dt_utc.strftime("%Y-%m-%d %H:%M:%S"))
        else:
            utc_times.append(dt)
    return utc_times


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
    eq_rad, pol_rad = MOON_EQ_RAD, MOON_POL_RAD
    if not ignore_bodvrd:
        # The ones obtained with bodvrd are not correct, not accurate
        _, radii_moon = spice.bodvrd("MOON", "RADII", 3)
        eq_rad, pol_rad = radii_moon[0], radii_moon[2]
    return eq_rad, pol_rad
