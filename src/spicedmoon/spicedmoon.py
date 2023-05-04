"""Spiced Moon

Calculation of lunar data using NASA's SPICE toolbox.

It exports the following functions:

    * get_moon_datas - Calculates needed MoonData from SPICE toolbox
    * get_moon_datas_from_extra_kernels - Calculates needed MoonData from SPICE toolbox
        and using data from extra kernels for the observer body
    * get_sun_moon_datas - Calculates solar selenographic coordinates.
    * get_moon_datas_from_moon - Calculates needed MoonData from SPICE toolbox from selenographic coordinates
"""

"""___Built-In Modules___"""
from dataclasses import dataclass
import os
import math
from typing import List, Union
from datetime import datetime, timezone
import time

"""___Third-Party Modules___"""
import numpy as np
import spiceypy as spice

"""___spiced_moon Modules___"""
# import here

"""___Authorship___"""
__author__ = 'Javier Gatón Herguedas, Juan Carlos Antuña Sánchez, Ramiro González Catón,\
Roberto Román, Carlos Toledano'
__created__ = "2022/03/03"
__maintainer__ = "Javier Gatón Herguedas"
__email__ = "gaton@goa.uva.es"
__status__ = "Development"

CUSTOM_KERNEL_NAME = "custom.bsp"
EARTH_ID_CODE = 399
MOON_ID_CODE = 301

_DEFAULT_OBSERVER_NAME = "Observer"
_DEFAULT_OBSERVER_FRAME = "Observer_LOCAL_LEVEL"
_DEFAULT_OBSERVER_ZENITH_NAME = "EARTH"

@dataclass
class MoonData:
    """
    Moon data needed to calculate Moon's irradiance, probably obtained from NASA's SPICE Toolbox

    Attributes
    ----------
    dist_sun_moon_au : float
        Distance between the Sun and the Moon (in astronomical units)
    dist_sun_moon_km : float
        Distance between the Sun and the Moon (in kilometers)
    dist_obs_moon : float
        Distance between the Observer and the Moon (in kilometers)
    lon_sun_rad : float
        Selenographic longitude of the Sun (in radians)
    lat_obs : float
        Selenographic latitude of the observer (in degrees)
    lon_obs : float
        Selenographic longitude of the observer (in degrees)
    mpa_deg : float
        Moon phase angle (in degrees)
    azimuth : float
        Azimuth angle (in degrees)
    zenith : float
        Zenith angle (in degrees)

    """
    dist_sun_moon_au: float
    dist_sun_moon_km: float
    dist_obs_moon: float
    lon_sun_rad: float
    lat_obs: float
    lon_obs: float
    mpa_deg: float
    azimuth: float
    zenith: float

@dataclass
class MoonSunData:
    """Dataclass with information of the relation between the Sun and the Moon.

    Attributes
    ----------
    lon_sun_rad: float
        Selenographic longitude of the Sun in radians
    lat_sun_rad: float
        Selenographic latitude of the Sun in radians
    dist_sun_moon_km: float
        Distance between the Sun and the Moon in km
    dist_sun_moon_au: float
        Distance between the Sun and the Moon in AU
    """
    lon_sun_rad: float
    lat_sun_rad: float
    dist_sun_moon_km: float
    dist_sun_moon_au: float

def _furnsh_safer(k_path: str):
    """
    Perfroms the SPICE furnsh_c function, but in case that it fails it tries again after a small time
    interval.

    Furnsh has produced an error once, very rarely, but it can be solved trying again, so this is what this
    function solves.

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

def _calculate_states(ets: np.ndarray, pos_iau: np.ndarray, delta_t: float,
                      source_frame: str, target_frame: str) -> np.ndarray:
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
            spice.pxform(source_frame, target_frame, et_value),
            pos_iau)

    for i in range(len(ets) - 1):
        states[i, num_coordinates:] = (states[i + 1, :num_coordinates] -
                                       states[i, :num_coordinates]) / delta_t

    pos_np1 = np.dot(
        spice.pxform(source_frame, target_frame, ets[-1] + delta_t),
        pos_iau)
    states[-1, num_coordinates:] = (pos_np1 - states[-1, :num_coordinates]) / delta_t
    return states

@dataclass
class _EarthLocation():
    """
    Data for the creation of an observer point on earth surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Earth's surface
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """
    __slots__ = ['point_id', 'states']
    def __init__(self, point_id: int, lat: float, lon: float, altitude: float, ets: np.ndarray,
                 delta_t: float, source_frame: str, target_frame: str):
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
        self.point_id = point_id
        eq_rad = 6378.1366 # Earth equatorial radius
        pol_rad = 6356.7519 # Earth polar radius
        alt_km = altitude/1000
        flattening = (eq_rad - pol_rad)/eq_rad
        pos_iau_earth = spice.pgrrec('EARTH', math.radians(lon), math.radians(lat), alt_km,
                                     eq_rad, flattening)
        self.states = _calculate_states(ets, pos_iau_earth, delta_t, source_frame, target_frame)


@dataclass
class _MoonLocation():
    """
    Data for the creation of an observer point on Moon's surface

    Attributes
    ----------
    point_id : int
        ID code that will be associated with the point on Moon's surface
    states : np.ndarray of float64
        Array of geometric states of body relative to center
    """
    __slots__ = ['point_id', 'states']
    def __init__(self, point_id: int, lat: float, lon: float, altitude: float, ets: np.ndarray,
                 delta_t: float, source_frame: str, target_frame: str):
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
        """
        self.point_id = point_id
        eq_rad = 1738.1 # Moon equatorial radius
        pol_rad = 1736 # Moon polar radius
        alt_km = altitude/1000
        flattening = (eq_rad - pol_rad)/eq_rad
        pos_iau_moon = spice.pgrrec('MOON', math.radians(lon), math.radians(lat), alt_km,
                                     eq_rad, flattening)
        self.states = _calculate_states(ets, pos_iau_moon, delta_t, source_frame, target_frame)


def _get_sun_moon_data(utc_time: str):
    et_date = spice.str2et(utc_time)
    m_eq_rad = 1738.1
    m_pol_rad = 1736
    flattening = (m_eq_rad-m_pol_rad)/m_eq_rad
    # Calculate selenographic longitude of sun
    sun_spoint, _, _ = spice.subslr("INTERCEPT/ELLIPSOID", "MOON", et_date, 'MOON_ME',
                                    "NONE", "SUN")
    lon_sun_rad, lat_sun_rad, _ = spice.recpgr("MOON", sun_spoint, m_eq_rad, flattening)

    # Calculate the distance between sun and moon (AU)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", "SUN")
    dist_sun_moon_km = math.sqrt(state[0]**2 + state[1]**2 + state[2]**2)
    dist_sun_moon_au = spice.convrt(dist_sun_moon_km, "KM", "AU")




    limit_lat_rad = math.pi/2
    limit_lon_rad = math.pi
    if lat_sun_rad > limit_lat_rad:
        lat_sun_rad = limit_lat_rad + (limit_lat_rad-lat_sun_rad)
        lon_sun_rad -= limit_lon_rad
    elif lat_sun_rad < -limit_lat_rad:
        lat_sun_rad = -limit_lat_rad - (limit_lat_rad+lat_sun_rad)
        lon_sun_rad += limit_lon_rad

    while lon_sun_rad > limit_lon_rad:
        lon_sun_rad -= limit_lon_rad*2
    while lon_sun_rad < -limit_lon_rad:
        lon_sun_rad += limit_lon_rad*2

    return MoonSunData(lon_sun_rad, lat_sun_rad, dist_sun_moon_km, dist_sun_moon_au)


def _get_moon_data(utc_time: str, observer_name: str = _DEFAULT_OBSERVER_NAME,
                   observer_frame: str = _DEFAULT_OBSERVER_FRAME,
                   observer_zenith_name: str = _DEFAULT_OBSERVER_ZENITH_NAME,
                   correct_zenith_azimuth: bool = False, longitude: float = 0,
                   colat: float = 0) -> MoonData:
    """Calculation of the moon data for the given utc_time for the loaded observer

    Parameters
    ----------
    utc_time : str
        Time at which the lunar data will be calculated, in a valid UTC DateTime format
    observer_name : str
        Name of the body of the observer that should be loaded from the extra kernels.
        By default is "Observer", in which case it shouldn't be loaded from the extra
        kernels but from the custom kernel.
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
    observer_zenith_name : str
        The observer used for the zenith and azimuth calculation. By default it's "EARTH".
    correct_zenith_azimuth : bool
        In case that it's calculated without using the extra kernels, the coordinates should be
        corrected rotating them into the correct location.
    longitude : float
        Geographic longitude of the observer point. Used if it's calculated without using the
        extra kernels.
    colat : float
        Geographic colatitude of the observer point. Used if it's calculated without using the
        extra kernels.
    Returns
    -------
    MoonData
        Moon data obtained from SPICE toolbox
    """
    et_date = spice.str2et(utc_time)

    #_, radios_luna = spice.bodvrd("MOON", "RADII", 3)
    # The ones obtained with bodvrd are not correct
    m_eq_rad = 1738.1 #radios_luna[0] # Moon Equatorial Radius
    m_pol_rad = 1736 #radios_luna[2] # Moon polar radius
    flattening = (m_eq_rad-m_pol_rad)/m_eq_rad

    # Calculate moon zenith and azimuth
    state_zenith, _ = spice.spkezr("MOON", et_date, observer_frame, "NONE", observer_zenith_name)
    rectan_zenith = np.split(state_zenith, 2)[0]
    if correct_zenith_azimuth:
        lon_rad = (longitude+180) * spice.rpd()
        colat_rad = colat * spice.rpd()
        bf2tp = spice.eul2m(-lon_rad,-colat_rad, 0, 3,2,3)
        rectan_zenith = spice.mtxv(bf2tp, rectan_zenith)

    _, longi, lati = spice.reclat(rectan_zenith)

    zenith = 90.0 - lati * spice.dpr()
    azimuth = 180 - longi * spice.dpr()

    # Calculate moon phase angle
    spoint, _, _ = spice.subpnt("INTERCEPT/ELLIPSOID", "MOON", et_date, 'MOON_ME',
                                "NONE", observer_name)
    phase = spice.phaseq(et_date, "MOON", "SUN", observer_name, "NONE")
    phase = phase * spice.dpr()

    et_date_2 = et_date + 1
    phase2 = spice.phaseq(et_date_2, "MOON", "SUN", observer_name, "NONE")
    phase2 = phase2 * spice.dpr()
    if phase2 < phase:
        phase = -phase

    # Calculate selenographic coordinates of the observer
    lon_obs, lat_obs, _ = spice.recpgr("MOON", spoint, m_eq_rad, flattening)
    lon_obs = lon_obs * spice.dpr()
    lat_obs = lat_obs * spice.dpr()

    # Calculate the distance between observer and moon (KM)
    state, _ = spice.spkezr("MOON", et_date, "MOON_ME", "NONE", observer_name)
    dist_obs_moon = math.sqrt(state[0]**2 + state[1]**2 + state[2]**2)

    smd = _get_sun_moon_data(utc_time)
    lon_sun_rad = smd.lon_sun_rad
    dist_sun_moon_km = smd.dist_sun_moon_km
    dist_sun_moon_au = smd.dist_sun_moon_au

    limit_lat = 90
    limit_lon = 180
    if lat_obs > limit_lat:
        lat_obs = limit_lat + (limit_lat-lat_obs)
        lon_obs -= limit_lon
    elif lat_obs < -limit_lat:
        lat_obs = -limit_lat - (limit_lat+lat_obs)
        lon_obs += limit_lon

    while lon_obs > limit_lon:
        lon_obs -= limit_lon*2
    while lon_obs < -limit_lon:
        lon_obs += limit_lon*2

    moon_data = MoonData(dist_sun_moon_au, dist_sun_moon_km, dist_obs_moon, lon_sun_rad,
                         lat_obs, lon_obs, phase, azimuth, zenith)
    return moon_data


def _get_moon_datas_id(utc_times: List[str], kernels_path: str,
                       observer_id: int, observer_frame: str,
                       custom_kernel_dir: str,
                       correct_zenith_azimuth: bool = False,
                       latitude: float = 0, longitude: float = 0,
                       earth_as_zenith_observer: bool = False) -> List[MoonData]:
    """Calculation of needed MoonDatas from SPICE toolbox

    Moon phase angle, selenographic coordinates and distance from observer point to moon.
    Selenographic longitude and distance from sun to moon.

    Parameters
    ----------
    utc_times : list of str
        Times at which the lunar data will be calculated, in a valid UTC DateTime format
    kernels_path : str
        Path where the SPICE kernels are stored
    observer_id : int
        Observer's body ID
    observer_frame : str
        Observer frame that will be used in the calculations of the azimuth and zenith.
    custom_kernel_dir: str
        Path where the writable kernel custom.bsp will be stored.
    correct_zenith_azimuth : bool
        In case that it's calculated without using the extra kernels, the coordinates should be
        corrected rotating them into the correct location.
    latitude : float
        Geographic latitude of the observer point.
    longitude : float
        Geographic longitude of the observer point.
    earth_as_zenith_observer : bool
        If True the Earth will be used as the observer for the zenith and azimuth calculation.
        Otherwise it will be the actual observer. By default is False.
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf",
               "pck00010.tpc", "naif0011.tls", "de421.bsp", "earth_assoc_itrf93.tf",
               "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]

    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        _furnsh_safer(k_path)
    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    _furnsh_safer(custom_kernel_path)

    observer_name = _DEFAULT_OBSERVER_NAME
    spice.boddef(observer_name, observer_id)
    if earth_as_zenith_observer:
        zenith_observer = "EARTH"
    else:
        zenith_observer = observer_name
    moon_datas = []
    colat = 90-(latitude%90)
    lon = longitude%180
    for utc_time in utc_times:
        new_md = _get_moon_data(utc_time, observer_name, observer_frame, zenith_observer,
            correct_zenith_azimuth, lon, colat)
        moon_datas.append(new_md)

    spice.kclear()

    return moon_datas


def _create_earth_point_kernel(utc_times: List[str], kernels_path: str, lat: int, lon: int,
                               altitude: float, id_code: int, custom_kernel_dir: str) -> None:
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
    """
    kernels = ["pck00010.tpc", "naif0011.tls", "earth_assoc_itrf93.tf",
               "de421.bsp", "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        _furnsh_safer(k_path)

    polynomial_degree = 5
    # Degree of the lagrange polynomials that will be used to interpolate the states
    delta_t = 1 # TDB seconds between states. Arbitrary.
    min_states_polynomial = polynomial_degree + 1
    # Min # states that are required to define a polynomial of that degree
    ets = np.array([])
    left_states = int(min_states_polynomial/2)
    right_states = left_states + min_states_polynomial%2
    for utc_time in utc_times:
        et0 = spice.str2et(utc_time)
        etprev = et0 - delta_t * left_states
        etf = et0 + delta_t * right_states
        ets_t = np.arange(etprev, etf, delta_t)
        for et_t in ets_t:
            if et_t not in ets:
                index = np.searchsorted(ets, et_t)
                ets = np.insert(ets, index, et_t)

    target_frame = source_frame = 'ITRF93'
    obs = _EarthLocation(id_code, lat, lon, altitude, ets, delta_t, source_frame, target_frame)

    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, 'SPK_file', 0)

    center = EARTH_ID_CODE
    spice.spkw09(handle, obs.point_id, center, target_frame,
                 ets[0], ets[-1], '0', polynomial_degree, len(ets),
                 obs.states.tolist(), ets.tolist())
    spice.spkcls(handle)

    spice.kclear()


def _create_moon_point_kernel(utc_times: List[str], kernels_path: str, lat: int, lon: int,
                               altitude: float, id_code: int, custom_kernel_dir: str) -> None:
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
    """
    kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf",
               "pck00010.tpc", "naif0011.tls", "earth_assoc_itrf93.tf",
               "de421.bsp", "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        _furnsh_safer(k_path)

    polynomial_degree = 5
    # Degree of the lagrange polynomials that will be used to interpolate the states
    delta_t = 1 # TDB seconds between states. Arbitrary.
    min_states_polynomial = polynomial_degree + 1
    # Min # states that are required to define a polynomial of that degree
    ets = np.array([])
    left_states = int(min_states_polynomial/2)
    right_states = left_states + min_states_polynomial%2
    for utc_time in utc_times:
        et0 = spice.str2et(utc_time)
        etprev = et0 - delta_t * left_states
        etf = et0 + delta_t * right_states
        ets_t = np.arange(etprev, etf, delta_t)
        for et_t in ets_t:
            if et_t not in ets:
                index = np.searchsorted(ets, et_t)
                ets = np.insert(ets, index, et_t)

    target_frame = source_frame = 'MOON_ME'
    obs = _MoonLocation(id_code, lat, lon, altitude, ets, delta_t, source_frame, target_frame)

    custom_kernel_path = os.path.join(custom_kernel_dir, CUSTOM_KERNEL_NAME)
    handle = spice.spkopn(custom_kernel_path, 'SPK_file', 0)

    center = MOON_ID_CODE
    spice.spkw09(handle, obs.point_id, center, target_frame,
                 ets[0], ets[-1], '0', polynomial_degree, len(ets),
                 obs.states.tolist(), ets.tolist())
    spice.spkcls(handle)

    spice.kclear()


def _remove_custom_kernel_file(kernels_path: str) -> None:
    """Remove the custom SPK kernel file if it exists

    Parameters
    ----------
    kernels_path : str
        Path where the SPICE kernels are stored
    """
    custom_kernel_path = os.path.join(kernels_path, CUSTOM_KERNEL_NAME)
    if os.path.exists(custom_kernel_path):
        os.remove(custom_kernel_path)


def _dt_to_str(dts: Union[List[datetime], List[str]]) -> List[str]:
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


def get_sun_moon_datas(times: Union[List[str], List[datetime]], kernels_path: str) -> List[MoonSunData]:
    """
    Obtain solar selenographic coordinates of at multiple times.

    times : list of str | list of datetime
        Times at which the lunar data will be calculated.
        If they are str, they must be in a valid UTC format allowed by SPICE, such as
        %Y-%m-%d %H:%M:%S.
        If they are datetimes they must be timezone aware, or they will be understood
        as computer local time.
    kernels_path : str
        Path where the SPICE kernels are stored
    """
    utc_times = _dt_to_str(times)
    if(len(utc_times) == 0):
        return []
    kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf",
               "pck00010.tpc", "naif0011.tls", "de421.bsp", "earth_assoc_itrf93.tf",
               "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]
    for kernel in kernels:
        k_path = os.path.join(kernels_path, kernel)
        _furnsh_safer(k_path)
    
    msds = []
    for utc_time in utc_times:
        msd = _get_sun_moon_data(utc_time)
        msds.append(msd)
    spice.kclear()
    return msds

def get_moon_datas_from_extra_kernels(times: Union[List[str], List[datetime]], kernels_path: str,
                                      extra_kernels: List[str], extra_kernels_path: str,
                                      observer_name: str, observer_frame: str,
                                      earth_as_zenith_observer: bool = False
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
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    base_kernels = ["moon_pa_de421_1900-2050.bpc", "moon_080317.tf",
               "pck00010.tpc", "naif0011.tls", "de421.bsp", "earth_assoc_itrf93.tf",
               "earth_latest_high_prec.bpc", "earth_070425_370426_predict.bpc"]
    for kernel in base_kernels:
        k_path = os.path.join(kernels_path, kernel)
        _furnsh_safer(k_path)
    for kernel in extra_kernels:
        k_path = os.path.join(extra_kernels_path, kernel)
        _furnsh_safer(k_path)

    if earth_as_zenith_observer:
        zenith_observer = "EARTH"
    else:
        zenith_observer = observer_name
    moon_datas = []
    utc_times = _dt_to_str(times)
    for utc_time in utc_times:
        moon_datas.append(_get_moon_data(utc_time, observer_name, observer_frame,
            zenith_observer))

    spice.kclear()

    return moon_datas


def get_moon_datas(lat: float, lon: float, altitude: float,
                   times: Union[List[str], List[datetime]],
                   kernels_path: str, correct_zenith_azimuth: bool = True,
                   observer_frame: str = "ITRF93",
                   earth_as_zenith_observer: bool = False, custom_kernel_path: str = None,
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
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    if custom_kernel_path == None:
        custom_kernel_path = kernels_path
    id_code = 399100
    utc_times = _dt_to_str(times)
    if(len(utc_times) == 0):
        return []
    _remove_custom_kernel_file(custom_kernel_path)
    _create_earth_point_kernel(utc_times, kernels_path, lat, lon, altitude, id_code, custom_kernel_path)
    return _get_moon_datas_id(utc_times, kernels_path, id_code, observer_frame, custom_kernel_path,
        correct_zenith_azimuth, lat, lon, earth_as_zenith_observer)


def get_moon_datas_from_moon(lat: float, lon: float, altitude: float,
                   times: Union[List[str], List[datetime]],
                   kernels_path: str, correct_zenith_azimuth: bool = True,
                   observer_frame: str = "MOON_ME", custom_kernel_path: str = None,
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
    custom_kernel_path: str
        Path of the kernel custom.bsp that will be edited by the library, not only read.
        If none, it will be the same as kernels_path.
    Returns
    -------
    list of MoonData
        Moon data obtained from SPICE toolbox
    """
    if custom_kernel_path == None:
        custom_kernel_path = kernels_path
    id_code = 301100
    utc_times = _dt_to_str(times)
    if(len(utc_times) == 0):
        return []
    _remove_custom_kernel_file(custom_kernel_path)
    _create_moon_point_kernel(utc_times, kernels_path, lat, lon, altitude, id_code, custom_kernel_path)
    return _get_moon_datas_id(utc_times, kernels_path, id_code, observer_frame, custom_kernel_path,
        correct_zenith_azimuth, lat, lon, False)
