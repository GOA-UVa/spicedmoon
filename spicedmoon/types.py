from dataclasses import dataclass


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
    lat_sun_rad: float
    lat_obs: float
    lon_obs: float
    mpa_deg: float
    azimuth: float
    zenith: float


@dataclass
class MoonSunData:
    """
    Dataclass with information of the relation between the Sun and the Moon.

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
