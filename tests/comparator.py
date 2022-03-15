#!/usr/bin/env python3
import spicedmoon as spm
import pylunar
from typing import Tuple
from datetime import datetime

def _decdeg2dms(dd: float) -> Tuple[int, int, int]:
    """
    Converts decimal degrees to degree, minute, second

    Parameters
    ----------
    dd : float
        Value to be transformed from decimal degrees.

    Returns
    -------
    deg : int
        Degrees.
    mnt : int
        Minutes.
    sec : int
        Seconds.
    """
    mnt, sec = divmod(dd * 3600, 60)
    deg, mnt = divmod(mnt, 60)
    return int(deg), int(mnt), int(sec)

def main():
    utc_times = ["2022-01-17 00:00:00", "2022-03-09 14:33:01"]
    kernels_path = "kernels"
    extra_kernels = ["EarthStations.tf", "EarthStations.bsp"]
    extra_kernels_path = "kernels"
    observer_name = "VALLADOLID"
    observer_frame = "VALLADOLID_LOCAL_LEVEL"
    oxf_lat = 51.759
    oxf_lon = -1.2560000
    oxf_alt = 87
    iz_lat = 28.309283
    iz_lon = -16.499143
    lon = -4.70583
    lat = 41.6636
    alt = 705
    frame = "ITRF93"
    correction = True
    moon_datas = spm.spicedmoon.get_moon_datas(lat, lon, alt, utc_times, kernels_path, correction, frame, False)
    moon_datas_extra = spm.spicedmoon.get_moon_datas_from_extra_kernels(utc_times, kernels_path,
        extra_kernels, extra_kernels_path, observer_name, observer_frame, False)

    md_izana = spm.spicedmoon.get_moon_datas(iz_lat, iz_lon, 2400, utc_times, kernels_path, correction, frame, False)
    mde_izana = spm.spicedmoon.get_moon_datas_from_extra_kernels(utc_times, kernels_path,
        extra_kernels, extra_kernels_path, "IZANA", "IZANA_LOCAL_LEVEL", False)

    md_oxf = spm.spicedmoon.get_moon_datas(oxf_lat, oxf_lon, oxf_alt, utc_times, kernels_path, correction, frame, False)
    mde_oxf = spm.spicedmoon.get_moon_datas_from_extra_kernels(utc_times, kernels_path,
        extra_kernels, extra_kernels_path, "OXFORD", "OXFORD_LOCAL_LEVEL", False)
    for i, md in enumerate(moon_datas):
        fecha = utc_times[i]
        mde = moon_datas_extra[i]
        mdi = md_izana[i]
        mdei = mde_izana[i]
        mdo = md_oxf[i]
        mdeo = mde_oxf[i]
        print(fecha)
        print("Valladolid")
        print(md.azimuth, md.zenith, md.mpa_deg)
        print(mde.azimuth, mde.zenith, mde.mpa_deg)
        print(get_pylunar(fecha, lat, lon))
        print("Izaña")
        print(mdi.azimuth, mdi.zenith, mdi.mpa_deg)
        print(mdei.azimuth, mdei.zenith, mdei.mpa_deg)
        print(get_pylunar(fecha, iz_lat, iz_lon))
        print("Oxford")
        print(mdo.azimuth, mdo.zenith, mdo.mpa_deg)
        print(mdeo.azimuth, mdeo.zenith, mdeo.mpa_deg)
        print(get_pylunar(fecha, oxf_lat, oxf_lon))

def get_pylunar(dt_str, lat, lon):
    dt = datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S')
    mi = pylunar.MoonInfo(_decdeg2dms(lat), _decdeg2dms(lon))
    mi.update(dt)
    az = mi.azimuth()
    ze = 90-mi.altitude()
    return az, ze

if __name__ == "__main__":
    main()