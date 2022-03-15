#!/usr/bin/env python3
from datetime import datetime, timedelta
from typing import Tuple
import math

import spicedmoon as spm
import pylunar
import ephem

def _decdeg2dms(dd: float) -> Tuple[int, int, int]:
    mnt, sec = divmod(dd * 3600, 60)
    deg, mnt = divmod(mnt, 60)
    return int(deg), int(mnt), int(sec)

def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta

def print_result(az, ze, phase):
    print("{},{},{}".format(az, ze, phase))

def print_pylunar(dts_str, lat, lon, alt):
    mi = pylunar.MoonInfo(_decdeg2dms(lat), _decdeg2dms(lon))
    for dt_s in dts_str:
        dt = datetime.strptime(dt_s, '%Y-%m-%d %H:%M:%S')
        mi.update(dt)
        az = mi.azimuth()
        ze = 90 - mi.altitude()
        ph = mi.libration_phase_angle()
        if ph > 180:
            ph = 360-ph
        print_result(az, ze, ph)

def print_spicedmoon(dts_str, lat, lon, alt):
    mds = spm.get_moon_datas(lat, lon, alt, dts_str, "./kernels")
    for md in mds:
        print_result(md.azimuth, md.zenith, md.mpa_deg)

def print_ephem(dts_str, lat, lon, alt):
    obs = ephem.Observer()
    obs.lat = math.radians(lat)
    obs.long = math.radians(lon)
    m = ephem.Moon()
    for dt_s in dts_str:
        dt = datetime.strptime(dt_s, '%Y-%m-%d %H:%M:%S')
        obs.date = dt
        m.compute(obs)
        az = math.degrees(m.az)
        ze = 90 - math.degrees(m.alt)
        ph = math.degrees(m.moon_phase)
        if ph > 180:
            ph = 360-ph
        print_result(az, ze, ph)


def main():
    dts = [dt.strftime('%Y-%m-%d %H:%M:%S') for dt in 
       datetime_range(datetime(2022, 1, 17, 0), datetime(2022, 1, 17, 23, 31), 
       timedelta(minutes=30))]
    # izana
    lat = 28.309283
    lon = -16.499143
    alt = 2400
    print_spicedmoon(dts, lat, lon, alt)

if __name__ == "__main__":
    main()
