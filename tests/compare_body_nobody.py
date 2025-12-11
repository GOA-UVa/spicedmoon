#!/usr/bin/env python3
from typing import List
from datetime import datetime, timedelta
from typing import Tuple
import math

import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("TkAgg")


import spicedmoon as spm


def _decdeg2dms(dd: float) -> Tuple[int, int, int]:
    mnt, sec = divmod(dd * 3600, 60)
    deg, mnt = divmod(mnt, 60)
    return int(deg), int(mnt), int(sec)


def datetime_range(start, end, delta):
    current = start
    while current < end:
        yield current
        current += delta


def get_spicedmoon_earth(dts_str, lat, lon, alt):
    mds = spm.get_moon_datas(
        lat, lon, alt, dts_str, "./kernels", earth_as_zenith_observer=True
    )
    return mds


def get_spicedmoon_llh(dts_str, lat, lon, alt):
    mds = spm.get_moon_datas_llhs(
        [(lat, lon, alt / 1000) for _ in range(len(dts_str))],
        dts_str,
        "./kernels",
        source_planetographic_frame="ITRF93",
    )
    return mds


def get_spicedmoon_obs(dts_str, lat, lon, alt):
    mds = spm.get_moon_datas(
        lat, lon, alt, dts_str, "./kernels", earth_as_zenith_observer=False
    )
    return mds


def get_reldif(a, b, percent=True):
    rd = (a - b) / a
    if percent:
        return 100 * rd
    return rd


def plot_reldifs(mds0: List[spm.MoonData], mds1: List[spm.MoonData]):
    vars = [v for v in mds0[0].__dir__() if not v.startswith("_")]
    mrds = {}
    for v in vars:
        rds = []
        for m0, m1 in zip(mds0, mds1):
            rd = get_reldif(m0.__dict__[v], m1.__dict__[v], False)
            rds.append(rd)
        mrds[v] = rds
    fig, axes = plt.subplots(3, 4)
    for i, k in enumerate(mrds):
        ax = axes[i // 4][i % 4]
        ax.hist(mrds[k], bins=15, color="skyblue", edgecolor="black")
        ax.set_title(k)
    plt.show()


def main():
    dts = [
        dt.strftime("%Y-%m-%d %H:%M:%S")
        for dt in datetime_range(
            datetime(2022, 4, 22, 0), datetime(2022, 4, 22, 6), timedelta(minutes=30)
        )
    ]
    # izana
    lat = 28.309283
    lon = -16.499143
    alt = 2400
    mds0 = get_spicedmoon_obs(dts, lat, lon, alt)
    mds1 = get_spicedmoon_llh(dts, lat, lon, alt)
    matplotlib.rc("axes.formatter", useoffset=False)
    plot_reldifs(mds0, mds1)


if __name__ == "__main__":
    main()
