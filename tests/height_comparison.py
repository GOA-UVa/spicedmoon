#!/usr/bin/env python3

import spicedmoon as spm

def main():
    utc_times = ["2022-01-17 01:00:00"]
    kernels_path = "kernels"
    iz_lat = 28.309283
    iz_lon = -16.499143
    alts = [0, 10, 100, 1000, 10000, 100000, 1000000]
    for alt in alts:
        md = spm.get_moon_datas(iz_lat, iz_lon, alt, utc_times, kernels_path)[0]
        print("Altitude: {}".format(alt))
        print("Azimuth: {}, Zenith: {}, Phase: {}".format(md.azimuth, md.zenith, md.mpa_deg))

if __name__ == "__main__":
    main()
