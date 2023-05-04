# spicedmoon

![Version 1.0.7](https://img.shields.io/badge/version-1.0.7-informational)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

Calculation of lunar data using NASA’s SPICE toolbox.

This data includes:
- Distance between the Sun and the Moon (in astronomical units)
- Distance between the Sun and the Moon (in kilometers)
- Distance between the Observer and the Moon (in kilometers)
- Selenographic longitude of the Sun (in radians)
- Selenographic latitude of the observer (in degrees)
- Selenographic longitude of the observer (in degrees)
- Moon phase angle (in degrees)
- Azimuth angle (in degrees)
- Zenith angle (in degrees)

It exports the following functions:
* get_moon_datas - Calculates needed MoonData from SPICE toolbox
* get_moon_datas_from_extra_kernels - Calculates needed MoonData from SPICE toolbox
and using data from extra kernels for the observer body

## Requirements

- python>=3.7
- numpy>=1.21.5
- spiceypy>=5.0.0

## Installation

```sh
pip install spicedmoon
```

### Kernels

In order to use the package, a directory with all the kernels must be downloaded.

That directory must contain the following kernels:
- [https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/de421.bsp](https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/de421.bsp)
- [https://naif.jpl.nasa.gov/pub/naif/pds/wgc/kernels/pck/earth_070425_370426_predict.bpc](https://naif.jpl.nasa.gov/pub/naif/pds/wgc/kernels/pck/earth_070425_370426_predict.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/planets/earth_assoc_itrf93.tf)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_080317.tf)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de421_1900-2050.bpc)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0011.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0011.tls)
- [https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc)

## Usage

If one wanted to obtain the lunar azimuth and zenith at Izaña Atmospheric Observatory, on
2022-01-17 at 00:00:00 (UTC), we could use the following snippet of code:

```python
import spicedmoon as spm

lat = 28.309283
lon = -16.499143
alt = 2373
dts = ["2022-01-17 00:00:00"]
mds = spm.get_moon_datas(lat, lon, alt, dts, "./kernels")
az = mds[0].azimuth
ze = mds[0].zenith
```
Note that the kernels directory is called "kernels" and is located in the directory where the
previous code is executed.
