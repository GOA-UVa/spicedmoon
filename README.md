# spicedmoon

![Version 0.1.0](https://img.shields.io/badge/version-0.1.0-informational)

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

- numpy>=1.22.2
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