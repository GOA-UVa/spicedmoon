[![Latest Tag][version-shield]][version-url]
[![License: LGPL 3.0][license-shield]][license-url]
[![Stargazers][stars-shield]][stars-url]

<br />
<div align="center">
  <a href="https://github.com/goa-uva/spicedmoon">
    <img src="https://raw.githubusercontent.com/GOA-UVa/spicedmoon/master/docs/images/spicedmoon_logo.png" alt="spicedmoon logo" height="80" style="height: 80px !important;">
  </a>
  <h3 align="center">spicedmoon</h3>
  <p align="center">
    Calculation of lunar data using NASA's SPICE toolkit.
  </p>
</div>


## About the project

`spicedmoon` is a Python toolkit built on top of **SPICE** for computing high-precision lunar observational geometry.  
It provides selenographic coordinates, Sun-Moon-observer geometry, lunar phase angle and sign, azimuth & zenith,
and more, from arbitrary observer locations on Earth or the Moon.

The library offers **two complementary computation methodologies**:

### 1. Custom-Body Method (Custom Kernel Based)

This method reproduces the workflow used in classic C/SPICE pipelines:
a **temporary SPK and frame kernel** is generated for the observer, which becomes a synthetic SPICE "body".
Geometry is then computed using NAIF high-level routines like `subpnt` (INTERCEPT/ELLIPSOID), `subslr`, `spkezr`,
or`pxform`, and local-level frames.

This method uses fully SPICE-native geometry, but needs to generate temporary kernel files for the custom body,
or use pre-existing ones.

### 2. Direct Geometry Method (No Custom Kernels)

A lightweight alternative that computes geometry via **pure vector mathematics**:

SPK ephemerides are don via `spkpos`, and frame transformations via `pxform`, but ellipsoid intersections
are done analytically and lunar phase is computed from vector angles. This methodology
works without temporary observer kernels.

This method yields results very similar to method 1's ones, and fits better in environments where generating
kernel files is undesirable.

## Getting started

### Prerequisites

- python>=3.7

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

### Installation

You can install `spicedmoon` either from the source directory (for development) or directly from PyPI.

From source (editable mode):
```sh
pip install -e .
```

From PyPI (recommended):
```sh
pip install spicedmoon
```

## Usage

`spicedmoon` exposes several high-level functions grouped into two computation
methodologies: the **Direct Geometry Method** (no custom kernels) and the
**Custom-Body Method** (kernel-based). Both approaches compute lunar geometry,
but differ in how observer locations are represented.

Below you will find the recommended functions, organized by methodology.

### 1. Direct Geometry Method (No Custom Kernels)

These functions compute all geometry analytically, without generating observer kernels.
They are lightweight, fast, and ideal for simulation or scripting environments.

#### 1.1 get_moon_datas_llhs

Computes lunar geometry from **planetographic coordinates** (latitude, longitude, height) in degrees and km.

##### Example
```python
from spicedmoon import get_moon_datas_llhs

llhs = [(40.0, -3.5, 0.7)]        # (lat, lon, height_km)
times = ["2025-01-10 00:00:00"]
kernels_path = "/path/to/kernels"

md = get_moon_datas_llhs(llhs, times, kernels_path)[0]

print(md.lon_obs, md.lat_obs) # Observer's selenographic longitude and latitude in degrees
# 1.4414527418916223 -5.396883424855335
print(md.azimuth, md.zenith) # Azimuth and zenith of the target in degrees. Az is measured from North towards East.
# 265.02364998605594 45.359639487448284
print(md.mpa_deg)  # signed moon phase angle in degrees
# -51.31044215492803
```

#### 1.2 get_moon_datas_xyzs

Same as above, but receives **rectangular coordinates** in the chosen source frame.

Useful when your application already works with Earth-fixed or inertial positions.

##### Example
```python
from spicedmoon import get_moon_datas_xyzs

xyz = [(4510.0, 3480.0, 4060.0)]     # km, in source frame (default: J2000)
times = ["2025-01-10 00:00:00"]
kernels_path = "/path/to/kernels"

md = get_moon_datas_xyzs(xyz, times, kernels_path)[0]
print(md.lon_obs, md.lat_obs) # Observer's selenographic longitude and latitude in degrees
# 2.321034649559275 -5.3589043558009095
print(md.lon_sun_rad, md.lat_sun_rad) # Sun's selenographic longitude and latitude in radians
# 0.9199491298659819 -0.02542042276583308
print(md.mpa_deg)  # signed moon phase angle in degrees
# -50.433894317010626
```


### 2. Custom-Body Method (Kernel-Based)

This method creates a temporary SPK + frame kernel for each observer, reproducing
the classic C/SPICE workflow used by NAIF.  
Use this method when absolute consistency with SPICE’s own geometry routines is required
(e.g., comparing with legacy C pipelines).

#### 2.1 get_moon_datas

Computes lunar geometry by:

- Creating a custom observer body,
- Generating an SPK for that body,
- Loading it into SPICE,
- Using high-level routines like `subpnt`, `subslr`, `spkezr`, and `pxform`.

##### Example
```python
from spicedmoon import get_moon_datas

lat = 40.0
lon = -3.5
alt_m = 700     # altitude in meters here
times = ["2025-01-10 00:00:00"]
kernels_path="/path/to/kernels"

md = get_moon_datas(lat, lon, alt_m, times, kernels_path)[0]
print(md.lon_obs, md.lat_obs) # Observer's selenographic longitude and latitude in degrees
# 1.4414527418916219 -5.396888519065468
print(md.azimuth, md.zenith) # Azimuth and zenith of the target in degrees. Az is measured from North towards East.
# 265.02364998605594 45.359639487448284
print(md.mpa_deg)  # signed moon phase angle in degrees
# -51.31044215492803
```

#### 2.2 get_moon_datas_from_extra_kernels

Allows using **pre-existing custom-body kernels** (e.g., Earth station networks).

Useful when you already have `.bsp` and `.tf` files describing station locations.

One must specify the local frame in `observer_frame` to correctly calculate zenith and azimuth.
This frame must be present in the custom extra kernel files.

```python
from spicedmoon import get_moon_datas_from_extra_kernels

times = ["2025-01-10 00:00:00"],
kernels_path = "/path/to/kernels"

md = get_moon_datas_from_extra_kernels(
    times,
    kernels_path,
    extra_kernels=["stations.bsp", "stations.tf"],
    extra_kernels_path="/path/to/extra/kernels",
    observer_name="DSS-14",
    observer_frame="DSS_LOCAL_LEVEL",
)[0]
```

#### 2.3 get_moon_datas_from_moon

Custom-body variant for **observers located on the Moon**.

This function is analogous to `get_moon_datas`, but assumes the body of reference
for planetographic coordinates is the Moon instead of the Earth.
It generates a custom observer body on the lunar surface and uses SPICE to compute
selenographic geometry from that location.

```python
from spicedmoon import get_moon_datas_from_moon

lat = 10
lon = 45
alt_m = 10000
times = ["2025-01-10 00:00:00"]
kernels_path = "/path/to/kernels"

md = get_moon_datas_from_moon(lat, lon, alt_m, times, kernels_path)[0]
print(md.lon_obs, md.lat_obs) # Observer's selenographic longitude and latitude in degrees
# 44.999999999999986 10.000144376768828
print(md.azimuth, md.zenith) # Azimuth and zenith of the target in degrees. Az is measured from North towards East.
# 180.0 19.97647138243788
print(md.mpa_deg)  # signed moon phase angle in degrees
# -13.767196747269022
```

### 3. Lunar + Solar Geometry

#### 3.1 get_sun_moon_datas

This helper computes solar geometry regarding to the moon and returns
a `MoonSunData` structure.

Useful when you need only Sun position in respect to the moon.

```python
from spicedmoon import get_sun_moon_datas

times = ["2025-01-10 00:00:00"]
kernels_path = "/path/to/kernels"

msd = get_sun_moon_datas(times, kernels_path)[0]
print(msd.lon_sun_rad, msd.lat_sun_rad) # Sun's selenographic longitude and latitude in radians
# 0.9199491298659819 -0.025420447501053357
print(msd.dist_sun_moon_km, msd.dist_sun_moon_au) # Distance between the Sun and the Moon in kilometers and in astronomical units
# 147349561.92390758 0.9849709846767327
```

### 4. Data Structures

#### 4.1 `MoonData`

Returned by all `get_moon_datas*` functions (both direct and custom-body).

Fields include:
- `dist_sun_moon_au` : Distance between the Sun and the Moon (in astronomical units)
- `dist_sun_moon_km` : Distance between the Sun and the Moon (in kilometers)
- `dist_obs_moon` : Distance between the Observer and the Moon (in kilometers)
- `lon_sun_rad` : Selenographic longitude of the Sun (in radians)
- `lat_obs` : Selenographic latitude of the observer (in degrees)
- `lon_obs` : Selenographic longitude of the observer (in degrees)
- `mpa_deg` : Moon phase angle (in degrees)
- `azimuth` : Azimuth angle (in degrees)
- `zenith` : Zenith angle (in degrees)

### 4.2 `MoonSunData`

Returned by `get_sun_moon_datas`.

Contains:
- `lon_sun_rad` : Selenographic longitude of the Sun in radians
- `lat_sun_rad` : Selenographic latitude of the Sun in radians
- `dist_sun_moon_km` : Distance between the Sun and the Moon in km
- `dist_sun_moon_au` : Distance between the Sun and the Moon in AU

## Comparing both methods

Under `tests/` one can see mulitple scripts useful to check result differences between methodologies
and against results from external library (`ephem`, `pylunar`, etc.).

Custom-Body and Direct-Geometry should the same results if performed as follows:
```python
from spicedmoon import get_moon_datas, get_moon_datas_llhs
# lat & lon are floats with the geographic coordinates in decimal degrees
# alt is a float with the altitude over sea-level in meters
# dts_str is a str with the date&time in a SPICE-compatible format
# kpath is a str with pointing to the kernels directory path
mds = get_moon_datas(
    lat, lon, alt, dts_str, kpath, earth_as_zenith_observer=False
)
mds = get_moon_datas_llhs(
    [(lat, lon, alt / 1000) for _ in range(len(dts_str))],
    dts_str,
    kpath,
)
```

In the following figure we can see that the relative differences are extremely low.

![Relative Differences: Custom-Body VS Direct-Geometry](https://raw.githubusercontent.com/GOA-UVa/spicedmoon/master/docs/images/custombody_directgeom.png)


## Structure

The package is divided in multiple submodules, each dealing with different calculations
and functionalities. Its structure can be represented in a UML diagram:

![UML diagram](https://raw.githubusercontent.com/GOA-UVa/spicedmoon/master/docs/images/package_structure.png)


## Authors

- [Javier Gatón Herguedas](mailto:gaton@goa.uva.es) - *Maintainer* - [GOA-UVa](https://goa.uva.es)
- [Juan Carlos Antuña Sánchez](mailto:juancarlos.antuna@grasp-earth.com) - *Contributor* - [GRASP-Earth](https://grasp-earth.com)
- [Stefan Adriaensen](mailto:stefan.adriaensen@vito.be) - *Contributor* - [VITO](https://vito.be)
- [Ramiro González Catón](mailto:ramiro@goa.uva.es) - *Contributor* - [GOA-UVa](https://goa.uva.es)
- [Roberto Román Diez](mailto:robertor@goa.uva.es) - *Contributor* - [GOA-UVa](https://goa.uva.es)
- [Carlos Toledano Olmeda](mailto:toledano@goa.uva.es) - *Contributor* - [GOA-UVa](https://goa.uva.es)


## License

Distributed under the LGPL-v3 License. See [LGPL v3](./LICENSE) for more information.


[stars-shield]: https://img.shields.io/github/stars/goa-uva/spicedmoon.svg?style=for-the-badge
[stars-url]: https://github.com/goa-uva/spicedmoon/stargazers
[version-shield]: https://img.shields.io/github/v/tag/goa-uva/spicedmoon?style=for-the-badge
[version-url]: https://github.com/goa-uva/spicedmoon/tags
[license-shield]: https://img.shields.io/github/license/goa-uva/spicedmoon.svg?style=for-the-badge
[license-url]: https://github.com/goa-uva/spicedmoon/blob/master/LICENSE
