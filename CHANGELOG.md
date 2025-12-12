# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # "## [unreleased] - yyyy-mm-dd"

## [1.1.0] - yyyy-mm-dd

### Added
- New modular package layout:
  - Core public functionality is now split into dedicated submodules: `geometry`, `basics`, `angular`,
    `coordinates`, `constants`, `types`, `heliac`, and the `custombody` package (`custombody.core`,
    `custombody.preexisting`, `custombody.geotic`, `custombody.selenic`).
  - This clarifies the separation between *direct geometry* and *custom-body (kernel-based)* methods.
- New top-level *direct geometry* function `get_moon_datas_llhs(...)`:
  - Computes lunar geometry directly from planetographic coordinates (latitude, longitude, height)
    without generating custom observer kernels.
- Extended *direct geometry* to support zenith and azimuth: Direct-geometry computations can now obtain zenith and azimuth values.
- New **intercept-ellipsoid** option for direct geometry:
  - Direct-geometry methods (`get_moon_datas_xyzs(...)`, `get_moon_datas_llhs(...)`) can now compute the
    intersection with the Moon surface using SPICE’s **ellipsoid intercept** formulation.
  - This yields results that closely match the custom-body approach (`get_moon_datas_from_moon(...)` /
    `get_moon_datas_from_extra_kernels(...)`) while keeping the simplicity and performance of direct geometry.
- New `spicedmoon.angular` module:
  - `get_zn_az(...)` computes zenith and azimuth from SEZ / rectangular coordinates.
  - `get_phase_sign(...)` encapsulates the sign logic for the phase angle.
- Documentation & tooling:
  - Added `docs/` with diagrams (e.g. package structure UML and images) and improved usage documentation.
  - Added `tests/` with comparison and regression scripts for lunar geometry (e.g. `comparator.py`, `night_calculator.py`, etc.).
  - Added configuration for code quality tools (Black, Ruff, mypy, etc.) and pre-commit hooks.

### Changed
- **License**:
  - Project license changed from **MIT** to **LGPL-3.0-only**.
- **Project configuration**:
  - Migrated project metadata from `setup.cfg` to a modern `pyproject.toml` configuration.
  - Version is now managed dynamically via `setuptools_scm`.
  - Updated project classifiers and metadata (scientific/astronomy topics, explicit Python version classifiers up to 3.14).
- **Public API organization**:
  - High-level functions (`get_moon_datas`, `get_moon_datas_from_extra_kernels`, `get_moon_datas_from_moon`,
    `get_sun_moon_datas`, `get_moon_datas_xyzs`, `get_moon_datas_llhs`) are now thin wrappers around the new modular
    internals in `geometry` and `custombody.*`.
  - `MoonData` and `MoonSunData` are now defined as `dataclass`es in `spicedmoon.types` and imported from there at package level.
- **Moon radii handling**:
  - Moon equatorial and polar radii are now obtained by default from internal constants (`constants.MOON_EQ_RAD`, `constants.MOON_POL_RAD`)
    via `basics.get_radii_moon(...)`.
  - Optionally, SPICE `bodvrd("MOON", "RADII", 3)` can still be used by toggling the `ignore_bodvrd` flag, but is no
    longer the default due to its lower accuracy for this use case.
- **Datetime handling**:
  - The internal datetime-to-string conversion has been promoted to the public `basics.dt_to_str(...)` helper and
    now emits a `RuntimeWarning` if any `datetime` is timezone-naive.
  - This makes timezone assumptions explicit instead of silently treating naive datetimes as local time, preventing subtle bugs.
- **Direct geometry behaviour**:
  - Direct-geometry functions now use **intercept-ellipsoid as the default** method for determining
    the intersection with the lunar reference ellipsoid.
- **Backwards compatibility & deprecation**:
  - Importing `spicedmoon.spicedmoon` is still supported, but now goes through the `spicedmoon.deprecated` shim.
  - Accessing `spicedmoon.spicedmoon` emits a `FutureWarning`, and users are encouraged to import directly
    from `spicedmoon` or the new submodules (`geometry`, `custombody.*`, etc.).
- **General cleanup**:
  - Improved type hints across the codebase.
  - Clarified docstrings and parameter names for the main public functions, especially distinguishing
    direct vs custom-kernel methodologies.

### Fixed
- **Phase sign**:
  - The sign of the lunar phase angle is now determined using `angular.get_phase_sign(...)`, based on the relative selenographic
    longitudes of the observer and the Sun.
  - This replaces the previous heuristic that compared the phase angle at a slightly later time to infer the sign,
    which could yield a wrong result around full and new Moon.

### Removed
- **Monolithic internals**:
  - The internal function `get_moon_datas_xyzs_no_zenith_azimuth(...)` has been removed.
  - All direct-geometry computations now go through `get_moon_datas_xyzs(...)` or `get_moon_datas_llhs(...)`,
    which always returns a fully populated `MoonData` (including zenith and azimuth).
  - Several private helpers tied to the old monolithic `spicedmoon.spicedmoon` implementation have been
    dropped in favor of the new modular design (`geometry`, `angular`, `custombody.core`, etc.).
  - Public behaviour is preserved via the new submodules and the `spicedmoon.deprecated` compatibility layer.

## [1.0.13] - 2023-10-12

Initial version that serves as the baseline for tracking changes in the change log.


[unreleased]: https://github.com/GOA-UVa/spicedmoon/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/GOA-UVa/spicedmoon/compare/v1.0.13...v1.1.0
[1.0.13]: https://github.com/GOA-UVa/spicedmoon/releases/tag/v1.0.13
