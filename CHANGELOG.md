# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

[//]: # "## [unreleased] - yyyy-mm-dd"

## [unreleased] - yyyy-mm-dd

### Added
- `dt_to_str` throws a runtime warning if any `datetime` is timezone-naive.
- Direct calculation can now obtain zenith and azimuth values, without using custom intermediate kernel.

### Fixed
- Phase sign is now determined using the relative selenographic longitudes of the observer and the Sun,
  replacing the previous method that inferred the sign by comparing the phase angle at a later time.

### Changed
- Changed project license from MIT to LGPL.
- Modernized project structure: consolidated configuration into `pyproject.toml`, added a `CHANGELOG`, and applied current best practices.

## [1.0.13] - 2023-10-12

Initial version that serves as the baseline for tracking changes in the change log.


[unreleased]: https://github.com/GOA-UVa/spicedmoon/compare/v1.0.13...HEAD
[1.0.13]: https://github.com/GOA-UVa/spicedmoon/releases/tag/v1.0.13
