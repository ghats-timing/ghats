# Changelog

All notable changes to GHATS are documented in this file.

## [Unreleased]

## [4.0.0] - 2026-08-29

### Added

#### Bispectral analysis

- Added a new `BISPECTRUM` module for higher-order Fourier timing analysis.
- Added one- and two-dimensional bispectrum calculations.
- Added bicoherence and normalized bispectral products.
- Added hypotenuse-domain bispectral analysis.
- Added cross-bispectrum calculations and segment-matching tools.
- Added frequency-band, Lorentzian-pair, and target-pair selection utilities.
- Added two-dimensional bispectrum rebinning and plotting.
- Added FITS input and output for bispectrum and cross-bispectrum products.

Principal new routines include:

- `gh_bispec`
- `gh_bispec2d`
- `gh_bispec_hypot`
- `gh_cross_bispec2d`
- `gh_cross_bispec2d_match`
- `gh_cross_bispec_target_segments`
- `gh_bispec2d_rebin`
- `gh_bispec_plot`
- `gh_read_bispec_fits`
- `gh_write_bispec_fits`
- `gh_write_cross_bispec_fits`

#### Mission support

- Added support for Einstein Probe.
- Added support for IXPE.
- Added mission-specific event readers, channel-selection routines, time-resolution detection, and FITS metadata handling for both missions.

#### Analysis utilities

- Added `ghx_read_pds_metafile` for reading and validating power-spectrum metafiles.
- Added FFT generation from ASCII light curves.
- Added utilities for splitting FFT files and power-spectrum files.
- Added a utility for concatenating power-spectrum products.
- Added Poisson-level estimation over selectable frequency intervals.
- Extended `gh_info` to inspect both GHATS PDS and FFT products.
- Added quick-look power-spectrum generation directly from GHATS FFT files in `ghats_all`.

#### Help system

- Added a central GHATS help routine.
- Added searchable help keywords to user-facing routines.
- Added alphabetical and task-oriented routine indices.
- Added a command-line help utility.
- Added `/HELP` support and expanded usage information for AstroSat CZTI.

#### Documentation

- Added a LaTeX-based GHATS User Manual.
- Added chapters covering installation, quick-start workflows, FFT and power-spectrum production, averaging and rebinning, dynamic power spectra, cross-spectral analysis, plotting, and XSPEC use.
- Added file-format, conversion, and revision-history appendices.
- Added a standalone GHATS File Format Specification.
- Added manual source files, styles, front matter, and figures.

### Changed

#### FFT and power-spectrum handling

- Standardized FFT-length values as IDL `LONG` integers throughout mission front ends and low-level product-writing routines.
- Enforced storage of the FFT-length `NFT` header field as an IDL `LONG`, as required by the GHATS binary-header layout.
- Improved validation and reporting of GHATS FFT and PDS headers.
- Added diagnostic handling for older products whose `NFT` header field was incorrectly written as a 16-bit integer.
- Improved FFT quick-look processing and product-specific labelling.
- Improved metadata reporting for FFT and power-spectrum products.
- Updated common FFT-writing routines.
- Updated RXTE and XMM-Newton FFT processing.

#### Core analysis

- Substantially updated `ghx`, `gh_info`, and `ghats_all`.
- Improved handling of power-spectrum metafiles.
- Updated averaging, rebinning, plotting, PHA production, response handling, and XSPEC export.
- Updated dynamic power-spectrum routines.
- Updated light-curve and housekeeping routines.
- Updated version reporting and internal routine information.
- Improved Poisson-noise handling.
- Added safeguards for quick-look plots containing no positive powers.

#### Cross-spectral analysis

- Updated cross-spectrum generation.
- Updated real and imaginary cross-spectrum calculations.
- Improved metadata and output handling.
- Improved integration with the common GHATS workflow.

#### Mission interfaces

- Updated support for RXTE.
- Updated AstroSat CZTI and LAXPC interfaces.
- Updated Insight-HXMT HE, ME, and LE interfaces.
- Updated NICER, NuSTAR, Swift, and XMM-Newton interfaces.
- Improved consistency of FFT generation and metadata handling across missions.

#### Build and repository organization

- Updated compilation lists and the default build configuration.
- Added repository-wide ignore rules for generated and local files.
- Removed internal smoke-test routines from the distributed help indices.

### Fixed

- Fixed malformed GHATS FFT and PDS headers produced when an FFT length supplied as an IDL `INT` was written in place of the required `LONG`, shifting subsequent header fields on disk.
- Fixed the corresponding FFT-length type handling across RXTE, AstroSat, Einstein Probe, Insight-HXMT, IXPE, NICER, NuSTAR, Swift, and XMM-Newton front ends.
- Fixed low-level FFT and power-spectrum writers so that the `NFT` header field is always stored as an IDL `LONG`.
- Fixed misleading or invalid header-derived quantities when inspecting legacy or malformed products.
- Fixed FFT quick-look processing so that powers are calculated from individual FFT records rather than from averaged complex amplitudes.
- Fixed failures in logarithmic quick-look plots when no positive powers are available.
- Fixed an incomplete obsolete CZTI interactive-default branch that prevented the distributed default configuration from compiling.
- Fixed internal diagnostic routines being included in the public help indices.

### Compatibility and upgrade notes

- The intended GHATS binary-header layout has not changed: the `NFT` field is stored as an IDL `LONG`.
- Some products generated through older command-line paths may have written this field as a 16-bit IDL `INT`, shifting the remaining header fields and producing malformed files.
- `gh_info` and `ghats_all` can recognize this specific malformed layout and issue a warning, but affected FFT or PDS products should be regenerated before scientific analysis.
- Correctly written products from earlier GHATS versions are not expected to be affected by this correction.
- Users with independent software that reads GHATS binary products should verify that the `NFT` field is interpreted as a 32-bit integer and consult the updated GHATS File Format Specification.

