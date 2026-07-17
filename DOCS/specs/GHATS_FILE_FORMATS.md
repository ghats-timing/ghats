# GHATS File Format Specification

Draft: 2026-07-17

This note describes the file formats written and read by the current GHATS routines for native timing products and bispectrum products. It is a description of the implementation, not a new definition of the science conventions.

The source of truth is the IDL code. In particular:

- `FFT/mu_write_single_power.pro`
- `FFT/mu_write_single_fft.pro`
- `FFT/mu_compute_fft.pro`
- `ANALYSIS/ghats_getheader.pro`
- `ANALYSIS/ghats_openpds.pro`
- `ANALYSIS/ghats_openfft.pro`
- `ANALYSIS/read_pds_line.pro`
- `ANALYSIS/read_fft_line.pro`
- `ANALYSIS/gh_dyn_fits.pro`
- `BISPECTRUM/gh_bispec_fits.pro`
- `BISPECTRUM/gh_write_bispec_fits.pro`
- `BISPECTRUM/gh_write_cross_bispec_fits.pro`
- `BISPECTRUM/gh_read_bispec_fits.pro`

## Scope

This document covers:

- GHATS binary `.pds` files.
- GHATS binary `.fft` files.
- GHATS dynamic-PDS FITS files written by `gh_dyn_fits`.
- GHATS bispectrum FITS files.

This document does not cover:

- Mission event files read as input.
- Log files, PostScript plots, or screen output.
- Text parameter files, file lists, or metafiles.
- XSPEC/OGIP products such as `.pha` and `.rmf`.

## IDL Binary Products

GHATS `.pds` and `.fft` files are IDL unformatted binary files written with `WRITEU` and read with `READU`. The files use the native binary layout of the platform that wrote them. Current GHATS products are normally little-endian, but external readers should verify the header fields rather than assuming a product is valid.

Both formats have the same fixed header and the same per-transform metadata record. The final payload in each transform differs:

- `.pds` stores one real power array per transform.
- `.fft` stores one complex Fourier-amplitude array per transform.

### Shared Header

The header is written once at the beginning of the file by `mu_write_single_power` or `mu_write_single_fft`. It is read by `ghats_getheader`.

The write order is:

| Field | IDL type | Meaning |
| --- | --- | --- |
| `version_id` | string, 16 bytes | GHATS/version identifier. |
| `osserv` | string, 16 bytes | Observatory name. |
| `strumento` | string, 16 bytes | Instrument name. |
| `sorgente` | string, 16 bytes | Source or target name. |
| `rmjd` | double | Start time, in the GHATS MJD/RMJD convention used by the writer. |
| `nft` | long | Number of input time bins in each transform. |
| `T` | double | Frequency resolution, equal to `1.0 / transform_duration`. |
| `ntotal_ffts` | long | Number of transforms written to the file. |
| `canali` | int[2] | Selected channel range. |
| `proliferation` | int | Proliferation factor used by the writer. |
| `baryflag` | int | Barycentric-correction flag. |
| `n_spectral_bins` | int | Number of stored spectral channels; current writers set this to 3. |
| `background_flag` | int | Background flag; current writers set this to 0. |
| `dummy` | byte[100] | Reserved flags and mission-dependent metadata. |

The number of transforms, `ntotal_ffts`, begins at byte offset 84 in the current binary layout. Several producers initially write a placeholder value and patch this field after the final number of transforms is known.

The first bytes of `dummy` have established meanings in current products:

| Field | Meaning |
| --- | --- |
| `dummy[0]` | Bit-coded turbo/GTI information. Current writers use `1` for turbo and add `2` when the GTI flag is set. |
| `dummy[1]` | Window-function code. |

Some readers also inspect later `dummy` bytes for mission-specific quantities. Those bytes should be treated as reserved unless a product-specific writer documents them.

### Frequency Grid

For both `.pds` and `.fft`, the stored frequency bins are reconstructed as:

```idl
frequency = (FINDGEN(nft / 2) + 1.0) * T
```

Thus the DC bin is not stored. The stored array contains bins 1 through `nft/2`, including the Nyquist bin.

The segment duration represented by one transform is:

```idl
duration = 1.0 / T
```

The time resolution of the original binned light curve is:

```idl
dt = duration / nft
```

The Nyquist frequency is:

```idl
nyquist = nft * T / 2.0
```

### Per-Transform Record

After the header, each transform is written as one binary record with common scalar metadata followed by the transform payload.

The common write order is:

| Field | IDL type | Meaning |
| --- | --- | --- |
| `rmjd` | double | Time assigned to this transform. |
| `cnts` | float | Total counts in the transform after the windowing step used by the writer. |
| `fndet` | float | Number of active detectors, stored as a float. |
| `a0` | float | DC-related normalization field. For PDS products this is normally `2 * ABS(DC)`. For FFT products written by `mu_write_single_fft`, this is currently written as 0. |
| `poisson` | float | Poisson level stored with the transform. |
| `current_vle_rate` | float | VLE rate or mission-equivalent housekeeping value. |
| `std21` | float | Standard-2 or mission-dependent auxiliary value. |
| `std22` | float | Standard-2 or mission-dependent auxiliary value. |
| `std23` | float | Standard-2 or mission-dependent auxiliary value. |
| payload | format-dependent | Power or complex FFT data. |

The scalar metadata are read by `read_pds_line` and `read_fft_line`.

## `.pds` Files

A `.pds` file stores a power-density spectrum for each selected transform.

The payload in each per-transform record is:

| Field | IDL type | Length | Meaning |
| --- | --- | --- | --- |
| `pwr` | float array | `nft / 2` | Power values for frequency bins 1 through `nft/2`. |

### PDS Normalization

For normal GHATS PDS production, `mu_compute_fft` computes the Fourier transform and then writes powers as:

```idl
pwr = ABS(rdata[1:np/2])^2 * 2.0 / cnts
```

where `cnts` is the total number of counts in the windowed input segment. This is Leahy normalization.

For RXTE/PCA products, the Poisson level may be computed with the Zhang et al. prescription through `mu_zhang`. For non-PCA products, the stored Poisson level is normally 2.0 unless the mission writer sets something else.

## `.fft` Files

A `.fft` file stores the complex Fourier amplitudes for each selected transform.

The payload in each per-transform record is:

| Field | IDL type | Length | Meaning |
| --- | --- | --- | --- |
| `rdata` | complex array | `nft / 2` | Complex Fourier amplitudes for frequency bins 1 through `nft/2`. |

IDL complex values are stored as two 32-bit floating-point numbers per element: real part followed by imaginary part.

### FFT Convention

The stored FFT array is produced in `mu_compute_fft` as:

```idl
rdata = FFT(rdata, 1, /OVERWRITE)
fftdata = rdata[1:np/2]
```

The DC bin is not stored. GHATS therefore stores the positive-frequency bins produced by IDL `FFT(input, 1)`.

The source comment in `mu_compute_fft` notes that the IDL call used here corresponds to the forward transform in the convention used by the original GHATS timing routines. External code that computes cross spectra or lags should match this implementation convention rather than changing the lag sign convention.

### FFT Units

The `.fft` payload is the raw complex Fourier amplitude. It is not saved in Leahy, rms, or fractional-rms normalization. Normalizations are applied later by analysis routines that read the FFT products.

## Dynamic-PDS FITS Files

`gh_dyn_fits` writes a dynamic PDS or time-frequency image using `WRITEFITS`.

The current HDU order is:

| HDU | Contents | Notes |
| --- | --- | --- |
| Primary | `dynima` | Dynamic PDS image. |
| Extension 1 | `nu` | Frequency array. |
| Extension 2 | `time` | Time array. |
| Extension 3 | `rate` | Rate array. |

The routine does not currently add custom `EXTNAME` keywords or detailed axis metadata. Readers of these files should therefore treat the HDU order above as part of the current format.

## Bispectrum FITS Files

The current GHATS bispectrum FITS format is image-extension based and is written by `gh_bispec_fits` through `gh_write_bispec_fits`.

The primary HDU contains a dummy one-byte image. This is intentional: it avoids problems with FITS writers that do not handle a zero-length primary image reliably.

### Primary Header

The primary header records product-level metadata. Common keywords include:

| Keyword | Meaning |
| --- | --- |
| `CREATOR` | Writer routine, normally `GH_BISPEC_FITS`. |
| `CONTENT` | Product description, normally `2D bispectrum products`. |
| `PRODUCT` | Product label, such as `BISPEC` or a cross-bispectrum label. |
| `SRCFILE1`, `SRCFILE2`, `SRCFILE3` | Optional source files used to make the product. |
| `NFREQ1`, `NFREQ2` | Number of frequency bins on the two bispectrum axes. |
| `F1MIN`, `F1MAX` | First-axis frequency range. |
| `F2MIN`, `F2MAX` | Second-axis frequency range. |
| `DF1`, `DF2` | Frequency spacing on each axis. |
| `BUNITB` | Unit string for `BREAL`, `BIMAG`, and `BMOD`. |
| `RAWSUMS` | Flag indicating whether raw-sum extensions are present. |
| `INTBICOH` | Flag indicating whether intrinsic-observed bicoherence extensions are present. |
| `NEXTEND` | Number of extensions expected by the writer. |
| `REBIN1`, `REBIN2` | Rebin factors for the two frequency axes. |
| `AXLOG1`, `AXLOG2` | Flags indicating logarithmic axis construction. |

For intrinsic-observed bicoherence products, the primary header may also include:

| Keyword | Meaning |
| --- | --- |
| `POIMETH` | Poisson correction method. |
| `POILEVEL` | Poisson level used. |
| `POIFMIN`, `POIFMAX` | Frequency range used for the Poisson estimate, when applicable. |
| `BICNUM` | Numerator convention. Current intrinsic-observed products use `UNCHANGED`. |
| `BICDEN` | Denominator convention. Current intrinsic-observed products use `POISSON_CORR`. |

Cross-bispectrum products may also include:

| Keyword | Meaning |
| --- | --- |
| `CORRMODE` | Correction mode. |
| `DIAGCORR` | Diagonal correction convention. |
| `DTCOVAR` | Dead-time covariance convention. |
| `NUMBIAS` | Numerator bias convention. |
| `BAND1`, `BAND2`, `BAND3` | Energy-band labels. |
| `BANDOVLP` | Band-overlap description. |
| `XORDER` | Cross-bispectrum ordering convention. |
| `FCONV` | Frequency or lag-sign convention note. |

History records from the writer are stored as FITS `HISTORY` cards.

### Required Extensions

The standard bispectrum product contains these extensions, in this order:

| EXTNAME | Type | Shape | Units |
| --- | --- | --- | --- |
| `BREAL` | image | `[NFREQ1, NFREQ2]` | `BUNITB` |
| `BIMAG` | image | `[NFREQ1, NFREQ2]` | `BUNITB` |
| `BMOD` | image | `[NFREQ1, NFREQ2]` | `BUNITB` |
| `BPHASE` | image | `[NFREQ1, NFREQ2]` | radians |
| `BICOH` | image | `[NFREQ1, NFREQ2]` | dimensionless |
| `NPROD_USED` | long image | `[NFREQ1, NFREQ2]` | count |
| `FREQ1` | double image | `[NFREQ1]` | Hz |
| `FREQ2` | double image | `[NFREQ2]` | Hz |

The two-dimensional image extensions include frequency-axis metadata:

| Keyword | Axis 1 meaning | Axis 2 meaning |
| --- | --- | --- |
| `CTYPE1`, `CTYPE2` | `FREQ1` | `FREQ2` |
| `CUNIT1`, `CUNIT2` | `Hz` | `Hz` |
| `CRPIX1`, `CRPIX2` | Reference pixel, normally 1 | Reference pixel, normally 1 |
| `CRVAL1`, `CRVAL2` | First frequency value | First frequency value |
| `CDELT1`, `CDELT2` | Frequency spacing | Frequency spacing |

The image extensions also carry the rebinning keywords `REBIN1`, `REBIN2`, `AXLOG1`, and `AXLOG2`.

The unit stored in `BUNITB` depends on the requested bispectrum normalization:

| Normalization | `BUNITB` |
| --- | --- |
| rms | `rms^3` |
| Leahy | `Leahy^(3/2)` |
| other/native | `native` |

### Optional Raw-Sum Extensions

If raw sums are present, the writer appends:

| EXTNAME | Type | Shape | Meaning |
| --- | --- | --- | --- |
| `RAW_BSUM_REAL` | image | `[NFREQ1, NFREQ2]` | Real part of the raw bispectrum sum. |
| `RAW_BSUM_IMAG` | image | `[NFREQ1, NFREQ2]` | Imaginary part of the raw bispectrum sum. |
| `RAW_DEN1` | image | `[NFREQ1, NFREQ2]` | First raw denominator sum. |
| `RAW_DEN2` | image | `[NFREQ1, NFREQ2]` | Second raw denominator sum. |

These extensions use the same frequency-axis metadata as the required two-dimensional image extensions.

### Optional Intrinsic-Observed Extensions

If intrinsic-observed bicoherence products are present, the writer appends:

| EXTNAME | Type | Shape | Meaning |
| --- | --- | --- | --- |
| `INT_BICOH` | image | `[NFREQ1, NFREQ2]` | Intrinsic-observed bicoherence estimate. |
| `RAW_DEN1_CORR` | image | `[NFREQ1, NFREQ2]` | Corrected first denominator sum. |
| `RAW_DEN2_CORR` | image | `[NFREQ1, NFREQ2]` | Corrected second denominator sum. |

These extensions carry the Poisson-correction metadata used for the intrinsic-observed product.

### Cross-Bispectrum Extensions

`gh_write_cross_bispec_fits` uses the same core format and may append additional extensions for cross-bispectrum diagnostics:

| EXTNAME | Type | Shape | Meaning |
| --- | --- | --- | --- |
| `INT_BICOH_INDEP` | image | `[NFREQ1, NFREQ2]` | Independent-noise intrinsic bicoherence estimate. |
| `RAW_DEN1_CORR_INDEP` | image | `[NFREQ1, NFREQ2]` | Independent-noise corrected first denominator. |
| `RAW_DEN2_CORR_INDEP` | image | `[NFREQ1, NFREQ2]` | Independent-noise corrected second denominator. |
| `INTR_VALID` | byte image | `[NFREQ1, NFREQ2]` | Validity mask for intrinsic quantities. |
| `CORR_FLAGS` | long image | `[NFREQ1, NFREQ2]` | Correction flags. |
| `DIAG_FLAG` | byte image | `[NFREQ1, NFREQ2]` | Diagonal-correction flag. |

These extensions include the same frequency-axis metadata as the main two-dimensional bispectrum images, plus correction-convention keywords when supplied by the caller.

### Reading Bispectrum FITS Products

`gh_read_bispec_fits` reads current products primarily by `EXTNAME`. It also contains fallback behavior for older files with less complete extension naming. New readers should prefer `EXTNAME` over extension number wherever possible.

## Compatibility Notes

The binary `.pds` and `.fft` formats are compact and historically tied to IDL's native binary representation. The FITS products are more self-describing, especially the bispectrum products, where extensions are named and frequency-axis metadata are written.

When checking compatibility between `.pds` files, routines should compare at least:

- `nft`
- frequency resolution `T`
- transform duration `1.0 / T`
- Nyquist frequency `nft * T / 2.0`
- channel range
- observatory and instrument metadata, when relevant to the analysis
- normalization and Poisson-handling assumptions implied by the producing routine

For `.fft` files, compatibility checks should additionally ensure that later cross-spectral or lag calculations use the same GHATS Fourier and lag conventions.
