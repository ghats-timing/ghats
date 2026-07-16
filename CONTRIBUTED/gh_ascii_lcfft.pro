;+
; NAME:
;      GH_ASCII_LCFFT
;
; PURPOSE:
;      Produce GHATS-compatible .pds or .fft files from an ASCII light curve.
;
; EXPLANATION:
;      GH_ASCII_LCFFT reads an ASCII light curve with columns
;
;          time(s)   counts_or_rate   error
;
;      and writes a GHATS-style power-spectrum file (.pds) or complex
;      Fourier-spectrum file (.fft).  The output is written through the
;      standard GHATS routines MU_COMPUTE_FFT, MU_WRITE_SINGLE_POWER, and
;      MU_WRITE_SINGLE_FFT, so that the binary file format is the same as
;      the one produced by GH_XTE.
;
;      The third column is read only to tolerate the input format. It is
;      not used in the FFT/PDS calculation.
;
;      By default, the second column is assumed to contain counts per
;      original light-curve bin. If the second column is instead a count
;      rate, use /RATE.
;
;      The default time resolution is 0.016 s. This is appropriate for the
;      barycentered HXMT/ULTRACAM light curves discussed here. For other
;      data sets, use DT=...
;
; CALLING SEQUENCE:
;      GH_ASCII_LCFFT, infile, outtype, treb, npds, outfile
;
;      GH_ASCII_LCFFT, infile, outtype, treb, npds, outfile, $
;                      DT=dt, MJD0=mjd0, /RATE
;
; INPUTS:
;      INFILE    : ASCII light-curve filename.
;
;      OUTTYPE   : Output type. Accepted values are:
;                    'POWER', 'PDS'  -> write .pds file
;                    'FFT', 'FOURIER' -> write .fft file
;
;      TREB      : Integer rebinning factor in time.
;                  Counts are summed, not averaged.
;
;      NPDS      : Number of rebinned light-curve points per FFT.
;                  Must be a power of two.
;
;      OUTFILE   : Output filename.
;
; KEYWORDS:
;      DT        : Original light-curve time resolution, in seconds.
;                  Default: 0.016d0.
;
;      MJD0      : MJD corresponding to time=0 in the ASCII file.
;                  Default: 60196.0d0.
;
;      RATE      : If set, the second column is interpreted as count rate.
;                  The routine converts it to counts per bin using DT before
;                  rebinning.
;
;      GAP_TOL   : Gap tolerance in units of DT.
;                  A gap is detected when
;
;                      ABS((time[i+1]-time[i]) - DT) GT GAP_TOL*DT
;
;                  Default: 0.5d0.
;
;      WIND      : Window type passed to GHATS_WINDOW.
;                  Default: 'Boxcar'.
;
;      WPAR      : Optional parameter for window functions that require one.
;
;      SOURCE    : Source name written in the GHATS header.
;                  Default: 'ASCII'.
;
;      OBSERVATORY :
;                  Observatory name written in the GHATS header.
;                  Default: 'ASCII'.
;
;      INSTRUMENT :
;                  Instrument name written in the GHATS header.
;                  Default: 'LC'.
;
;      HELP      : Print concise help and return.
;
; OUTPUTS:
;      A GHATS-compatible .pds or .fft file.
;
; DESIGN NOTES:
;      1. Frequency range:
;         MU_COMPUTE_FFT writes Fourier bins 1:NPDS/2. The DC bin is not
;         written. The Nyquist bin is included, as in GH_XTE. The frequency
;         resolution is 1/T, where T = NPDS * DT * TREB.
;
;      2. NaN masking:
;         Rows with non-finite time or non-finite counts/rates are treated
;         as boundaries. They are not interpolated over. This avoids mixing
;         bad data into FFT segments and keeps behaviour close to GH_XTE,
;         where invalid intervals are skipped.
;
;      3. Gap handling:
;         The routine never rebins or FFTs across gaps. If a continuous
;         stretch does not contain enough points for an integer number of
;         rebinned bins or FFT segments, the leftover points are dropped.
;
;      4. Rebinning:
;         This is strictly one-dimensional time rebinning. Counts are summed.
;         Future two-dimensional rebinning, for example simultaneous
;         time/energy rebinning before Fourier products, should not be added
;         here without checking consistency with the cross-spectral routines.
;
;      5. Consistency with GH_CS and GH_CROSS_RE_IM_NEW:
;         The routine preserves the GHATS convention of equal-length FFT
;         segments, positive Fourier bins only, raw complex FFT output for
;         .fft files, and counts-based input to MU_COMPUTE_FFT. For cross
;         spectra, the different input light curves must have the same time
;         grid and the same accepted/gapped intervals, otherwise the segment
;         numbers will not correspond one-to-one.
;
; EXAMPLES:
;      ; Make a PDS, assuming counts per 16-ms bin:
;      GH_ASCII_LCFFT, 'lc_HE.txt', 'POWER', 1, 4096, 'lc_HE.pds'
;
;      ; Make an FFT file:
;      GH_ASCII_LCFFT, 'lc_HE.txt', 'FFT', 1, 4096, 'lc_HE.fft'
;
;      ; Data with 1 ms bins:
;      GH_ASCII_LCFFT, 'lc.txt', 'FFT', 1, 8192, 'lc.fft', DT=0.001d0
;
;      ; If the second column is count rate rather than counts/bin:
;      GH_ASCII_LCFFT, 'lc_rate.txt', 'POWER', 4, 4096, 'lc.pds', $
;                      DT=0.016d0, /RATE
;
;      ; Print help:
;      GH_ASCII_LCFFT, /HELP
;
; COMMON BLOCKS:
;      sis
;      versione
;      barycentered
;      finestre
;
; ROUTINES USED:
;      MU_COMPUTE_FFT
;      MU_WRITE_SINGLE_POWER
;      MU_WRITE_SINGLE_FFT
;      GHATS_WINDOW
;      MU_POWER_OF_TWO
;
; MODIFICATION HISTORY:
;      M. Mendez  14 May 2026
;          First version for ASCII barycentered light curves.
;
;      Developed with assistance from ChatGPT, based on GH_XTE,
;      MU_COMPUTE_FFT, MU_WRITE_SINGLE_POWER, MU_WRITE_SINGLE_FFT,
;      GHATS_WINDOW, and GHATS_GETHEADER.
;-

PRO gh_ascii_lcfft, infile, outtype, treb, npds, outfile, $
                    dt=dt, mjd0=mjd0, rate=rate, $
                    gap_tol=gap_tol, wind=wind, wpar=wpar, $
                    source=source, observatory=observatory, $
                    instrument=instrument, help=help

  ;------------------------------------------------------------
  ; Common blocks used by the standard GHATS writing routines.
  ;------------------------------------------------------------
  COMMON sis, sistema
  COMMON versione, version_id
  COMMON barycentered, baryflag
  COMMON finestre, finestra, winn

  ;------------------------------------------------------------
  ; Help.
  ;------------------------------------------------------------
  IF KEYWORD_SET(help) THEN BEGIN
     PRINT, ''
     PRINT, 'GH_ASCII_LCFFT'
     PRINT, ''
     PRINT, 'Usage:'
     PRINT, '  GH_ASCII_LCFFT, infile, outtype, treb, npds, outfile'
     PRINT, ''
     PRINT, 'Examples:'
     PRINT, "  GH_ASCII_LCFFT, 'lc_HE.txt', 'POWER', 1, 4096, 'lc_HE.pds'"
     PRINT, "  GH_ASCII_LCFFT, 'lc_HE.txt', 'FFT',   1, 4096, 'lc_HE.fft'"
     PRINT, "  GH_ASCII_LCFFT, 'lc.txt', 'FFT', 1, 8192, 'lc.fft', DT=0.001d0"
     PRINT, "  GH_ASCII_LCFFT, 'lc_rate.txt', 'POWER', 4, 4096, 'lc.pds', /RATE"
     PRINT, ''
     PRINT, 'Input ASCII columns: time(s), counts_or_rate, error'
     PRINT, 'The third column is ignored.'
     PRINT, ''
     PRINT, 'Defaults: DT=0.016 s, MJD0=60196.0, GAP_TOL=0.5'
     PRINT, ''
     RETURN
  ENDIF

  ;------------------------------------------------------------
  ; Basic argument check.
  ;------------------------------------------------------------
  IF N_PARAMS() NE 5 THEN BEGIN
     PRINT, 'Usage: GH_ASCII_LCFFT, infile, outtype, treb, npds, outfile'
     PRINT, '       Use /HELP for examples.'
     RETURN
  ENDIF

  ;------------------------------------------------------------
  ; Defaults.
  ;------------------------------------------------------------
  IF N_ELEMENTS(dt) EQ 0 THEN dt = 0.016d0
  IF N_ELEMENTS(mjd0) EQ 0 THEN mjd0 = 60196.0d0
  IF N_ELEMENTS(gap_tol) EQ 0 THEN gap_tol = 0.5d0
  IF N_ELEMENTS(wind) EQ 0 THEN wind = 'Boxcar'

  IF N_ELEMENTS(source) EQ 0 THEN source = 'ASCII'
  IF N_ELEMENTS(observatory) EQ 0 THEN observatory = 'ASCII'
  IF N_ELEMENTS(instrument) EQ 0 THEN instrument = 'LC'

  ; The data are barycentered. This flag is written in the GHATS header.
  baryflag = 1

  ; If the external GHATS version string is not already defined, create one.
  ; The writers expect a 16-character GHATS-style string.
  IF N_ELEMENTS(version_id) EQ 0 THEN version_id = 'GHATSLC0001     '
  IF STRMID(version_id,0,5) NE 'GHATS' THEN version_id = 'GHATSLC0001     '

  ; Some GHATS routines keep this common block. It is not critical here,
  ; but setting it avoids undefined common-block surprises.
  IF N_ELEMENTS(sistema) EQ 0 THEN sistema = 'IDL'

  ;------------------------------------------------------------
  ; Normalise output type.
  ;------------------------------------------------------------
  outlow = STRLOWCASE(STRMID(STRTRIM(outtype,2),0,1))

  CASE outlow OF
     'p': chout = 'POWER'
     'f': chout = 'FOURIER'
     ELSE: BEGIN
        PRINT, 'Unknown output type: ', outtype
        PRINT, 'Use POWER/PDS or FFT/FOURIER.'
        RETURN
     END
  ENDCASE

  ;------------------------------------------------------------
  ; Check rebinning and FFT length.
  ;------------------------------------------------------------
  treb_l = LONG(treb)
  IF treb_l LT 1L THEN BEGIN
     PRINT, 'TREB must be >= 1.'
     RETURN
  ENDIF

  np = LONG(npds)
  IF np LT 1L THEN BEGIN
     PRINT, 'NPDS must be >= 1.'
     RETURN
  ENDIF

  IF MU_POWER_OF_TWO(np) NE 1 THEN BEGIN
     PRINT, 'NPDS must be a power of two.'
     RETURN
  ENDIF

  ; Effective time resolution after rebinning.
  tres_fft = dt * DOUBLE(treb_l)
  T = DOUBLE(np) * tres_fft

  ;------------------------------------------------------------
  ; Read the ASCII file.
  ; We use two passes so that the arrays can be allocated once.
  ;------------------------------------------------------------
  line = ''
  nrows = 0L

  OPENR, lun, infile, /GET_LUN, ERROR=err
  IF err NE 0 THEN BEGIN
     PRINT, 'Cannot open input file: ', infile
     RETURN
  ENDIF

  WHILE ~EOF(lun) DO BEGIN
     READF, lun, line
     line2 = STRTRIM(line,2)
     IF line2 NE '' THEN BEGIN
        tok = STRSPLIT(line2, /EXTRACT)
        IF N_ELEMENTS(tok) GE 2 THEN nrows = nrows + 1L
     ENDIF
  ENDWHILE

  FREE_LUN, lun

  IF nrows LT 2L THEN BEGIN
     PRINT, 'Not enough rows in input light curve.'
     RETURN
  ENDIF

  time = DBLARR(nrows)
  value = DBLARR(nrows)
  errcol = DBLARR(nrows)

  OPENR, lun, infile, /GET_LUN
  irow = 0L

  WHILE ~EOF(lun) DO BEGIN
     READF, lun, line
     line2 = STRTRIM(line,2)

     IF line2 NE '' THEN BEGIN
        tok = STRSPLIT(line2, /EXTRACT)

        IF N_ELEMENTS(tok) GE 2 THEN BEGIN
           time[irow] = DOUBLE(tok[0])
           value[irow] = DOUBLE(tok[1])

           IF N_ELEMENTS(tok) GE 3 THEN BEGIN
              errcol[irow] = DOUBLE(tok[2])
           ENDIF ELSE BEGIN
              errcol[irow] = !VALUES.D_NAN
           ENDELSE

           irow = irow + 1L
        ENDIF
     ENDIF
  ENDWHILE

  FREE_LUN, lun

  ; In case the second pass accepted fewer rows than the first pass.
  IF irow LT nrows THEN BEGIN
     time = time[0:irow-1L]
     value = value[0:irow-1L]
     errcol = errcol[0:irow-1L]
     nrows = irow
  ENDIF

  ;------------------------------------------------------------
  ; Convert rates to counts per original bin if requested.
  ; MU_COMPUTE_FFT assumes counts per bin, because it uses
  ; CNTS = TOTAL(RDATA).
  ;------------------------------------------------------------
  counts0 = DBLARR(nrows)

  IF KEYWORD_SET(rate) THEN BEGIN
     counts0 = value * dt
  ENDIF ELSE BEGIN
     counts0 = value
  ENDELSE

  ;------------------------------------------------------------
  ; Identify finite data.
  ; The error column is deliberately not used. Non-finite time or counts
  ; define a break in the light curve.
  ;------------------------------------------------------------
  valid = BYTE(FINITE(time) * FINITE(counts0))

  ;------------------------------------------------------------
  ; Count how many complete FFTs will be written.
  ; This must be known before the first call to MU_COMPUTE_FFT, because
  ; it is written in the GHATS file header.
  ;------------------------------------------------------------
  ntotal_ffts = 0L
  istart = 0L

  WHILE istart LT nrows DO BEGIN

     ; Skip invalid rows.
     WHILE (istart LT nrows) AND (valid[istart] EQ 0B) DO istart = istart + 1L
     IF istart GE nrows THEN BREAK

     iend = istart

     ; Extend current continuous stretch.
     WHILE iend LT nrows-1L DO BEGIN
        IF valid[iend+1L] EQ 0B THEN BREAK
        dti = time[iend+1L] - time[iend]
        IF ABS(dti - dt) GT gap_tol*dt THEN BREAK
        iend = iend + 1L
     ENDWHILE

     nraw = iend - istart + 1L
     nreb = nraw / treb_l
     ntotal_ffts = ntotal_ffts + LONG(nreb / np)

     istart = iend + 1L
  ENDWHILE

  IF ntotal_ffts LT 1L THEN BEGIN
     PRINT, 'No complete FFT segments found.'
     PRINT, 'Check DT, GAP_TOL, TREB, and NPDS.'
     Print,DT, GAP_TOL, TREB, NPDS
     RETURN
  ENDIF

  ;------------------------------------------------------------
  ; Window function.
  ; The window index WINN is stored in the GHATS dummy header.
  ;------------------------------------------------------------
  IF ((wind EQ 'Hamming') OR (wind EQ 'Triplet') OR $
      (wind EQ 'Gauss')   OR (wind EQ 'Kaiser')) THEN BEGIN

     IF N_ELEMENTS(wpar) EQ 0 THEN BEGIN
        PRINT, 'The selected window function requires WPAR.'
        RETURN
     ENDIF

     finestra = GHATS_WINDOW(np, wind, PAR=wpar, WINN=winn)

  ENDIF ELSE BEGIN
     finestra = GHATS_WINDOW(np, wind, WINN=winn)
  ENDELSE

  IF N_ELEMENTS(finestra) NE np THEN BEGIN
     PRINT, 'Problem creating window function.'
     RETURN
  ENDIF

  IF TOTAL(finestra LT 0.0) GT 0 THEN BEGIN
     PRINT, 'Invalid window function.'
     RETURN
  ENDIF

  ;------------------------------------------------------------
  ; Metadata for the GHATS header.
  ; These are placeholders because the ASCII light curve has no FITS header.
  ; Keep them short because the writers use fixed 16-character strings.
  ;------------------------------------------------------------
  canali = INTARR(2)
  canali[0] = 0
  canali[1] = 0

  ndet = 1
  i_vle = 3
  current_vle_rate = 0.0
  std21 = 0.0
  std22 = 0.0
  std23 = 0.0

  BANDE = [0,0,0,0,0,0]
  turbo = 0
  gti_flag = 0
  proliferation = 1

  new_output_flag = 1
  output_unit = 9

  pwr = FLTARR(np/2)
  rdata = FLTARR(np)

  nffts = 0L

  PRINT, 'Input file                  : ', infile
  PRINT, 'Output file                 : ', outfile
  PRINT, 'Output type                 : ', chout
  PRINT, 'Original dt (s)             : ', STRTRIM(STRING(dt),2)
  PRINT, 'Rebin factor                : ', STRTRIM(STRING(treb_l),2)
  PRINT, 'Effective dt (s)            : ', STRTRIM(STRING(tres_fft),2)
  PRINT, 'Points per FFT              : ', STRTRIM(STRING(np),2)
  PRINT, 'FFT length (s)              : ', STRTRIM(STRING(T),2)
  PRINT, 'Total FFTs                  : ', STRTRIM(STRING(ntotal_ffts),2)
  PRINT, 'MJD at time zero            : ', STRTRIM(STRING(mjd0),2)
  PRINT, 'Gap tolerance               : ', STRTRIM(STRING(gap_tol),2), ' * dt'
  IF KEYWORD_SET(rate) THEN PRINT, 'Second column interpreted as count rate.'
  IF ~KEYWORD_SET(rate) THEN PRINT, 'Second column interpreted as counts/bin.'
  PRINT, 'Third column ignored.'
  PRINT, ''

  ;------------------------------------------------------------
  ; Main processing loop.
  ; We repeat the same stretch finding used above, now actually filling
  ; rebinned light curves and calling MU_COMPUTE_FFT.
  ;------------------------------------------------------------
  istart = 0L

  WHILE istart LT nrows DO BEGIN

     ; Skip invalid rows.
     WHILE (istart LT nrows) AND (valid[istart] EQ 0B) DO istart = istart + 1L
     IF istart GE nrows THEN BREAK

     iend = istart

     WHILE iend LT nrows-1L DO BEGIN
        IF valid[iend+1L] EQ 0B THEN BREAK
        dti = time[iend+1L] - time[iend]
        IF ABS(dti - dt) GT gap_tol*dt THEN BREAK
        iend = iend + 1L
     ENDWHILE

     nraw = iend - istart + 1L
     nreb = nraw / treb_l

     ; Drop incomplete rebinned bins at the end of the stretch.
     IF nreb GE np THEN BEGIN

        reb_counts = DBLARR(nreb)
        reb_time = DBLARR(nreb)

        FOR ir=0L, nreb-1L DO BEGIN
           j1 = istart + ir*treb_l
           j2 = j1 + treb_l - 1L

           ; Counts are summed. Do not average.
           reb_counts[ir] = TOTAL(counts0[j1:j2])
           reb_time[ir] = time[j1]
        ENDFOR

        nseg = nreb / np

        FOR iseg=0L, nseg-1L DO BEGIN
           k1 = iseg*np
           k2 = k1 + np - 1L

           rdata = FLOAT(reb_counts[k1:k2])

           ; The ASCII times are seconds from midnight MJD0.
           ; This gives the segment start time in MJD, matching the
           ; convention expected by the GHATS writers.
           tstart_fft = reb_time[k1]
           start_time = mjd0 + DOUBLE(tstart_fft)/86400.0d0

           MU_COMPUTE_FFT, rdata, pwr, np, T, start_time, $
                           source, observatory, instrument, $
                           outfile, chout, ndet, i_vle, $
                           current_vle_rate, tstart_fft, $
                           new_output_flag, output_unit, ntotal_ffts, $
                           canali, std21, std22, std23, $
                           BANDE, turbo, gti_flag, proliferation

           nffts = nffts + 1L
           PRINT, 'FFT# ', STRTRIM(STRING(nffts),2), '/', $
                 STRTRIM(STRING(ntotal_ffts),2), $
                 '  Start time(s): ', STRTRIM(STRING(tstart_fft),2)

        ENDFOR
     ENDIF

     istart = iend + 1L
  ENDWHILE

  ;------------------------------------------------------------
  ; Sanity check and header correction.
  ; This mirrors the protection in GH_XTE. It should normally not trigger.
  ; Header position 84 is where GH_XTE rewrites the number of FFTs.
  ;------------------------------------------------------------
  IF nffts NE ntotal_ffts THEN BEGIN
     PRINT, 'Readjusting number of FFTs to ', nffts
     OPENU, uu, outfile, /GET_LUN
     POINT_LUN, uu, 84
     WRITEU, uu, nffts
     FREE_LUN, uu
  ENDIF

  CLOSE, /ALL

  PRINT, 'gh_ascii_lcfft: normal termination'

END
