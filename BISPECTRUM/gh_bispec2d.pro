pro gh_bispec2d,filename,frequency1,frequency2, $
                breal,bimag,bmod,bphase,bicoh,nprod_used,nseg, $
                f1range=f1range,f2range=f2range, $
                index=index,time=time,sel=sel,rate=rate, $
                rms=back,leahy=leahy, $
	                raw_bsum_real=raw_bsum_real, $
	                raw_bsum_imag=raw_bsum_imag, $
	                raw_den1=raw_den1, $
	                raw_den2=raw_den2, $
	                int_bicoh=int_bicoh, $
	                den1_corr=den1_corr_kw, $
	                den2_corr=den2_corr_kw, $
	                bicoh_indep=int_bicoh_indep, $
	                indep_den1_corr=den1_corr_indep, $
	                indep_den2_corr=den2_corr_indep, $
	                intr_valid=intr_valid, $
	                corr_flags=corr_flags, $
	                diag_flag=diag_flag, $
	                poisson=poisson_key, $
	                manual_poisson=manual_poisson, $
		                poi_method_used=poisson_method_used, $
		                poi_level_used=poisson_level_used, $
		                diag_half_width=diag_half_width, $
		                diag_bseg_real=diag_bseg_real, $
		                diag_bseg_imag=diag_bseg_imag, $
		                diag_d1seg=diag_d1seg, $
		                diag_d2seg=diag_d2seg, $
		                diag_d1cseg=diag_d1cseg, $
		                diag_d2cseg=diag_d2cseg, $
		                diag_freq=diag_freq, $
		                diag_npix=diag_npix, $
		                percent_step=percent_step,print_every=print_every, $
	                help=help
;+
; NAME:
;      GH_BISPEC2D
;
; PURPOSE:
;      Computes the 2D bispectrum and squared bicoherence from a GHATS FFT file.
;
; EXPLANATION:
;      This routine reads a GHATS FFT file and computes
;
;          B(f1,f2) = < X(f1) X(f2) X*(f1+f2) >
;
;      over a selected rectangular region in the positive-frequency
;      (f1,f2) plane.
;
;      The GHATS FFT files do not include the DC bin:
;
;          rdata[0] = X(df)
;          rdata[1] = X(2df)
;
;      Therefore the array index of f1+f2 is:
;
;          isum = i1 + i2 + 1
;
;      Pixels for which f1+f2 is above the Nyquist frequency are invalid
;      and remain NaN in the derived products.
;
;      Optionally, the raw accumulated sums can be returned through keyword
;      outputs. These are useful for exact 2D rebinning later:
;
;          RAW_BSUM_REAL = real part of sum X1 X2 X3*
;          RAW_BSUM_IMAG = imaginary part of sum X1 X2 X3*
;          RAW_DEN1      = sum |X1 X2|^2
;          RAW_DEN2      = sum |X3|^2
;
;      A later rebinning routine should rebin these sums, not BMOD, BPHASE
;      or BICOH directly.
;
; CALLING SEQUENCE:
;      GH_BISPEC2D,FILENAME,FREQ1,FREQ2, $
;                  BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD_USED,NSEG, $
;                  F1RANGE=F1RANGE,F2RANGE=F2RANGE, $
;                  [,INDEX=INDEX][,TIME=TIME][,SEL=SEL][,RATE=RATE] $
;                  [,RMS=BACK][,/LEAHY] $
;                  [,RAW_BSUM_REAL=RAW_BSUM_REAL] $
;                  [,RAW_BSUM_IMAG=RAW_BSUM_IMAG] $
;                  [,RAW_DEN1=RAW_DEN1] $
;                  [,RAW_DEN2=RAW_DEN2] $
;                  [,PERCENT_STEP=PERCENT_STEP] $
;                  [,PRINT_EVERY=PRINT_EVERY] $
;                  [,/HELP]
;
; INPUTS:
;      FILENAME = name of the input GHATS FFT file
;
; OUTPUTS:
;      FREQUENCY1  = frequency array for the first bispectral axis
;      FREQUENCY2  = frequency array for the second bispectral axis
;      BREAL       = real part of the averaged 2D bispectrum
;      BIMAG       = imaginary part of the averaged 2D bispectrum
;      BMOD        = modulus of the averaged 2D bispectrum
;      BPHASE      = biphase, in radians
;      BICOH       = measured squared bicoherence
;      NPROD_USED  = number of bispectral products contributing to each pixel
;      NSEG        = number of selected FFT segments
;
; KEYWORDS:
;      F1RANGE = two-element frequency range [f1min,f1max] for axis 1.
;                Mandatory.
;
;      F2RANGE = two-element frequency range [f2min,f2max] for axis 2.
;                Mandatory.
;
;      INDEX   = selected FFT index range
;      TIME    = selected time range
;      SEL     = explicit array of selected FFT indices
;      RATE    = selected count-rate range
;
;      RMS     = background rate for rms^3 conversion. If set, BREAL, BIMAG
;                and BMOD are converted to rms^3 units. BPHASE and BICOH are
;                unchanged. Use RMS=0.001 for effectively zero background.
;
;      LEAHY   = convert raw FFT amplitudes to Leahy-amplitude units before
;                forming the bispectrum:
;
;                    X_L = X_raw * sqrt(2/cnts)
;
;                The bispectrum is then in Leahy^(3/2). RMS and LEAHY are
;                mutually exclusive.
;
;      RAW_BSUM_REAL = optional output. Real part of accumulated complex
;                      bispectrum sums before division by NPROD_USED.
;
;      RAW_BSUM_IMAG = optional output. Imaginary part of accumulated complex
;                      bispectrum sums before division by NPROD_USED.
;
;      RAW_DEN1      = optional output. Accumulated bicoherence denominator
;                      term sum |X(f1) X(f2)|^2.
;
;      RAW_DEN2      = optional output. Accumulated bicoherence denominator
;                      term sum |X(f1+f2)|^2.
;
;      INT_BICOH     = optional output. Poisson-denominator-corrected
;                      intrinsic observed squared bicoherence. The complex
;                      bispectrum numerator is unchanged.
;
;      DEN1_CORR     = optional output. Poisson-corrected accumulated
;                      denominator term corresponding to RAW_DEN1.
;
;      DEN2_CORR     = optional output. Poisson-corrected accumulated
;                      denominator term corresponding to RAW_DEN2.
;
;      MANUAL_POISSON = scalar Poisson level in the same power normalization
;                       as |X|^2 after applying /LEAHY or native scaling.
;
;      POISSON=[fmin,fmax].
;                Estimate a constant Poisson level from the selected FFTs'
;                average |X(f)|^2 in this range, excluding the Nyquist bin.
;
;      POISSON = /POISSON. Use Zhang-style frequency-dependent Poisson level.
;                This is mutually exclusive with MANUAL_POISSON and range mode.
;
;      POI_METHOD_USED = optional output string naming the Poisson method used
;                        for INT_BICOH.
;
;      POI_LEVEL_USED  = optional output scalar containing the measured or
;                        supplied scalar/mean Poisson level used for INT_BICOH.
;
;      PERCENT_STEP  = progress-print interval in percent. Default is 10.
;
;      PRINT_EVERY   = progress-print interval in FFTs read. Overrides
;                      PERCENT_STEP if supplied.
;
;      HELP    = print short usage message and return.
;
; EXAMPLES:
;      Basic 2D bispectrum in native FFT units:
;
;      IDL> gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $
;             f1range=[0.1,20.0],f2range=[0.1,20.0]
;
;      Same, converted to rms^3 units:
;
;      IDL> gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $
;             f1range=[0.1,20.0],f2range=[0.1,20.0],rms=0.001
;
;      Same, converted to Leahy^(3/2) units:
;
;      IDL> gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $
;             f1range=[0.1,20.0],f2range=[0.1,20.0],/leahy
;
;      Also return raw accumulated sums for later exact 2D rebinning:
;
;      IDL> gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $
;             f1range=[0.1,20.0],f2range=[0.1,20.0], $
;             raw_bsum_real=bsr,raw_bsum_imag=bsi,raw_den1=den1,raw_den2=den2
;
; COMMON BLOCKS:
;      None
;
; ROUTINES USED:
;      GHATS_OPENFFT
;      GHATS_GETHEADER
;      READ_FFT_LINE
;
; NOTES:
;      This routine computes only the requested frequency rectangle. It does
;      not compute the full Nyquist x Nyquist plane and crop afterwards.
;
;      For exact rebinning, rebin RAW_BSUM_REAL, RAW_BSUM_IMAG, RAW_DEN1,
;      RAW_DEN2 and NPROD_USED, then recompute BREAL, BIMAG, BMOD, BPHASE
;      and BICOH. Do not average BMOD, BPHASE or BICOH directly.
;
; MODIFICATION HISTORY:
;      M. Mendez  21 May 2026  first version, with help from ChatGPT
;      M. Mendez  21 May 2026  added /LEAHY and safer empty-range checks,
;                              with help from ChatGPT
;      M. Mendez  22 May 2026  added optional raw accumulated-sum outputs
;                              for future exact 2D rebinning, with help
;                              from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC2D'
   print,'Compute the 2D bispectrum and squared bicoherence from a GHATS FFT file.'
   print,' '
   print,'Usage:'
   print,"  gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $"
   print,'              f1range=[0.1,20.0],f2range=[0.1,20.0]'
   print,' '
   print,'Example (minimal call):'
   print,"  gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $"
   print,'              f1range=[0.1,20.0],f2range=[0.1,20.0]'
   print,' '
   print,'Example (full call with raw sums):'
   print,"  gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $"
   print,'              f1range=[0.1,20.0],f2range=[0.1,20.0], $'
   print,'              index=[0,99],time=[0,1000],rate=[100,5000], $'
   print,'              /leahy, $'
   print,'              raw_bsum_real=bsr,raw_bsum_imag=bsi, $'
   print,'              raw_den1=den1,raw_den2=den2'
   print,' '
   print,'Optional keywords:'
   print,'  INDEX          selected FFT index range'
   print,'  TIME           selected time range'
   print,'  SEL            selected FFT indices'
   print,'  RATE           selected count-rate range'
   print,'  RMS            background rate for rms^3 conversion'
   print,'  LEAHY          convert FFT amplitudes to Leahy-amplitude units'
   print,'  RAW_BSUM_REAL  output real accumulated bispectrum sum'
   print,'  RAW_BSUM_IMAG  output imaginary accumulated bispectrum sum'
   print,'  RAW_DEN1       output accumulated sum |X1 X2|^2'
   print,'  RAW_DEN2       output accumulated sum |X3|^2'
   print,'  INT_BICOH      output denominator-corrected bicoherence'
   print,'  DEN1_CORR      output corrected sum corresponding to RAW_DEN1'
   print,'  DEN2_CORR      output corrected sum corresponding to RAW_DEN2'
   print,'  MANUAL_POISSON=p        scalar Poisson level in current |X|^2 units'
   print,'  POISSON=[a,b] estimate scalar Poisson level from frequency range'
   print,'  /POISSON       Zhang-style frequency-dependent Poisson estimate'
   print,'  POI_METHOD_USED output Poisson method label'
   print,'  POI_LEVEL_USED  output scalar/mean Poisson level'
   print,'  PERCENT_STEP   progress print interval in percent, default 10'
   print,'  PRINT_EVERY    progress print interval in FFTs read'
   print,' '
   print,'Raw sums are for later exact 2D rebinning. Do not rebin BMOD,'
   print,'BPHASE or BICOH directly.'
   print,' '
   return
endif

;--------------------------------------------------------------------------
; Required frequency ranges.
;--------------------------------------------------------------------------
if(~keyword_set(f1range)) then begin
   massage,'F1RANGE must be specified'
   retall
endif

if(~keyword_set(f2range)) then begin
   massage,'F2RANGE must be specified'
   retall
endif

if(n_elements(f1range) ne 2) then begin
   massage,'F1RANGE must have two elements: [f1min,f1max]'
   retall
endif

if(n_elements(f2range) ne 2) then begin
   massage,'F2RANGE must have two elements: [f2min,f2max]'
   retall
endif

; RMS and LEAHY are two different normalizations. Applying both would be
; ambiguous, so fail explicitly.
if(keyword_set(back) and keyword_set(leahy)) then begin
   massage,'Use either RMS or LEAHY, not both'
   retall
endif

;--------------------------------------------------------------------------
; Optional Poisson correction mode for intrinsic observed bicoherence.
; Explicit user-supplied methods are mutually exclusive. Header Poisson is
; used only if no explicit method is supplied and the header value is valid.
;--------------------------------------------------------------------------
poimode = 0L
n_explicit_poi = 0L

if(n_elements(manual_poisson) gt 0) then begin
   n_explicit_poi = n_explicit_poi + 1L
   poimode = 1L
endif

if(n_elements(poisson_key) eq 2) then begin
   poisson_freq_range = double(poisson_key)
   n_explicit_poi = n_explicit_poi + 1L
   poimode = 2L
endif

if(keyword_set(poisson_key) and n_elements(poisson_key) ne 2) then begin
   n_explicit_poi = n_explicit_poi + 1L
   poimode = 3L
endif

if(n_explicit_poi gt 1L) then begin
   massage,'Choose only one explicit Poisson method: MANUAL_POISSON, POISSON=[range], or /POISSON'
   retall
endif

if(n_elements(poisson_freq_range) gt 0) then begin
   if(n_elements(poisson_freq_range) ne 2) then begin
      massage,'POISSON range must have two elements: [fmin,fmax]'
      retall
   endif
   if(double(poisson_freq_range[1]) lt double(poisson_freq_range[0])) then begin
      massage,'POISSON range invalid: fmax < fmin'
      retall
   endif
endif

;--------------------------------------------------------------------------
; Open FFT file and read the GHATS FFT header.
;--------------------------------------------------------------------------
ghats_openfft,filename,unit,/dialog

ntrafos             = 0l
dummy               = bytarr(100)
gh_version_string   = '                '
observatory         = '                '
instrument          = '                '
target              = '                '
rmjd0               = 0.0D0

ghats_getheader,unit,gh_version_string,observatory,instrument,target,rmjd0, $
                     nft_header,T,ntrafos,e,proliferation,baryflag, $
                     n_spectral_bins,background_flag,dummy

muflag = dummy[0]

;--------------------------------------------------------------------------
; nft_header is the number of time bins in the original light-curve segment.
; The stored FFT has nft_header/2 positive-frequency bins, excluding DC.
;--------------------------------------------------------------------------
nft = nft_header/2
df  = T

allfreq = (findgen(nft)+1.0) * df
fnyq    = allfreq[nft-1]

last_non_nyq = nft - 2L
if(last_non_nyq lt 0L) then begin
   free_lun,unit
   massage,'Frequency array too short to exclude Nyquist bin'
   retall
endif

if(poimode eq 2L) then begin
   fpoi1_user = double(poisson_freq_range[0])
   fpoi2_user = double(poisson_freq_range[1])
   fpoi1 = fpoi1_user > double(allfreq[0])
   fpoi2 = fpoi2_user < double(allfreq[last_non_nyq])

   if(fpoi2 lt fpoi1) then begin
      free_lun,unit
      massage,'POISSON range has no overlap with available non-Nyquist frequencies'
      retall
   endif

   if((fpoi1 ne fpoi1_user) or (fpoi2 ne fpoi2_user)) then begin
      print,'WARNING: POISSON range clipped to available non-Nyquist frequencies.'
      print,'         Requested: ',fpoi1_user,' - ',fpoi2_user,' Hz'
      print,'         Used     : ',fpoi1,' - ',fpoi2,' Hz'
   endif

   wpoi = where((allfreq ge fpoi1) and $
                (allfreq le fpoi2) and $
                (indgen(nft) le last_non_nyq), npoi)

   if(npoi le 0) then begin
      free_lun,unit
      massage,'POISSON range has no valid non-Nyquist Fourier bins'
      retall
   endif
endif

if(poimode eq 3L) then begin
   i_vle = fix([dummy[17:18]],0)+1
   tdead = 1.0d-5
endif

;--------------------------------------------------------------------------
; Select requested frequency ranges. Inclusive limits are used.
;--------------------------------------------------------------------------
w1 = where(allfreq ge f1range[0] and allfreq le f1range[1], nf1)
w2 = where(allfreq ge f2range[0] and allfreq le f2range[1], nf2)

if(nf1 le 0) then begin
   free_lun,unit
   massage,'No Fourier bins found in F1RANGE'
   retall
endif

if(nf2 le 0) then begin
   free_lun,unit
   massage,'No Fourier bins found in F2RANGE'
   retall
endif

frequency1 = allfreq[w1]
frequency2 = allfreq[w2]

;--------------------------------------------------------------------------
; Allocate output and accumulator arrays.
;
; bsum, den1 and den2 are the raw accumulated sums. These are retained so
; that they can optionally be returned for exact 2D rebinning later.
;
; Array convention: [nf1,nf2], first index is f1, second index is f2.
;--------------------------------------------------------------------------
nan = !VALUES.F_NAN

bsum = complexarr(nf1,nf2)
den1 = dblarr(nf1,nf2)
den2 = dblarr(nf1,nf2)

nprod_used = lonarr(nf1,nf2)

breal  = fltarr(nf1,nf2) + nan
bimag  = fltarr(nf1,nf2) + nan
bmod   = fltarr(nf1,nf2) + nan
bphase = fltarr(nf1,nf2) + nan
bicoh  = fltarr(nf1,nf2) + nan
int_bicoh = fltarr(nf1,nf2) + nan
den1_corr = dblarr(nf1,nf2) + !VALUES.D_NAN
den2_corr = dblarr(nf1,nf2) + !VALUES.D_NAN
compute_indep = (arg_present(int_bicoh_indep) or arg_present(den1_corr_indep) or $
                 arg_present(den2_corr_indep) or arg_present(intr_valid) or $
                 arg_present(corr_flags) or arg_present(diag_flag))
if(compute_indep) then begin
   int_bicoh_indep = fltarr(nf1,nf2) + nan
   den1_corr_indep = dblarr(nf1,nf2) + !VALUES.D_NAN
   den2_corr_indep = dblarr(nf1,nf2) + !VALUES.D_NAN
   intr_valid = bytarr(nf1,nf2)
   corr_flags = lonarr(nf1,nf2)
   diag_flag = bytarr(nf1,nf2)
endif

diag_segments = 0B
if(n_elements(diag_half_width) gt 0) then begin
   diag_segments = 1B
   diag_hw = long(diag_half_width)
   if((diag_hw ne 0L) and (diag_hw ne 1L) and $
      (diag_hw ne 2L) and (diag_hw ne 4L)) then begin
      free_lun,unit
      massage,'DIAG_HALF_WIDTH must be 0, 1, 2, or 4'
      retall
   endif
   ndiag = nf1 < nf2
   diag_freq = dblarr(ndiag)
   diag_npix = lonarr(ndiag)
   for idg=0L,ndiag-1L do begin
      diag_freq[idg] = 0.5d0*(double(frequency1[idg]) + double(frequency2[idg]))
      jlo = 0L > (idg - diag_hw)
      jhi = (nf2 - 1L) < (idg + diag_hw)
      npix = 0L
      for jdg=jlo,jhi do begin
         i_native = w1[idg]
         j_native = w2[jdg]
         isum = i_native + j_native + 1L
         if(isum lt nft) then npix = npix + 1L
      endfor
      diag_npix[idg] = npix
   endfor
   diag_bseg_real_full = fltarr(ndiag,ntrafos)
   diag_bseg_imag_full = fltarr(ndiag,ntrafos)
   diag_d1seg_full = dblarr(ndiag,ntrafos)
   diag_d2seg_full = dblarr(ndiag,ntrafos)
   diag_seg_index = lonarr(ntrafos)
endif

power_sum = dblarr(nft)
poisson_header_sum = 0.0d0
n_poisson_header = 0L
poisson_zhang_sum = dblarr(nft)

;--------------------------------------------------------------------------
; Data selection.
;--------------------------------------------------------------------------
if(keyword_set(index)) then begin
   firstfft = index[0]
   lastfft  = index[1]
endif else begin
   firstfft = 0L
   lastfft  = 10000000L
endelse

if(keyword_set(time)) then begin
   firsttime = time[0]
   lasttime  = time[1]
endif else begin
   firsttime = 0.0D0
   lasttime  = 1.0d10
endelse

if(keyword_set(rate)) then begin
   firstrate = rate[0]
   lastrate  = rate[1]
endif else begin
   firstrate = 0.0
   lastrate  = 2000000.0
endelse

if(keyword_set(sel)) then begin
   goodarray = sel
endif else begin
   goodarray = lindgen(ntrafos)
endelse

; Do not index directly with WHERE unless we know the result is non-empty.
wgood = where((goodarray ge firstfft) and (goodarray le lastfft), ngood)

if(ngood le 0) then begin
   free_lun,unit
   massage,'No FFTs found in selected INDEX/SEL range'
   retall
endif

goodarray = goodarray[wgood]
lastindex = max(goodarray)

;--------------------------------------------------------------------------
; Read FFTs and accumulate bispectral sums.
;--------------------------------------------------------------------------
rdata = complexarr(nft)

nseg = 0L
flux = 0.0d0

if(n_elements(print_every) gt 0) then begin
   progress_every = long(print_every)
endif else begin
   progress_every = 0L
endelse

if(n_elements(percent_step) gt 0) then begin
   progress_percent_step = double(percent_step)
endif else begin
   progress_percent_step = 10.0d0
endelse

if(progress_percent_step le 0.0d0) then progress_percent_step = 10.0d0
next_progress_percent = progress_percent_step

print,'Computing 2D bispectrum'
print,'F1RANGE: ',frequency1[0],' - ',frequency1[nf1-1],' Hz (',nf1,' bins)'
print,'F2RANGE: ',frequency2[0],' - ',frequency2[nf2-1],' Hz (',nf2,' bins)'
print,'Total FFTs in file: ',ntrafos

for itrafos=0L,ntrafos-1L do begin

   read_fft_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,rdata

   t_1    = ((rmjd-rmjd0)*86400.0d0)
   t_2    = t_1 + 1.0d0/df
   cratem = cnts * df
   gotcha = where(goodarray eq itrafos)

   if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (cratem ge firstrate) and (cratem le lastrate) and $
      (gotcha[0] ge 0)) then begin

      ; Optional per-segment Leahy-amplitude normalization.
      ; This is applied to the complex FFT amplitudes before forming z.
      if(keyword_set(leahy)) then begin
         if(cnts le 0.0) then begin
            free_lun,unit
            massage,'Cannot apply LEAHY normalization: segment has cnts <= 0'
            retall
         endif
         ampnorm = sqrt(2.0d0/double(cnts))
      endif else begin
         ampnorm = 1.0d0
      endelse

		      iseg_diag = nseg
		      nseg = nseg + 1L
		      if(diag_segments) then diag_seg_index[iseg_diag] = itrafos
		      flux = flux + cnts * df

	      power_norm = double(ampnorm)^2
	      power_sum = power_sum + double(abs(rdata*ampnorm)^2)

	      if(poimode eq 0L) then begin
	         if(finite(poisson) and double(poisson) gt 0.0d0) then begin
	            poisson_header_sum = poisson_header_sum + double(poisson)*power_norm
	            n_poisson_header = n_poisson_header + 1L
	         endif
	      endif

	      if(poimode eq 3L) then begin
	         poisson_estimate,poitmp,cnts,1.0d0/df,current_vle_rate,nft,fndet, $
	                          i_vle,tdead,/differential
	         leahy_power_norm = 2.0d0/double(cnts)
	         poisson_zhang_sum = poisson_zhang_sum + double(poitmp) * $
	                              (power_norm/leahy_power_norm)
	      endif

	      ;-------------------------------------------------------------------
      ; Only the requested frequency rectangle is computed.
      ; Invalid pixels with f1+f2 above Nyquist are skipped.
      ;-------------------------------------------------------------------
      for ii1=0L,nf1-1L do begin

         i1 = w1[ii1]

         for ii2=0L,nf2-1L do begin

            i2 = w2[ii2]

            ; GHATS FFT files exclude DC:
            ; rdata[0]=X(df), rdata[1]=X(2df), therefore isum=i1+i2+1.
            isum = i1 + i2 + 1L

            if(isum lt nft) then begin

               x1 = rdata[i1]    * ampnorm
               x2 = rdata[i2]    * ampnorm
               x3 = rdata[isum] * ampnorm

               z = x1 * x2 * conj(x3)

               ; Raw accumulated complex bispectrum numerator.
               bsum[ii1,ii2] = bsum[ii1,ii2] + z

	               ; Raw accumulated bicoherence denominator terms.
	               den1[ii1,ii2] = den1[ii1,ii2] + double(abs(x1*x2)^2)
	               den2[ii1,ii2] = den2[ii1,ii2] + double(abs(x3)^2)

	               nprod_used[ii1,ii2] = nprod_used[ii1,ii2] + 1L

		               if(diag_segments) then begin
		                  if((ii1 lt ndiag) and $
		                     (abs(long(ii1)-long(ii2)) le diag_hw)) then begin
		                     diag_bseg_real_full[ii1,iseg_diag] = $
		                        diag_bseg_real_full[ii1,iseg_diag] + float(z)
		                     diag_bseg_imag_full[ii1,iseg_diag] = $
		                        diag_bseg_imag_full[ii1,iseg_diag] + imaginary(z)
		                     diag_d1seg_full[ii1,iseg_diag] = $
		                        diag_d1seg_full[ii1,iseg_diag] + double(abs(x1*x2)^2)
		                     diag_d2seg_full[ii1,iseg_diag] = $
		                        diag_d2seg_full[ii1,iseg_diag] + double(abs(x3)^2)
		                  endif
		               endif

	            endif

         endfor
      endfor

   endif

   done_count = itrafos + 1L
   progress_percent = 100.0d0 * double(done_count) / double(ntrafos)

   if(progress_every gt 0L) then begin
      if((done_count mod progress_every) eq 0L or done_count eq ntrafos) then begin
         print,'  progress: ',strtrim(string(done_count),2),' / ',strtrim(string(ntrafos),2),' FFTs read (',strtrim(string(progress_percent,format='(F6.1)'),2),'%), selected ',strtrim(string(nseg),2)
      endif
   endif else begin
      if(progress_percent ge next_progress_percent or done_count eq ntrafos) then begin
         print,'  progress: ',strtrim(string(progress_percent,format='(F6.1)'),2),'% (',strtrim(string(done_count),2),' / ',strtrim(string(ntrafos),2),' FFTs read), selected ',strtrim(string(nseg),2)
         next_progress_percent = next_progress_percent + progress_percent_step
      endif
   endelse

   if(itrafos gt (lastindex-1)) then goto,finished_reading

endfor

finished_reading:

free_lun,unit

if(nseg le 0) then begin
   massage,'No FFTs retrieved for input selection'
   retall
endif

flux = flux / nseg

print,'  ',strtrim(string(nseg),1),' FFTs selected'

power_avg = power_sum / double(nseg)
poisson_power = dblarr(nft) + !VALUES.D_NAN
compute_intobs = 0B
poisson_method_used = 'NONE'
poisson_level_used = !VALUES.D_NAN

case poimode of
   1L: begin
      poisson_power = double(manual_poisson) + dblarr(nft)
      poisson_method_used = 'MANUAL_POISSON'
      poisson_level_used = double(manual_poisson)
      compute_intobs = 1B
   end

   2L: begin
      poisson_level = total(power_avg[wpoi],/double)/double(npoi)
      poisson_power = poisson_level + dblarr(nft)
      poisson_method_used = 'POISSON_RANGE'
      poisson_level_used = poisson_level
      compute_intobs = 1B
      print,'Using POISSON=[range] for INT_BICOH: ',fpoi1,' - ',fpoi2,' Hz'
      print,'  estimated Poisson level: ',poisson_level
   end

   3L: begin
      poisson_power = poisson_zhang_sum / double(nseg)
      poisson_method_used = 'ZHANG'
      wpn = where(finite(poisson_power),npn)
      if(npn gt 0) then poisson_level_used = total(poisson_power[wpn],/double)/double(npn)
      compute_intobs = 1B
      print,'Using /POISSON Zhang-style denominator correction for INT_BICOH.'
   end

   else: begin
      if(n_poisson_header gt 0L) then begin
         poisson_level = poisson_header_sum / double(n_poisson_header)
         poisson_power = poisson_level + dblarr(nft)
         poisson_method_used = 'HEADER_POISSON'
         poisson_level_used = poisson_level
         compute_intobs = 1B
         print,'Using header Poisson level for INT_BICOH: ',poisson_level
      endif else begin
         print,'No Poisson correction available; INT_BICOH will be NaN.'
      endelse
   end
endcase

source_power = power_avg - poisson_power

if(compute_indep and (compute_intobs eq 0)) then begin
   massage,'Diagonal-aware independent-noise auto correction requested but no Poisson/noise correction is available'
   retall
endif

if(diag_segments) then begin
   diag_bseg_real = diag_bseg_real_full[*,0:nseg-1L]
   diag_bseg_imag = diag_bseg_imag_full[*,0:nseg-1L]
   diag_d1seg = diag_d1seg_full[*,0:nseg-1L]
   diag_d2seg = diag_d2seg_full[*,0:nseg-1L]
   diag_d1cseg = dblarr(ndiag,nseg) + !VALUES.D_NAN
   diag_d2cseg = dblarr(ndiag,nseg) + !VALUES.D_NAN

   if(compute_intobs) then begin
      for idg=0L,ndiag-1L do begin
         jlo = 0L > (idg - diag_hw)
         jhi = (nf2 - 1L) < (idg + diag_hw)
         corr1_strip = 0.0d0
         corr2_strip = 0.0d0
         for jdg=jlo,jhi do begin
            i1 = w1[idg]
            i2 = w2[jdg]
            isum = i1 + i2 + 1L
            if(isum lt nft) then begin
               pn1 = poisson_power[i1]
               pn2 = poisson_power[i2]
               ps1 = source_power[i1]
               ps2 = source_power[i2]
               ps3 = source_power[isum]
               if(finite(pn1) and finite(pn2) and finite(ps1) and $
                  finite(ps2) and finite(ps3)) then begin
                  corr1_strip = corr1_strip + (pn1*ps2 + pn2*ps1 + pn1*pn2)
                  corr2_strip = corr2_strip + ps3
               endif
            endif
         endfor
         for isg=0L,nseg-1L do begin
            diag_d1cseg[idg,isg] = diag_d1seg[idg,isg] - corr1_strip
            diag_d2cseg[idg,isg] = corr2_strip
         endfor
      endfor
   endif
endif

;--------------------------------------------------------------------------
; Average bispectrum and compute derived quantities.
;
; B = bsum / nprod_used
;
; bicoherence is computed from the raw sums, not from the averaged B.
;--------------------------------------------------------------------------
for ii1=0L,nf1-1L do begin
   for ii2=0L,nf2-1L do begin

      if(nprod_used[ii1,ii2] gt 0L) then begin

         bavg = bsum[ii1,ii2] / double(nprod_used[ii1,ii2])

         breal[ii1,ii2]  = float(bavg)
         bimag[ii1,ii2]  = imaginary(bavg)
         bmod[ii1,ii2]   = abs(bavg)
         bphase[ii1,ii2] = atan(imaginary(bavg),float(bavg))

	         if((den1[ii1,ii2] gt 0.0d0) and (den2[ii1,ii2] gt 0.0d0)) then begin
	            bicoh[ii1,ii2] = float(abs(bsum[ii1,ii2])^2 / $
	                              (den1[ii1,ii2]*den2[ii1,ii2]))
	         endif else begin
	            bicoh[ii1,ii2] = nan
	         endelse

	         if(compute_intobs) then begin
	            i1 = w1[ii1]
	            i2 = w2[ii2]
	            isum = i1 + i2 + 1L
	            if(isum lt nft) then begin
	               pn1 = poisson_power[i1]
	               pn2 = poisson_power[i2]
	               ps1 = source_power[i1]
	               ps2 = source_power[i2]
	               ps3 = source_power[isum]

	               if(finite(pn1) and finite(pn2) and finite(ps1) and $
	                  finite(ps2) and finite(ps3)) then begin
	                  nprod_d = double(nprod_used[ii1,ii2])
	                  d1c = den1[ii1,ii2] - nprod_d*(pn1*ps2 + pn2*ps1 + pn1*pn2)
	                  d2c = nprod_d*ps3
	                  den1_corr[ii1,ii2] = d1c
	                  den2_corr[ii1,ii2] = d2c
		                  if((d1c gt 0.0d0) and (d2c gt 0.0d0)) then begin
		                     int_bicoh[ii1,ii2] = float(abs(bsum[ii1,ii2])^2 / (d1c*d2c))
		                  endif else begin
		                     int_bicoh[ii1,ii2] = nan
		                  endelse
		               endif
		            endif
		         endif

		         if(compute_indep) then begin
		            i1 = w1[ii1]
		            i2 = w2[ii2]
		            isum = i1 + i2 + 1L
		            if(isum lt nft) then begin
		               diagonal = (i1 eq i2)
		               d1ci = gh_bispec_den1_corr_indep(den1[ii1,ii2], $
		                         nprod_used[ii1,ii2],source_power[i1], $
		                         poisson_power[i1],source_power[i2], $
		                         poisson_power[i2],1B,diagonal,flag=flag1)
		               d2ci = gh_bispec_den2_corr_indep(nprod_used[ii1,ii2], $
		                         source_power[isum],flag=flag2)
		               den1_corr_indep[ii1,ii2] = d1ci
		               den2_corr_indep[ii1,ii2] = d2ci
		               cflag = flag1
		               if((flag2 and 4L) ne 0L) then cflag = cflag + 4L
		               if(((flag2 and 8L) ne 0L) and ((cflag and 8L) eq 0L)) then cflag = cflag + 8L
		               corr_flags[ii1,ii2] = cflag
		               if(diagonal) then diag_flag[ii1,ii2] = 1B
		               if((d1ci gt 0.0d0) and (d2ci gt 0.0d0)) then begin
		                  int_bicoh_indep[ii1,ii2] = float(abs(bsum[ii1,ii2])^2 / (d1ci*d2ci))
		                  intr_valid[ii1,ii2] = 1B
		               endif
		            endif
		         endif

		      endif

	   endfor
endfor

;--------------------------------------------------------------------------
; Optional rms^3 conversion.
;
; Important: if RMS is used, the derived BREAL/BIMAG/BMOD are converted
; after the sums are accumulated. The raw sums returned by RAW_* remain the
; raw accumulated sums in the normalization used internally above.
;--------------------------------------------------------------------------
if(keyword_set(back)) then begin

   if(back ge flux) then begin
      massage,'Error: background flux higher than source+bkg flux!'
      retall
   endif

   norm_bispec = (2.0d0^1.5d0) * (df^1.5d0) / ((flux-back)^3)

   breal = breal * norm_bispec
   bimag = bimag * norm_bispec
   bmod  = bmod  * norm_bispec

   print,'Bispectrum units: rms^3'

endif else begin

   if(keyword_set(leahy)) then begin
      print,'Bispectrum units: Leahy^(3/2)'
   endif else begin
      print,'Bispectrum units: native FFT units'
   endelse

endelse

;--------------------------------------------------------------------------
; Optional raw accumulated-sum outputs for future exact 2D rebinning.
;
; These are deliberately returned after the derived products are computed,
; but they are the accumulated sums, not averaged BREAL/BIMAG.
;--------------------------------------------------------------------------
raw_bsum_real = float(bsum)
raw_bsum_imag = imaginary(bsum)
raw_den1      = den1
raw_den2      = den2
den1_corr_kw = den1_corr
den2_corr_kw = den2_corr

end
