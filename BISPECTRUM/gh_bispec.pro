;==========================================================================
; Helper: create rebin groups for GH_BISPEC
;==========================================================================
pro gh_bispec_make_bins,frequency,irf,bin1,bin2,freq_reb

nf = n_elements(frequency)

if(nf le 0) then begin
   bin1 = -1L
   bin2 = -1L
   freq_reb = !values.f_nan
   return
endif

if(irf eq 0) then irf = 1

if(irf gt 0) then begin

   nb = long(ceil(float(nf)/float(irf)))
   bin1 = lonarr(nb)
   bin2 = lonarr(nb)
   freq_reb = dblarr(nb)

   for ib=0L,nb-1L do begin
      i1 = ib*long(irf)
      i2 = ((ib+1L)*long(irf)-1L) < (nf-1L)
      bin1[ib] = i1
      bin2[ib] = i2
      freq_reb[ib] = mean(frequency[i1:i2])
   endfor

endif else begin

   factor = exp(1.0d0/abs(double(irf)))

   b1 = lonarr(nf)
   b2 = lonarr(nf)
   fr = dblarr(nf)

   ib = 0L
   i1 = 0L

   while(i1 lt nf) do begin

      fhi = double(frequency[i1]) * factor
      w = where(frequency ge frequency[i1] and frequency lt fhi,nw)

      if(nw gt 0) then begin
         i2 = w[nw-1]
      endif else begin
         i2 = i1
      endelse

      if(i2 lt i1) then i2 = i1
      if(i2 ge nf) then i2 = nf-1L

      b1[ib] = i1
      b2[ib] = i2
      fr[ib] = mean(frequency[i1:i2])

      ib = ib + 1L
      i1 = i2 + 1L

   endwhile

   bin1 = b1[0:ib-1L]
   bin2 = b2[0:ib-1L]
   freq_reb = fr[0:ib-1L]

endelse

end


;==========================================================================
; Main routine
;==========================================================================
pro gh_bispec,filename,frequency, $
              breal,breal_err,bimag,bimag_err, $
              bmod,bmod_err,bphase,bphase_err, $
              bicoh,nprod_used,nseg, $
              f1range=f1range,f2range=f2range, $
              irf=irf,index=index,time=time,sel=sel,rate=rate, $
              rms=back,leahy=leahy,hfcorr=hfcorr,modfloor=modfloor, $
              help=help,xhelp=xhelp
;+
; NAME:
;      GH_BISPEC
;
; PURPOSE:
;      Computes a 1D ridge/strip of the bispectrum and bicoherence from
;      a GHATS FFT file.
;
; EXPLANATION:
;      This routine computes
;
;          B(f1,f2) = < X(f1) X(f2) X*(f1+f2) >
;
;      averaged over all f1 bins in F1RANGE, as a function of f2.
;
;      This is intended for extracting the bispectrum along a QPO ridge:
;      F1RANGE should usually be the FWHM interval of the QPO.
;
; CALLING SEQUENCE:
;      GH_BISPEC,FILENAME,FREQ, $
;                BREAL,BREAL_ERR,BIMAG,BIMAG_ERR, $
;                BMOD,BMOD_ERR,BPHASE,BPHASE_ERR, $
;                BICOH,NPROD_USED,NSEG, $
;                F1RANGE=F1RANGE,F2RANGE=F2RANGE, $
;                [,IRF=IRF][,INDEX=INDEX][,TIME=TIME][,SEL=SEL] $
;                [,RATE=RATE][,RMS=BACK][,/LEAHY] $
;                [,HFCORR=HFCORR][,MODFLOOR=MODFLOOR] $
;                [,/HELP][,/XHELP]
;
; KEYWORDS:
;      F1RANGE  = two-element range [f1min,f1max]. Mandatory.
;                 These f1 bins are averaged over.
;
;      F2RANGE  = two-element range [f2min,f2max]. Mandatory.
;
;      IRF      = rebin factor in f2. Positive = linear rebinning.
;                 Negative = logarithmic rebinning. Default is 1.
;
;      RMS      = background rate for rms^3 conversion.
;
;      LEAHY    = if set, convert raw FFT amplitudes to Leahy-amplitude
;                 units before forming the bispectrum. LEAHY and RMS are
;                 mutually exclusive.
;
;      HFCORR   = two-element range [fmin,fmax]. If set, subtract the
;                 complex high-frequency mean from BREAL/BIMAG and recompute
;                 BMOD/BPHASE. This is not a Poisson-noise subtraction.
;
;      MODFLOOR = two-element range [fmin,fmax]. If set, subtract the mean
;                 BMOD in that range from BMOD only. This is an empirical
;                 modulus-floor removal for comparison plots, not a Poisson
;                 correction.
;
;      XHELP    = print extended notes about HFCORR and MODFLOOR.
;
; EXAMPLE:
;      IDL> gh_bispec,'file.fft',freq, $
;             breal,breal_err,bimag,bimag_err, $
;             bmod,bmod_err,bphase,bphase_err, $
;             bicoh,nprod,nseg, $
;             f1range=[0.8,1.2],f2range=[0.1,500.0], $
;             irf=-20,rms=0.001
;
;      IDL> gh_bispec,'file.fft',freq, $
;             breal,breal_err,bimag,bimag_err, $
;             bmod,bmod_err,bphase,bphase_err, $
;             bicoh,nprod,nseg, $
;             f1range=[0.8,1.2],f2range=[0.1,500.0], $
;             hfcorr=[100,500]
;
; MODIFICATION HISTORY:
;      M. Mendez  21 May 2026  first version, developed from the GHATS
;                              cross-spectrum routines, with help from ChatGPT
;      M. Mendez  21 May 2026  added /LEAHY, HFCORR and MODFLOOR options,
;                              with help from ChatGPT
;      M. Mendez  22 May 2026  cleaned /HELP and /XHELP text, fixed LEAHY
;                              block, added safer IRF=0 handling and
;                              comments throughout, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC'
   print,'Compute a 1D ridge/strip of the bispectrum from a GHATS FFT file.'
   print,' '
   print,'Usage:'
   print,"  gh_bispec,'file.fft',freq,breal,breal_err,bimag,bimag_err, $"
   print,'            bmod,bmod_err,bphase,bphase_err,bicoh,nprod,nseg, $'
   print,'            f1range=[0.8,1.2],f2range=[0.1,500],irf=-20,rms=0.001'
   print,' '

   print,'Example (minimal call):'
   print,"  gh_bispec,'file.fft',freq,breal,breal_err,bimag,bimag_err, $"
   print,'            bmod,bmod_err,bphase,bphase_err,bicoh,nprod,nseg, $'
   print,'            f1range=[0.8,1.2],f2range=[0.1,500]'
   print,' '

   print,'Example (/LEAHY):'
   print,"  gh_bispec,'file.fft',freq,breal,breal_err,bimag,bimag_err, $"
   print,'            bmod,bmod_err,bphase,bphase_err,bicoh,nprod,nseg, $'
   print,'            f1range=[0.8,1.2],f2range=[0.1,500],/leahy'
   print,' '

   print,'Example (empirical corrections):'
   print,"  gh_bispec,'file.fft',freq,breal,breal_err,bimag,bimag_err, $"
   print,'            bmod,bmod_err,bphase,bphase_err,bicoh,nprod,nseg, $'
   print,'            f1range=[0.8,1.2],f2range=[0.1,500], $'
   print,'            hfcorr=[100,500],modfloor=[100,500]'
   print,' '

   print,'Main keywords:'
   print,'  F1RANGE   frequency interval averaged over, normally QPO FWHM'
   print,'  F2RANGE   output frequency interval'
   print,'  IRF       rebin factor (positive=linear, negative=log)'
   print,'  RMS       convert to rms^3 units'
   print,'  LEAHY     convert to Leahy^(3/2) units'
   print,'  HFCORR    subtract complex high-frequency baseline'
   print,'  MODFLOOR  subtract empirical BMOD floor'
   print,' '

   print,'IMPORTANT: Use /XHELP before using HFCORR or MODFLOOR.'
   print,' '

   return
endif

if(keyword_set(xhelp)) then begin

   print,' '
   print,'GH_BISPEC extended notes'
   print,' '

   print,'HFCORR=[fmin,fmax]:'
   print,'  Subtracts mean(BREAL)+i mean(BIMAG) measured over the range.'
   print,'  BMOD and BPHASE are then recomputed.'
   print,'  This is a linear complex baseline correction.'
   print,'  It is NOT a Poisson-noise subtraction.'
   print,' '

   print,'MODFLOOR=[fmin,fmax]:'
   print,'  Subtracts mean(BMOD) measured over the range from BMOD only.'
   print,'  This can be useful for visual comparison with a PDS.'
   print,'  However, it is only an empirical modulus-floor removal.'
   print,'  It may not remove the full positive-definite noise bias'
   print,'  in BMOD and can produce negative values.'
   print,' '

   print,'IMPORTANT:'
   print,' '
   print,'  None of this is equivalent to the Poisson correction in a PDS.'
   print,' '
   print,'  HFCORR subtracts any bias in the complex bispectrum to make'
   print,'  the complex average approximately zero at high frequencies.'
   print,'  However, since BMOD is positive definite, this does NOT'
   print,'  generally eliminate any bias in BMOD.'
   print,' '
   print,'  MODFLOOR may help compare the SHAPE of BMOD to a PDS,'
   print,'  but it is not equivalent to subtracting Poisson noise.'
   print,' '

   return
endif

;--------------------------------------------------------------------------
; Required inputs and keyword consistency checks.
;--------------------------------------------------------------------------

; F1RANGE is mandatory.
if(~keyword_set(f1range)) then begin
   massage,'F1RANGE must be specified'
   retall
endif

; F2RANGE is mandatory.
if(~keyword_set(f2range)) then begin
   massage,'F2RANGE must be specified'
   retall
endif

; Require valid two-element ranges.
if(n_elements(f1range) ne 2) then begin
   massage,'F1RANGE must have two elements: [f1min,f1max]'
   retall
endif

if(n_elements(f2range) ne 2) then begin
   massage,'F2RANGE must have two elements: [f2min,f2max]'
   retall
endif

; RMS and LEAHY are mutually exclusive.
if(keyword_set(back) and keyword_set(leahy)) then begin
   massage,'Use either RMS or LEAHY, not both'
   retall
endif

; Validate optional HFCORR range.
if(keyword_set(hfcorr)) then begin
   if(n_elements(hfcorr) ne 2) then begin
      massage,'HFCORR must have two elements: [fmin,fmax]'
      retall
   endif
endif

; Validate optional MODFLOOR range.
if(keyword_set(modfloor)) then begin
   if(n_elements(modfloor) ne 2) then begin
      massage,'MODFLOOR must have two elements: [fmin,fmax]'
      retall
   endif
endif

; IRF=0 is not meaningful. Default to no rebinning.
if(~keyword_set(irf)) then irf = 1
if(irf eq 0) then irf = 1

;--------------------------------------------------------------------------
; Open FFT file and read header
;--------------------------------------------------------------------------
ghats_openfft,filename,unit,/dialog

ntrafos             = 0L
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

nft = nft_header/2
df  = T

allfreq = (findgen(nft)+1.0) * df

;--------------------------------------------------------------------------
; Select f1 and f2 bins with inclusive limits.
;--------------------------------------------------------------------------
wf1 = where(allfreq ge f1range[0] and allfreq le f1range[1], nf1)
wf2 = where(allfreq ge f2range[0] and allfreq le f2range[1], nf2_native)

if(nf1 le 0) then begin
   free_lun,unit
   massage,'No Fourier bins found in F1RANGE'
   retall
endif

if(nf2_native le 0) then begin
   free_lun,unit
   massage,'No Fourier bins found in F2RANGE'
   retall
endif

freq_native = allfreq[wf2]

;--------------------------------------------------------------------------
; Define output/rebinned f2 bins.
;--------------------------------------------------------------------------
gh_bispec_make_bins,freq_native,irf,bin1,bin2,frequency

nfb = n_elements(frequency)

;--------------------------------------------------------------------------
; Selection of FFT segments.
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
   lasttime  = 1.0D10
endelse

if(keyword_set(rate)) then begin
   firstrate = rate[0]
   lastrate  = rate[1]
endif else begin
   firstrate = 0.0
   lastrate  = 2000000.0
endelse

if(keyword_set(sel)) then begin

   ; User supplied an explicit list of FFT indices.
   goodarray = sel

endif else begin

   ; Default: all FFTs in the file are initially allowed.
   goodarray = lindgen(ntrafos)

endelse

; Apply INDEX limits safely. Do not index directly with WHERE unless
; the result is known to be non-empty.
wgood = where((goodarray ge firstfft) and (goodarray le lastfft), ngood)

if(ngood le 0) then begin
   free_lun,unit
   massage,'No FFTs found in selected INDEX/SEL range'
   retall
endif

goodarray = goodarray[wgood]
lastindex = max(goodarray)

;--------------------------------------------------------------------------
; Allocate accumulators.
;--------------------------------------------------------------------------
rdata = complexarr(nft)

; Segment-by-segment rebinned complex bispectrum.
; Used only for estimating errors from segment scatter.
seg_reb = complexarr(nfb,ntrafos)

; Pooled raw sums used for the bicoherence normalization.
; These are not segment-averaged quantities.
bsum_reb = complexarr(nfb)
den1_reb = dblarr(nfb)
den2_reb = dblarr(nfb)

nprod_used = lonarr(nfb)

nseg = 0L
flux = 0.0D0

print,'Computing bispectrum ridge'
print,'F1RANGE: ',allfreq[wf1[0]],' - ',allfreq[wf1[nf1-1]],' Hz (',nf1,' bins)'
print,'F2RANGE: ',freq_native[0],' - ',freq_native[nf2_native-1],' Hz (',nf2_native,' bins)'
print,'Rebin factor: ',irf

;--------------------------------------------------------------------------
; Main loop over FFT segments.
;--------------------------------------------------------------------------
for itrafos=0L,ntrafos-1L do begin

   read_fft_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,rdata

   t_1    = ((rmjd-rmjd0)*86400.0D0)
   t_2    = t_1 + 1.0D0/df
   cratem = cnts * df
   gotcha = where(goodarray eq itrafos)

   if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (cratem ge firstrate) and (cratem le lastrate) and $
      (gotcha[0] ge 0)) then begin

      ; Optional per-segment Leahy-amplitude normalization.
      ; This must be applied before forming the third-order product.
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

      strip_native = complexarr(nf2_native)
      nvalid_native = lonarr(nf2_native)

      for j2=0L,nf2_native-1L do begin

         i2 = wf2[j2]

         zsum = complex(0.0,0.0)
         nv = 0L

         for j1=0L,nf1-1L do begin

            i1 = wf1[j1]
            isum = i1 + i2 + 1L

            if(isum lt nft) then begin

               x1 = rdata[i1]    * ampnorm
               x2 = rdata[i2]    * ampnorm
               x3 = rdata[isum] * ampnorm

               z = x1 * x2 * conj(x3)

               zsum = zsum + z
               nv = nv + 1L

            endif

         endfor

         if(nv gt 0L) then begin
            strip_native[j2] = zsum / double(nv)
            nvalid_native[j2] = nv
         endif

      endfor

      ; Rebin this segment's complex strip.
      for ib=0L,nfb-1L do begin

         i1b = bin1[ib]
         i2b = bin2[ib]

         wv = where(nvalid_native[i1b:i2b] gt 0L,nwv)

         if(nwv gt 0) then begin
            wv = wv + i1b
            seg_reb[ib,nseg] = total(strip_native[wv]) / double(nwv)
         endif

      endfor

      ; Accumulate bicoherence sums using all native products.
      for ib=0L,nfb-1L do begin

         i1b = bin1[ib]
         i2b = bin2[ib]

         for j2=i1b,i2b do begin

            i2 = wf2[j2]

            for j1=0L,nf1-1L do begin

               i1 = wf1[j1]
               isum = i1 + i2 + 1L

               if(isum lt nft) then begin

                  x1 = rdata[i1]    * ampnorm
                  x2 = rdata[i2]    * ampnorm
                  x3 = rdata[isum] * ampnorm

                  z = x1 * x2 * conj(x3)

                  bsum_reb[ib] = bsum_reb[ib] + z
                  den1_reb[ib] = den1_reb[ib] + double(abs(x1*x2)^2)
                  den2_reb[ib] = den2_reb[ib] + double(abs(x3)^2)

                  nprod_used[ib] = nprod_used[ib] + 1L

               endif

            endfor

         endfor

      endfor

      nseg = nseg + 1L
      flux = flux + cnts * df

   endif

   if(itrafos gt (lastindex-1)) then goto,finished_reading

endfor

finished_reading:

free_lun,unit

if(nseg lt 2) then begin
   massage,'At least two FFT segments are required to compute errors'
   retall
endif

flux = flux / nseg

print,'  ',strtrim(string(nseg),1),' FFTs selected'

;--------------------------------------------------------------------------
; Compute averages and errors from segment scatter.
;--------------------------------------------------------------------------
nan = !values.f_nan

breal      = fltarr(nfb) + nan
bimag      = fltarr(nfb) + nan
breal_err  = fltarr(nfb) + nan
bimag_err  = fltarr(nfb) + nan
bmod       = fltarr(nfb) + nan
bmod_err   = fltarr(nfb) + nan
bphase     = fltarr(nfb) + nan
bphase_err = fltarr(nfb) + nan
bicoh      = fltarr(nfb) + nan

for ib=0L,nfb-1L do begin

   vals = seg_reb[ib,0:nseg-1L]

   xr = float(vals)
   xi = imaginary(vals)

   breal[ib] = mean(xr)
   bimag[ib] = mean(xi)

   if(nseg gt 1) then begin
      breal_err[ib] = stddev(xr) / sqrt(float(nseg))
      bimag_err[ib] = stddev(xi) / sqrt(float(nseg))
   endif

   bmod[ib] = sqrt(breal[ib]^2 + bimag[ib]^2)
   bphase[ib] = atan(bimag[ib],breal[ib])

   if(bmod[ib] gt 0.0) then begin
      bmod_err[ib] = sqrt((breal[ib]/bmod[ib])^2 * breal_err[ib]^2 + $
                          (bimag[ib]/bmod[ib])^2 * bimag_err[ib]^2)

      bphase_err[ib] = sqrt((bimag[ib]^2 * breal_err[ib]^2 + $
                             breal[ib]^2 * bimag_err[ib]^2) / bmod[ib]^4)
   endif

   if((den1_reb[ib] gt 0.0D0) and (den2_reb[ib] gt 0.0D0)) then begin
      bicoh[ib] = float(abs(bsum_reb[ib])^2 / (den1_reb[ib]*den2_reb[ib]))
   endif

endfor

;--------------------------------------------------------------------------
; Optional rms^3 conversion.
;--------------------------------------------------------------------------
if(keyword_set(back)) then begin

   if(back ge flux) then begin
      massage,'Error: background flux higher than source+bkg flux!'
      retall
   endif

   norm_bispec = (2.0D0^1.5D0) * (df^1.5D0) / ((flux-back)^3)

   breal     = breal     * norm_bispec
   bimag     = bimag     * norm_bispec
   bmod      = bmod      * norm_bispec
   breal_err = breal_err * norm_bispec
   bimag_err = bimag_err * norm_bispec
   bmod_err  = bmod_err  * norm_bispec

   print,'Bispectrum units: rms^3'

endif else begin

   if(keyword_set(leahy)) then begin
      print,'Bispectrum units: Leahy^(3/2)'
   endif else begin
      print,'Bispectrum units: native FFT units'
   endelse

endelse

;--------------------------------------------------------------------------
; Optional complex high-frequency baseline correction.
; This is not a Poisson-noise subtraction.
;--------------------------------------------------------------------------
if(keyword_set(hfcorr)) then begin

   whf = where(frequency ge hfcorr[0] and frequency le hfcorr[1], nhf)

   if(nhf le 0) then begin
      massage,'No frequency bins found in HFCORR range'
      retall
   endif

   ; Estimate and subtract a complex offset.
   ; This can remove a real/imaginary baseline, but not the positive-definite
   ; finite-sample floor in BMOD.
   breal0 = mean(breal[whf])
   bimag0 = mean(bimag[whf])

   breal = breal - breal0
   bimag = bimag - bimag0

   for ib=0L,nfb-1L do begin

      bmod[ib] = sqrt(breal[ib]^2 + bimag[ib]^2)
      bphase[ib] = atan(bimag[ib],breal[ib])

      if(bmod[ib] gt 0.0) then begin
         bmod_err[ib] = sqrt((breal[ib]/bmod[ib])^2 * breal_err[ib]^2 + $
                             (bimag[ib]/bmod[ib])^2 * bimag_err[ib]^2)

         bphase_err[ib] = sqrt((bimag[ib]^2 * breal_err[ib]^2 + $
                                breal[ib]^2 * bimag_err[ib]^2) / bmod[ib]^4)
      endif else begin
         bmod_err[ib] = nan
         bphase_err[ib] = nan
      endelse

   endfor

   print,'Subtracted complex HF baseline over ',frequency[whf[0]], $
         ' - ',frequency[whf[nhf-1]],' Hz'
   print,'  Re baseline: ',breal0
   print,'  Im baseline: ',bimag0
   print,'  This is not a Poisson-noise subtraction.'

endif

;--------------------------------------------------------------------------
; Optional empirical modulus-floor removal.
; This is for comparison plots, not a statistically clean correction.
;--------------------------------------------------------------------------
if(keyword_set(modfloor)) then begin

   wmf = where(frequency ge modfloor[0] and frequency le modfloor[1], nmf)

   if(nmf le 0) then begin
      massage,'No frequency bins found in MODFLOOR range'
      retall
   endif

   mod0 = mean(bmod[wmf])
   bmod = bmod - mod0

   print,'Subtracted empirical BMOD floor over ',frequency[wmf[0]], $
         ' - ',frequency[wmf[nmf-1]],' Hz'
   print,'  BMOD floor: ',mod0
   print,'  This is not a Poisson-noise subtraction.'
   print,'  Negative BMOD values can occur after this empirical subtraction.'

endif

end

