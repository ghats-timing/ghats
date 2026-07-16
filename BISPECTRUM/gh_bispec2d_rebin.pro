;==========================================================================
; Helper: create 1D rebin groups
;==========================================================================
pro gh_bispec2d_make_bins,frequency,irf,bin1,bin2,freq_reb

nf = n_elements(frequency)

if(nf le 0) then begin
   bin1 = -1L
   bin2 = -1L
   freq_reb = !values.d_nan
   return
endif

if(~keyword_set(irf)) then irf = 1
if(irf eq 0) then irf = 1

;--------------------------------------------------------------------------
; Linear rebinning.
; Positive IRF means combine IRF adjacent native bins.
;--------------------------------------------------------------------------
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

   ;-----------------------------------------------------------------------
   ; Logarithmic rebinning.
   ; Negative IRF gives approximately constant fractional bin width:
   ;
   ;     f_hi/f_lo = 10^(1/abs(IRF))
   ;
   ; Thus IRF=-100 means 100 bins per decade.
   ;
   ;-----------------------------------------------------------------------
   factor = 10.0d0^(1.0d0/abs(double(irf)))

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
pro gh_bispec2d_rebin,frequency1,frequency2, $
                       raw_bsum_real,raw_bsum_imag,raw_den1,raw_den2,nprod_used, $
                       freq1_reb,freq2_reb, $
                       breal,bimag,bmod,bphase,bicoh,nprod_reb, $
                       irf1=irf1,irf2=irf2, $
                       f1range=f1range,f2range=f2range, $
                       raw_bsum_real_reb=raw_bsum_real_reb, $
                       raw_bsum_imag_reb=raw_bsum_imag_reb, $
                       raw_den1_reb=raw_den1_reb, $
                       raw_den2_reb=raw_den2_reb, $
                       d1corr=d1corr, $
                       d2corr=d2corr, $
                       den1_corr_reb=den1_corr_reb, $
                       den2_corr_reb=den2_corr_reb, $
                       int_bicoh=int_bicoh, $
                       help=help
;+
; NAME:
;      GH_BISPEC2D_REBIN
;
; PURPOSE:
;      Rebins a 2D bispectrum exactly from raw accumulated sums.
;
; EXPLANATION:
;      This routine rebins a 2D bispectrum using the raw accumulated
;      quantities returned by GH_BISPEC2D:
;
;          RAW_BSUM_REAL
;          RAW_BSUM_IMAG
;          RAW_DEN1
;          RAW_DEN2
;          NPROD_USED
;
;      It does NOT average BMOD, BPHASE or BICOH. Instead it first combines
;      the raw complex sums and bicoherence denominators, then recomputes:
;
;          BREAL
;          BIMAG
;          BMOD
;          BPHASE
;          BICOH
;
;      This is the correct way to rebin the bispectrum.
;
; CALLING SEQUENCE:
;      GH_BISPEC2D_REBIN,F1,F2, $
;                         RAW_BSUM_REAL,RAW_BSUM_IMAG,RAW_DEN1,RAW_DEN2,NPROD, $
;                         F1_REB,F2_REB, $
;                         BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD_REB, $
;                         [,IRF1=IRF1][,IRF2=IRF2] $
;                         [,F1RANGE=F1RANGE][,F2RANGE=F2RANGE][,/HELP]
;
; INPUTS:
;      FREQUENCY1     = native f1 axis
;      FREQUENCY2     = native f2 axis
;      RAW_BSUM_REAL  = real part of accumulated bispectrum sums
;      RAW_BSUM_IMAG  = imaginary part of accumulated bispectrum sums
;      RAW_DEN1       = accumulated sum |X(f1) X(f2)|^2
;      RAW_DEN2       = accumulated sum |X(f1+f2)|^2
;      NPROD_USED     = number of products contributing to each native pixel
;
; OUTPUTS:
;      FREQ1_REB   = rebinned f1 axis
;      FREQ2_REB   = rebinned f2 axis
;      BREAL       = real part of rebinned bispectrum
;      BIMAG       = imaginary part of rebinned bispectrum
;      BMOD        = modulus of rebinned bispectrum
;      BPHASE      = biphase, radians
;      BICOH       = rebinned squared bicoherence
;      NPROD_REB   = number of products contributing to each rebinned pixel
;
; KEYWORDS:
;      IRF1    = rebin factor for f1.
;                Positive = linear rebinning.
;                Negative = logarithmic rebinning, with abs(IRF1)
;                bins per decade.
;
;      IRF2    = rebin factor for f2.
;                Positive = linear rebinning.
;                Negative = logarithmic rebinning, with abs(IRF2)
;                bins per decade.
;
;      F1RANGE = optional f1 range to keep before rebinning.
;
;      F2RANGE = optional f2 range to keep before rebinning.
;
;      RAW_BSUM_REAL_REB = optional output rebinned real raw sum.
;
;      RAW_BSUM_IMAG_REB = optional output rebinned imaginary raw sum.
;
;      RAW_DEN1_REB      = optional output rebinned denominator term 1.
;
;      RAW_DEN2_REB      = optional output rebinned denominator term 2.
;
;      D1CORR            = optional input corrected denominator term 1.
;                          This deliberately does not begin with DEN1_CORR
;                          to avoid IDL keyword-abbreviation ambiguity with
;                          DEN1_CORR_REB.
;
;      D2CORR            = optional input corrected denominator term 2.
;
;      DEN1_CORR_REB     = optional output rebinned corrected denominator 1.
;
;      DEN2_CORR_REB     = optional output rebinned corrected denominator 2.
;
;      INT_BICOH         = optional output rebinned intrinsic observed
;                          squared bicoherence, computed from rebinned raw
;                          numerator and rebinned corrected denominators.
;
;      HELP    = print this usage message.
;
; EXAMPLES:
;      Linear 4 x 4 rebinning:
;
;      IDL> gh_bispec2d_rebin,f1,f2,bsr,bsi,den1,den2,nprod, $
;             f1r,f2r,breal,bimag,bmod,bphase,bicoh,nprod_r,irf1=4,irf2=4
;
;      Logarithmic rebinning on both axes:
;
;      IDL> gh_bispec2d_rebin,f1,f2,bsr,bsi,den1,den2,nprod, $
;             f1r,f2r,breal,bimag,bmod,bphase,bicoh,nprod_r,irf1=-20,irf2=-20
;
;      Different rebinning on each axis:
;
;      IDL> gh_bispec2d_rebin,f1,f2,bsr,bsi,den1,den2,nprod, $
;             f1r,f2r,breal,bimag,bmod,bphase,bicoh,nprod_r, $
;             irf1=4,irf2=-20
;
; NOTES:
;      This routine assumes the input arrays follow the GH_BISPEC2D
;      convention:
;
;          array[f1_index,f2_index]
;
;      Invalid native pixels should have NPROD_USED=0. They are ignored.
;
; MODIFICATION HISTORY:
;      M. Mendez  22 May 2026  first version, with help from ChatGPT
;      M. Mendez  22 May 2026  changed logarithmic rebinning to use
;                              10^(1/abs(IRF)) and removed IRF shortcut
;                              keyword to avoid ambiguity, with help
;                              from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC2D_REBIN'
   print,'Rebin a 2D bispectrum exactly from raw accumulated sums.'
   print,' '
   print,'Usage:'
   print,'  gh_bispec2d_rebin,f1,f2,bsr,bsi,den1,den2,nprod, $'
   print,'                     f1r,f2r,breal,bimag,bmod,bphase,bicoh,nprod_r, $'
   print,'                     irf1=4,irf2=4'
   print,' '
   print,'Examples:'
   print,'  irf1=-100,irf2=-100  same log rebinning on f1 and f2'
   print,'  irf1=20,irf2=20      same lin rebinning on f1 and f2'
   print,'  irf1=20,irf2=-20     different rebinning on f1 (lin) and f2 (log)'
   print,' '
   print,'Optional ranges:'
   print,'  f1range=[0.1,20.0],f2range=[0.1,20.0]'
   print,' '
   print,'Optional complete-FITS outputs:'
   print,'  raw_bsum_real_reb=, raw_bsum_imag_reb=, raw_den1_reb=, raw_den2_reb='
   print,'  d1corr=, d2corr=, den1_corr_reb=, den2_corr_reb=, int_bicoh='
   print,' '
   print,'Important: this routine rebins raw sums, not BMOD/BPHASE/BICOH.'
   print,' '
   return
endif

;--------------------------------------------------------------------------
; Basic input checks.
;--------------------------------------------------------------------------
nf1 = n_elements(frequency1)
nf2 = n_elements(frequency2)

if(nf1 le 0 or nf2 le 0) then begin
   massage,'Input frequency arrays are empty'
   retall
endif

s = size(raw_bsum_real)

if(s[0] ne 2) then begin
   massage,'RAW_BSUM_REAL must be a 2D array'
   retall
endif

if((s[1] ne nf1) or (s[2] ne nf2)) then begin
   massage,'Input array dimensions do not match frequency axes'
   retall
endif

; IDL allows abbreviated keyword names. Keep the input corrected-denominator
; keyword names distinct from DEN1_CORR_REB/DEN2_CORR_REB to avoid ambiguity.
have_corr = (n_elements(d1corr) gt 0) or (n_elements(d2corr) gt 0)
if(have_corr) then begin
   if((n_elements(d1corr) eq 0) or (n_elements(d2corr) eq 0)) then begin
      massage,'Give both D1CORR and D2CORR'
      retall
   endif
   sc1 = size(d1corr)
   sc2 = size(d2corr)
   if((sc1[0] ne 2) or (sc1[1] ne nf1) or (sc1[2] ne nf2) or $
      (sc2[0] ne 2) or (sc2[1] ne nf1) or (sc2[2] ne nf2)) then begin
      massage,'Corrected denominator arrays do not match frequency axes'
      retall
   endif
endif

;--------------------------------------------------------------------------
; Choose rebin factors.
; IRF1 controls the f1 axis.
; IRF2 controls the f2 axis.
; Positive values give linear rebinning.
; Negative values give logarithmic rebinning, with abs(IRF) bins per decade.
;--------------------------------------------------------------------------
if(keyword_set(irf1)) then begin
   rbf1 = irf1
endif else begin
   rbf1 = 1
endelse

if(keyword_set(irf2)) then begin
   rbf2 = irf2
endif else begin
   rbf2 = 1
endelse

if(rbf1 eq 0) then rbf1 = 1
if(rbf2 eq 0) then rbf2 = 1

;--------------------------------------------------------------------------
; Select optional f1/f2 ranges before rebinning.
;--------------------------------------------------------------------------
if(keyword_set(f1range)) then begin

   if(n_elements(f1range) ne 2) then begin
      massage,'F1RANGE must have two elements: [f1min,f1max]'
      retall
   endif

   w1 = where(frequency1 ge f1range[0] and frequency1 le f1range[1], nsel1)

   if(nsel1 le 0) then begin
      massage,'No bins found in F1RANGE'
      retall
   endif

endif else begin

   w1 = lindgen(nf1)
   nsel1 = nf1

endelse

if(keyword_set(f2range)) then begin

   if(n_elements(f2range) ne 2) then begin
      massage,'F2RANGE must have two elements: [f2min,f2max]'
      retall
   endif

   w2 = where(frequency2 ge f2range[0] and frequency2 le f2range[1], nsel2)

   if(nsel2 le 0) then begin
      massage,'No bins found in F2RANGE'
      retall
   endif

endif else begin

   w2 = lindgen(nf2)
   nsel2 = nf2

endelse

freq1_sel = frequency1[w1]
freq2_sel = frequency2[w2]

;--------------------------------------------------------------------------
; Build rebin groups on the selected axes.
;--------------------------------------------------------------------------
gh_bispec2d_make_bins,freq1_sel,rbf1,bin1a,bin1b,freq1_reb
gh_bispec2d_make_bins,freq2_sel,rbf2,bin2a,bin2b,freq2_reb

nrb1 = n_elements(freq1_reb)
nrb2 = n_elements(freq2_reb)

;--------------------------------------------------------------------------
; Allocate outputs.
;--------------------------------------------------------------------------
nan = !values.f_nan

breal     = dblarr(nrb1,nrb2) + nan
bimag     = dblarr(nrb1,nrb2) + nan
bmod      = dblarr(nrb1,nrb2) + nan
bphase    = fltarr(nrb1,nrb2) + !values.f_nan
bicoh     = fltarr(nrb1,nrb2) + !values.f_nan
nprod_reb = lonarr(nrb1,nrb2)
raw_bsum_real_reb = dblarr(nrb1,nrb2) + !values.d_nan
raw_bsum_imag_reb = dblarr(nrb1,nrb2) + !values.d_nan
raw_den1_reb = dblarr(nrb1,nrb2) + !values.d_nan
raw_den2_reb = dblarr(nrb1,nrb2) + !values.d_nan
if(have_corr) then begin
   den1_corr_reb = dblarr(nrb1,nrb2) + !values.d_nan
   den2_corr_reb = dblarr(nrb1,nrb2) + !values.d_nan
   int_bicoh = fltarr(nrb1,nrb2) + !values.f_nan
endif

;--------------------------------------------------------------------------
; Main 2D rebinning.
;
; For each output pixel, combine all valid native pixels in the corresponding
; f1 x f2 block:
;
;     BSUM_reb = sum BSUM_i
;     DEN1_reb = sum DEN1_i
;     DEN2_reb = sum DEN2_i
;     N_reb    = sum N_i
;
; Then recompute:
;
;     B_reb    = BSUM_reb / N_reb
;     bicoherence = |BSUM_reb|^2 / (DEN1_reb DEN2_reb)
;
;--------------------------------------------------------------------------
for iout=0L,nrb1-1L do begin

   i1a = bin1a[iout]
   i1b = bin1b[iout]

   native_i1 = w1[i1a:i1b]

   for jout=0L,nrb2-1L do begin

      i2a = bin2a[jout]
      i2b = bin2b[jout]

      native_i2 = w2[i2a:i2b]

      bsum_re = 0.0d0
      bsum_im = 0.0d0
      den1_sum = 0.0d0
      den2_sum = 0.0d0
      den1c_sum = 0.0d0
      den2c_sum = 0.0d0
      nsum = 0L

      for ii=0L,n_elements(native_i1)-1L do begin

         i_native = native_i1[ii]

         for jj=0L,n_elements(native_i2)-1L do begin

            j_native = native_i2[jj]

            npix = nprod_used[i_native,j_native]

            if(npix gt 0L and $
               finite(raw_bsum_real[i_native,j_native]) and $
               finite(raw_bsum_imag[i_native,j_native]) and $
               finite(raw_den1[i_native,j_native]) and $
               finite(raw_den2[i_native,j_native])) then begin

               bsum_re = bsum_re + double(raw_bsum_real[i_native,j_native])
               bsum_im = bsum_im + double(raw_bsum_imag[i_native,j_native])

               den1_sum = den1_sum + double(raw_den1[i_native,j_native])
               den2_sum = den2_sum + double(raw_den2[i_native,j_native])

               if(have_corr) then begin
                  if(finite(d1corr[i_native,j_native]) and $
                     finite(d2corr[i_native,j_native])) then begin
                     den1c_sum = den1c_sum + double(d1corr[i_native,j_native])
                     den2c_sum = den2c_sum + double(d2corr[i_native,j_native])
                  endif
               endif

               nsum = nsum + long(npix)

            endif

         endfor

      endfor

      if(nsum gt 0L) then begin

         nprod_reb[iout,jout] = nsum
         raw_bsum_real_reb[iout,jout] = bsum_re
         raw_bsum_imag_reb[iout,jout] = bsum_im
         raw_den1_reb[iout,jout] = den1_sum
         raw_den2_reb[iout,jout] = den2_sum
         if(have_corr) then begin
            den1_corr_reb[iout,jout] = den1c_sum
            den2_corr_reb[iout,jout] = den2c_sum
         endif

         bre = bsum_re / double(nsum)
         bim = bsum_im / double(nsum)

         breal[iout,jout]  = bre
         bimag[iout,jout]  = bim
         bmod[iout,jout]   = sqrt(bre^2 + bim^2)
         bphase[iout,jout] = atan(bim,bre)

         if((den1_sum gt 0.0d0) and (den2_sum gt 0.0d0)) then begin
            bicoh[iout,jout] = float((bsum_re^2 + bsum_im^2) / $
                               (den1_sum*den2_sum))
         endif

         if(have_corr) then begin
            if((den1c_sum gt 0.0d0) and (den2c_sum gt 0.0d0)) then begin
               int_bicoh[iout,jout] = float((bsum_re^2 + bsum_im^2) / $
                                      (den1c_sum*den2c_sum))
            endif
         endif

      endif

   endfor

endfor

print,'Rebinned 2D bispectrum'
print,'F1 bins: ',nf1,' -> ',nrb1
print,'F2 bins: ',nf2,' -> ',nrb2
print,'IRF1: ',rbf1
print,'IRF2: ',rbf2

end
