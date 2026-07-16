;+
; NAME:
;      GH_BISPEC_CORR_DEN
;
; PURPOSE:
;      Shared scalar helpers for source/noise-corrected bicoherence
;      denominator accumulators.
;
; EXPLANATION:
;      These routines deliberately operate on one native pixel at a time.
;      Rebinning code should first compute corrected native-pixel
;      accumulators, then sum those corrected accumulators before deriving
;      intrinsic bicoherence.
;
;      For independent-noise cross products,
;
;          B_123(f1,f2) = < X1(f1) X2(f2) X3*(f1+f2) >
;
;      with X=s+n and independent zero-mean Fourier noise, the off-diagonal
;      source term in RAW_DEN1=sum |X1 X2|^2 is estimated by subtracting
;
;          Pn1 Ps2 + Pn2 Ps1 + Pn1 Pn2 .
;
;      If X1 and X2 are the same band and f1=f2, RAW_DEN1=sum |X|^4.
;      Under circular complex Gaussian noise, without assuming Gaussian
;      source variability, subtract
;
;          4 Ps Pn + 2 Pn^2
;
;      to estimate the source fourth moment E|s|^4.
;
;      Flags returned by GH_BISPEC_DEN1_CORR_INDEP:
;          1 = same-band diagonal fourth-moment correction used
;          2 = corrected DEN1 is non-positive
;          4 = corrected DEN2 is non-positive (set by caller if desired)
;          8 = invalid or unavailable input
;-

function gh_bispec_den1_corr_indep,raw_den1,nprod,ps1,pn1,ps2,pn2, $
                                      same12,diagonal,flag=flag

flag = 0L
nan = !values.d_nan

if((finite(raw_den1) eq 0) or (finite(ps1) eq 0) or $
   (finite(pn1) eq 0) or (finite(ps2) eq 0) or $
   (finite(pn2) eq 0) or long(nprod) le 0L) then begin
   flag = 8L
   return,nan
endif

np = double(nprod)

if(keyword_set(same12) and keyword_set(diagonal)) then begin
   corr = double(raw_den1) - np*(4.0d0*double(ps1)*double(pn1) + $
                                 2.0d0*double(pn1)^2)
   flag = flag + 1L
endif else begin
   corr = double(raw_den1) - np*(double(pn1)*double(ps2) + $
                                 double(pn2)*double(ps1) + $
                                 double(pn1)*double(pn2))
endelse

if(finite(corr) eq 0) then begin
   flag = flag + 8L
   return,nan
endif

if(corr le 0.0d0) then flag = flag + 2L
return,corr

end


function gh_bispec_den2_corr_indep,nprod,ps3,flag=flag

flag = 0L
nan = !values.d_nan

if((finite(ps3) eq 0) or long(nprod) le 0L) then begin
   flag = 8L
   return,nan
endif

corr = double(nprod)*double(ps3)
if(finite(corr) eq 0) then begin
   flag = 8L
   return,nan
endif

if(corr le 0.0d0) then flag = flag + 4L
return,corr

end
