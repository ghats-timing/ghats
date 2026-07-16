;==========================================================================
; Helper: make 1D rebin groups for GH_BISPEC_HYPOT
;==========================================================================
pro gh_bispec_hypot_make_bins,frequency,irf,bin1,bin2,freq_reb

nf = n_elements(frequency)

if(nf le 0) then begin
   bin1 = -1L
   bin2 = -1L
   freq_reb = !values.d_nan
   return
endif

if(~keyword_set(irf)) then irf = 1
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
pro gh_bispec_hypot,frequency1,frequency2, $
                    raw_bsum_real,raw_bsum_imag,raw_den1,raw_den2,nprod_used, $
                    f0,width, $
                    freq_out, $
                    breal,bimag,bmod,bphase,bicoh,nprod,npix, $
                    breal_avg,bimag_avg,bmod_avg,bphase_avg,bicoh_avg, $
                    nprod_avg,npix_avg, $
                    xtype=xtype,sym=sym,irf=irf, $
                    f1range=f1range,f2range=f2range,mask=mask, $
                    all_f1=all_f1,all_f2=all_f2,all_freq=all_freq, $
                    all_breal=all_breal,all_bimag=all_bimag, $
                    all_bmod=all_bmod,all_bphase=all_bphase, $
                    all_bicoh=all_bicoh,all_nprod=all_nprod, $
                    help=help
;+
; NAME:
;      GH_BISPEC_HYPOT
;
; PURPOSE:
;      Extracts bispectral products along a hypotenuse f1+f2=f0.
;
; EXPLANATION:
;      This routine extracts a finite-width strip around
;
;          f1 + f2 = f0
;
;      from a 2D bispectrum, using the raw accumulated sums returned by
;      GH_BISPEC2D:
;
;          RAW_BSUM_REAL
;          RAW_BSUM_IMAG
;          RAW_DEN1
;          RAW_DEN2
;          NPROD_USED
;
;      The routine does not average BMOD, BPHASE or BICOH directly.
;      Instead, it combines the raw complex sums and denominators, then
;      recomputes:
;
;          BREAL
;          BIMAG
;          BMOD
;          BPHASE
;          BICOH
;
;      This is the correct way to collapse pixels along the hypotenuse.
;
; CALLING SEQUENCE:
;      GH_BISPEC_HYPOT,F1,F2, $
;                      RAW_BSUM_REAL,RAW_BSUM_IMAG,RAW_DEN1,RAW_DEN2,NPROD2D, $
;                      F0,WIDTH, $
;                      FREQ_OUT, $
;                      BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD,NPIX, $
;                      BREAL_AVG,BIMAG_AVG,BMOD_AVG,BPHASE_AVG,BICOH_AVG, $
;                      NPROD_AVG,NPIX_AVG, $
;                      [,XTYPE=XTYPE][,SYM=SYM][,IRF=IRF] $
;                      [,F1RANGE=F1RANGE][,F2RANGE=F2RANGE][,MASK=MASK] $
;                      [,ALL_F1=ALL_F1][,ALL_F2=ALL_F2][,ALL_FREQ=ALL_FREQ] $
;                      [,ALL_BREAL=ALL_BREAL][,ALL_BIMAG=ALL_BIMAG] $
;                      [,ALL_BMOD=ALL_BMOD][,ALL_BPHASE=ALL_BPHASE] $
;                      [,ALL_BICOH=ALL_BICOH][,ALL_NPROD=ALL_NPROD] $
;                      [,/HELP]
;
; INPUTS:
;      FREQUENCY1      = f1 axis of the 2D bispectrum
;      FREQUENCY2      = f2 axis of the 2D bispectrum
;      RAW_BSUM_REAL   = real part of accumulated complex bispectrum sums
;      RAW_BSUM_IMAG   = imaginary part of accumulated complex bispectrum sums
;      RAW_DEN1        = accumulated sum |X(f1)X(f2)|^2
;      RAW_DEN2        = accumulated sum |X(f1+f2)|^2
;      NPROD_USED      = number of products in each 2D pixel
;      F0              = hypotenuse sum frequency, f1+f2=f0
;      WIDTH           = full width of the strip in Hz.
;                        Pixels satisfy abs(f1+f2-f0) <= WIDTH/2.
;
; OUTPUTS:
;      FREQ_OUT        = output frequency coordinate of the extracted curve
;      BREAL           = real part along the hypotenuse
;      BIMAG           = imaginary part along the hypotenuse
;      BMOD            = modulus along the hypotenuse
;      BPHASE          = biphase along the hypotenuse, radians
;      BICOH           = squared bicoherence along the hypotenuse
;      NPROD           = total number of products per output bin
;      NPIX            = number of 2D pixels per output bin
;
;      BREAL_AVG       = real part averaged over the whole selected hypotenuse
;      BIMAG_AVG       = imaginary part averaged over the whole selected hypotenuse
;      BMOD_AVG        = modulus averaged over the whole selected hypotenuse
;      BPHASE_AVG      = biphase averaged over the whole selected hypotenuse
;      BICOH_AVG       = squared bicoherence over the whole selected hypotenuse
;      NPROD_AVG       = total number of products in the whole selected strip
;      NPIX_AVG        = total number of 2D pixels in the whole selected strip
;
; KEYWORDS:
;      XTYPE           = coordinate used for FREQ_OUT.
;                        'F1'  = f1, default
;                        'F2'  = f2
;                        'SUM' = f1+f2
;                        'MIN' = min(f1,f2)
;
;      SYM             = symmetry selection.
;                        'UPPER' = keep f2 >= f1, default
;                        'LOWER' = keep f2 <= f1
;                        'BOTH'  = keep both halves
;
;      IRF             = optional 1D rebin factor after extraction.
;                        Positive = linear rebinning.
;                        Negative = logarithmic rebinning, with abs(IRF)
;                        bins per decade.
;                        Default is 1.
;
;      F1RANGE         = optional f1 range [f1min,f1max]
;      F2RANGE         = optional f2 range [f2min,f2max]
;
;      MASK            = optional 2D mask. MASK=1 keeps a pixel, MASK=0 rejects.
;
;      ALL_*           = optional selected-pixel outputs for debugging.
;
; NOTES:
;      Empty output bins are kept and filled with NaN, with NPROD=0 and NPIX=0.
;
;      This V1 routine does not compute statistically rigorous errors.
;      The available raw sums are accumulated over segments, so segment-scatter
;      errors cannot be recovered here.
;
; MODIFICATION HISTORY:
;      M. Mendez  22 May 2026  first version, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC_HYPOT'
   print,'Extract bispectrum products along f1+f2=f0.'
   print,' '
   print,'Usage:'
   print,'  gh_bispec_hypot,f1,f2,bsr,bsi,den1,den2,nprod2d, $'
   print,'                  f0,width,freq_out, $'
   print,'                  breal,bimag,bmod,bphase,bicoh,nprod,npix, $'
   print,'                  breal_avg,bimag_avg,bmod_avg,bphase_avg,bicoh_avg, $'
   print,'                  nprod_avg,npix_avg'
   print,' '
   print,'Main keywords:'
   print,"  xtype='F1'|'F2'|'SUM'|'MIN'    default 'F1'"
   print,"  sym='UPPER'|'LOWER'|'BOTH'      default 'UPPER'"
   print,'  irf=IRF                         optional 1D rebinning'
   print,'  f1range=[fmin,fmax]'
   print,'  f2range=[fmin,fmax]'
   print,'  mask=mask2d'
   print,' '
   print,'Hypotenuse condition: abs(f1+f2-f0) <= width/2'
   print,' '
   return
endif

;--------------------------------------------------------------------------
; Basic checks.
;--------------------------------------------------------------------------
nf1 = n_elements(frequency1)
nf2 = n_elements(frequency2)

if(nf1 le 0 or nf2 le 0) then begin
   massage,'Input frequency arrays are empty'
   retall
endif

if(n_params() lt 24) then begin
   massage,'Usage: gh_bispec_hypot,f1,f2,bsr,bsi,den1,den2,nprod2d,f0,width,...'
   retall
endif

if(width le 0.0d0) then begin
   massage,'WIDTH must be positive'
   retall
endif

s = size(raw_bsum_real)

if(s[0] ne 2) then begin
   massage,'RAW_BSUM_REAL must be a 2D array'
   retall
endif

if((s[1] ne nf1) or (s[2] ne nf2)) then begin
   massage,'RAW_BSUM_REAL dimensions do not match frequency axes'
   retall
endif

arrays_ok = 1

sb = size(raw_bsum_imag)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(raw_den1)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(raw_den2)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(nprod_used)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0

if(arrays_ok eq 0) then begin
   massage,'One or more raw arrays do not match frequency axes'
   retall
endif

if(keyword_set(mask)) then begin
   sm = size(mask)
   if((sm[0] ne 2) or (sm[1] ne nf1) or (sm[2] ne nf2)) then begin
      massage,'MASK must have same dimensions as raw arrays'
      retall
   endif
endif

if(keyword_set(f1range)) then begin
   if(n_elements(f1range) ne 2) then begin
      massage,'F1RANGE must be [f1min,f1max]'
      retall
   endif
endif

if(keyword_set(f2range)) then begin
   if(n_elements(f2range) ne 2) then begin
      massage,'F2RANGE must be [f2min,f2max]'
      retall
   endif
endif

;--------------------------------------------------------------------------
; Defaults.
;--------------------------------------------------------------------------
if(keyword_set(xtype)) then begin
   xt = strupcase(strtrim(xtype,2))
endif else begin
   xt = 'F1'
endelse

if((xt ne 'F1') and (xt ne 'F2') and (xt ne 'SUM') and (xt ne 'MIN')) then begin
   massage,"XTYPE must be 'F1', 'F2', 'SUM' or 'MIN'"
   retall
endif

if(keyword_set(sym)) then begin
   sy = strupcase(strtrim(sym,2))
endif else begin
   sy = 'UPPER'
endelse

if((sy ne 'UPPER') and (sy ne 'LOWER') and (sy ne 'BOTH')) then begin
   massage,"SYM must be 'UPPER', 'LOWER' or 'BOTH'"
   retall
endif

if(~keyword_set(irf)) then irf = 1
if(irf eq 0) then irf = 1

halfwidth = 0.5d0 * double(width)

;--------------------------------------------------------------------------
; First pass: count selected pixels.
;--------------------------------------------------------------------------
nsel = 0L

for i=0L,nf1-1L do begin

   f1v = double(frequency1[i])

   if(keyword_set(f1range)) then begin
      if((f1v lt f1range[0]) or (f1v gt f1range[1])) then goto,next_i_count
   endif

   for j=0L,nf2-1L do begin

      f2v = double(frequency2[j])

      if(keyword_set(f2range)) then begin
         if((f2v lt f2range[0]) or (f2v gt f2range[1])) then goto,next_j_count
      endif

      if(sy eq 'UPPER') then begin
         if(f2v lt f1v) then goto,next_j_count
      endif

      if(sy eq 'LOWER') then begin
         if(f2v gt f1v) then goto,next_j_count
      endif

      if(abs(f1v + f2v - double(f0)) gt halfwidth) then goto,next_j_count

      if(nprod_used[i,j] le 0L) then goto,next_j_count

      if(keyword_set(mask)) then begin
         if(mask[i,j] eq 0) then goto,next_j_count
      endif

      if(~finite(raw_bsum_real[i,j])) then goto,next_j_count
      if(~finite(raw_bsum_imag[i,j])) then goto,next_j_count
      if(~finite(raw_den1[i,j])) then goto,next_j_count
      if(~finite(raw_den2[i,j])) then goto,next_j_count

      nsel = nsel + 1L

      next_j_count:

   endfor

   next_i_count:

endfor

nan = !values.f_nan

if(nsel le 0L) then begin

   print,'No valid hypotenuse pixels found for f0=',f0,' width=',width

   freq_out = [!values.d_nan]
   breal    = [nan]
   bimag    = [nan]
   bmod     = [nan]
   bphase   = [nan]
   bicoh    = [nan]
   nprod    = [0L]
   npix     = [0L]

   breal_avg  = nan
   bimag_avg  = nan
   bmod_avg   = nan
   bphase_avg = nan
   bicoh_avg  = nan
   nprod_avg  = 0L
   npix_avg   = 0L

   return

endif

;--------------------------------------------------------------------------
; Second pass: store selected pixels.
;--------------------------------------------------------------------------
pix_f1     = dblarr(nsel)
pix_f2     = dblarr(nsel)
pix_freq   = dblarr(nsel)
pix_bsr    = dblarr(nsel)
pix_bsi    = dblarr(nsel)
pix_den1   = dblarr(nsel)
pix_den2   = dblarr(nsel)
pix_nprod  = lonarr(nsel)

k = 0L

for i=0L,nf1-1L do begin

   f1v = double(frequency1[i])

   if(keyword_set(f1range)) then begin
      if((f1v lt f1range[0]) or (f1v gt f1range[1])) then goto,next_i_store
   endif

   for j=0L,nf2-1L do begin

      f2v = double(frequency2[j])

      if(keyword_set(f2range)) then begin
         if((f2v lt f2range[0]) or (f2v gt f2range[1])) then goto,next_j_store
      endif

      if(sy eq 'UPPER') then begin
         if(f2v lt f1v) then goto,next_j_store
      endif

      if(sy eq 'LOWER') then begin
         if(f2v gt f1v) then goto,next_j_store
      endif

      if(abs(f1v + f2v - double(f0)) gt halfwidth) then goto,next_j_store

      if(nprod_used[i,j] le 0L) then goto,next_j_store

      if(keyword_set(mask)) then begin
         if(mask[i,j] eq 0) then goto,next_j_store
      endif

      if(~finite(raw_bsum_real[i,j])) then goto,next_j_store
      if(~finite(raw_bsum_imag[i,j])) then goto,next_j_store
      if(~finite(raw_den1[i,j])) then goto,next_j_store
      if(~finite(raw_den2[i,j])) then goto,next_j_store

      pix_f1[k] = f1v
      pix_f2[k] = f2v

      case xt of
         'F1':  pix_freq[k] = f1v
         'F2':  pix_freq[k] = f2v
         'SUM': pix_freq[k] = f1v + f2v
         'MIN': pix_freq[k] = f1v < f2v
      endcase

      pix_bsr[k]   = double(raw_bsum_real[i,j])
      pix_bsi[k]   = double(raw_bsum_imag[i,j])
      pix_den1[k]  = double(raw_den1[i,j])
      pix_den2[k]  = double(raw_den2[i,j])
      pix_nprod[k] = long(nprod_used[i,j])

      k = k + 1L

      next_j_store:

   endfor

   next_i_store:

endfor

;--------------------------------------------------------------------------
; Optional all-pixel outputs.
;--------------------------------------------------------------------------
all_f1    = pix_f1
all_f2    = pix_f2
all_freq  = pix_freq
all_nprod = pix_nprod

all_breal  = dblarr(nsel) + nan
all_bimag  = dblarr(nsel) + nan
all_bmod   = dblarr(nsel) + nan
all_bphase = dblarr(nsel) + nan
all_bicoh  = dblarr(nsel) + nan

for k=0L,nsel-1L do begin

   if(pix_nprod[k] gt 0L) then begin

      bre = pix_bsr[k] / double(pix_nprod[k])
      bim = pix_bsi[k] / double(pix_nprod[k])

      all_breal[k]  = bre
      all_bimag[k]  = bim
      all_bmod[k]   = sqrt(bre^2 + bim^2)
      all_bphase[k] = atan(bim,bre)

      if((pix_den1[k] gt 0.0d0) and (pix_den2[k] gt 0.0d0)) then begin
         all_bicoh[k] = float((pix_bsr[k]^2 + pix_bsi[k]^2) / $
                        (pix_den1[k]*pix_den2[k]))
      endif

   endif

endfor

;--------------------------------------------------------------------------
; Sort selected pixels by output frequency.
;--------------------------------------------------------------------------
isort = sort(pix_freq)

pix_f1    = pix_f1[isort]
pix_f2    = pix_f2[isort]
pix_freq  = pix_freq[isort]
pix_bsr   = pix_bsr[isort]
pix_bsi   = pix_bsi[isort]
pix_den1  = pix_den1[isort]
pix_den2  = pix_den2[isort]
pix_nprod = pix_nprod[isort]

;--------------------------------------------------------------------------
; Collapse pixels with identical output frequency.
;--------------------------------------------------------------------------
tmp_freq  = dblarr(nsel)
tmp_breal = dblarr(nsel) + nan
tmp_bimag = dblarr(nsel) + nan
tmp_bmod  = dblarr(nsel) + nan
tmp_bpha  = dblarr(nsel) + nan
tmp_bicoh = dblarr(nsel) + nan
tmp_nprod = lonarr(nsel)
tmp_npix  = lonarr(nsel)

nbin0 = 0L
k0 = 0L

while(k0 lt nsel) do begin

   kval = pix_freq[k0]

   k1 = k0
   while((k1+1L) lt nsel) do begin
      if(pix_freq[k1+1L] eq kval) then begin
         k1 = k1 + 1L
      endif else begin
         goto,done_group
      endelse
   endwhile

   done_group:

   bsr_sum = total(pix_bsr[k0:k1],/double)
   bsi_sum = total(pix_bsi[k0:k1],/double)
   d1_sum  = total(pix_den1[k0:k1],/double)
   d2_sum  = total(pix_den2[k0:k1],/double)
   np_sum  = long(total(pix_nprod[k0:k1],/double))
   nx_sum  = k1 - k0 + 1L

   tmp_freq[nbin0]  = kval
   tmp_nprod[nbin0] = np_sum
   tmp_npix[nbin0]  = nx_sum

   if(np_sum gt 0L) then begin

      bre = bsr_sum / double(np_sum)
      bim = bsi_sum / double(np_sum)

      tmp_breal[nbin0] = bre
      tmp_bimag[nbin0] = bim
      tmp_bmod[nbin0]  = sqrt(bre^2 + bim^2)
      tmp_bpha[nbin0]  = atan(bim,bre)

      if((d1_sum gt 0.0d0) and (d2_sum gt 0.0d0)) then begin
         tmp_bicoh[nbin0] = float((bsr_sum^2 + bsi_sum^2) / (d1_sum*d2_sum))
      endif

   endif

   nbin0 = nbin0 + 1L
   k0 = k1 + 1L

endwhile

tmp_freq  = tmp_freq[0:nbin0-1L]
tmp_breal = tmp_breal[0:nbin0-1L]
tmp_bimag = tmp_bimag[0:nbin0-1L]
tmp_bmod  = tmp_bmod[0:nbin0-1L]
tmp_bpha  = tmp_bpha[0:nbin0-1L]
tmp_bicoh = tmp_bicoh[0:nbin0-1L]
tmp_nprod = tmp_nprod[0:nbin0-1L]
tmp_npix  = tmp_npix[0:nbin0-1L]

;--------------------------------------------------------------------------
; Optional 1D rebinning of the collapsed curve.
; To keep the rebinning exact, use the pixel-level raw sums again.
;--------------------------------------------------------------------------
if(irf ne 1) then begin

   gh_bispec_hypot_make_bins,tmp_freq,irf,bin1,bin2,freq_out

   nb = n_elements(freq_out)

   breal  = dblarr(nb) + nan
   bimag  = dblarr(nb) + nan
   bmod   = dblarr(nb) + nan
   bphase = dblarr(nb) + nan
   bicoh  = dblarr(nb) + nan
   nprod  = lonarr(nb)
   npix   = lonarr(nb)

   for ib=0L,nb-1L do begin

      fminb = tmp_freq[bin1[ib]]
      fmaxb = tmp_freq[bin2[ib]]

      w = where(pix_freq ge fminb and pix_freq le fmaxb,nw)

      if(nw gt 0) then begin

         bsr_sum = total(pix_bsr[w],/double)
         bsi_sum = total(pix_bsi[w],/double)
         d1_sum  = total(pix_den1[w],/double)
         d2_sum  = total(pix_den2[w],/double)
         np_sum  = long(total(pix_nprod[w],/double))

         nprod[ib] = np_sum
         npix[ib]  = nw

         if(np_sum gt 0L) then begin

            bre = bsr_sum / double(np_sum)
            bim = bsi_sum / double(np_sum)

            breal[ib]  = bre
            bimag[ib]  = bim
            bmod[ib]   = sqrt(bre^2 + bim^2)
            bphase[ib] = atan(bim,bre)

            if((d1_sum gt 0.0d0) and (d2_sum gt 0.0d0)) then begin
               bicoh[ib] = float((bsr_sum^2 + bsi_sum^2) / (d1_sum*d2_sum))
            endif

         endif

      endif

   endfor

endif else begin

   freq_out = tmp_freq
   breal    = tmp_breal
   bimag    = tmp_bimag
   bmod     = tmp_bmod
   bphase   = tmp_bpha
   bicoh    = tmp_bicoh
   nprod    = tmp_nprod
   npix     = tmp_npix

endelse

;--------------------------------------------------------------------------
; Whole-hypotenuse average from all selected raw sums.
;--------------------------------------------------------------------------
bsr_tot = total(pix_bsr,/double)
bsi_tot = total(pix_bsi,/double)
d1_tot  = total(pix_den1,/double)
d2_tot  = total(pix_den2,/double)

nprod_avg = long(total(pix_nprod,/double))
npix_avg  = nsel

breal_avg  = nan
bimag_avg  = nan
bmod_avg   = nan
bphase_avg = nan
bicoh_avg  = nan

if(nprod_avg gt 0L) then begin

   bre = bsr_tot / double(nprod_avg)
   bim = bsi_tot / double(nprod_avg)

   breal_avg  = bre
   bimag_avg  = bim
   bmod_avg   = sqrt(bre^2 + bim^2)
   bphase_avg = atan(bim,bre)

   if((d1_tot gt 0.0d0) and (d2_tot gt 0.0d0)) then begin
      bicoh_avg = float((bsr_tot^2 + bsi_tot^2) / (d1_tot*d2_tot))
   endif

endif

print,'Extracted hypotenuse'
print,'  f0: ',f0
print,'  width: ',width
print,'  xtype: ',xt
print,'  sym: ',sy
print,'  selected 2D pixels: ',nsel
print,'  output bins: ',n_elements(freq_out)
print,'  total products: ',nprod_avg

end
