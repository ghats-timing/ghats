pro gh_wcb_append_image,filename,arr,extname,frequency1,frequency2, $
                         corr_mode=corr_mode,diag_corr=diag_corr, $
                         dtcovar=dtcovar,numbias=numbias

df1 = !values.d_nan
df2 = !values.d_nan
if(n_elements(frequency1) gt 1) then df1 = double(frequency1[1])-double(frequency1[0])
if(n_elements(frequency2) gt 1) then df2 = double(frequency2[1])-double(frequency2[0])

mkhdr,hdr,arr,/image
fxaddpar,hdr,'EXTNAME',strtrim(string(extname),2)
if(n_elements(corr_mode) gt 0) then fxaddpar,hdr,'CORRMODE',strtrim(string(corr_mode),2)
if(n_elements(diag_corr) gt 0) then fxaddpar,hdr,'DIAGCORR',strtrim(string(diag_corr),2)
if(n_elements(dtcovar) gt 0) then fxaddpar,hdr,'DTCOVAR',strtrim(string(dtcovar),2)
if(n_elements(numbias) gt 0) then fxaddpar,hdr,'NUMBIAS',strtrim(string(numbias),2)
fxaddpar,hdr,'CTYPE1','FREQ1'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
mwrfits,arr,filename,hdr

end


pro gh_write_cross_bispec_fits,filename,frequency1,frequency2, $
                               breal,bimag,bmod,bphase,bicoh,nprod_used, $
                               raw_bsum_real=raw_bsum_real, $
                               raw_bsum_imag=raw_bsum_imag, $
                               raw_den1=raw_den1,raw_den2=raw_den2, $
                               int_bicoh_indep=int_bicoh_indep, $
                               den1_corr_indep=den1_corr_indep, $
                               den2_corr_indep=den2_corr_indep, $
                               intr_valid=intr_valid, $
                               corr_flags=corr_flags, $
                               diag_flag=diag_flag, $
                               poisson_method=poisson_method, $
                               poisson_level=poisson_level, $
                               poisson_freq_range=poisson_freq_range, $
                               product=product, $
                               srcfile1=srcfile1,srcfile2=srcfile2,srcfile3=srcfile3, $
                               band1=band1,band2=band2,band3=band3, $
                               bandovlp=bandovlp,history=history, $
                               leahy=leahy,rms=rms,irf1=irf1,irf2=irf2,help=help
;+
; NAME:
;      GH_WRITE_CROSS_BISPEC_FITS
;
; PURPOSE:
;      Write matched cross-bispectrum products with explicit independent-noise
;      corrected bicoherence extensions.
;
; EXPLANATION:
;      This wrapper preserves the standard GH_WRITE_BISPEC_FITS raw-product
;      schema, but does not write an unqualified INT_BICOH extension for cross
;      products.  Instead it appends:
;
;          INT_BICOH_INDEP
;          RAW_DEN1_CORR_INDEP
;          RAW_DEN2_CORR_INDEP
;          INTR_VALID
;          CORR_FLAGS
;          DIAG_FLAG
;
;      Metadata explicitly states that dead-time cross-band covariance is not
;      included and that the complex numerator is not corrected.
;-

if(keyword_set(help)) then begin
   print,'GH_WRITE_CROSS_BISPEC_FITS'
   print,'Write cross-bispectrum FITS with explicit INDEP_NOISE extensions.'
   return
endif

hist = history
if(n_elements(hist) eq 0) then hist = ['Matched cross-bispectrum product']
hist = [hist, $
        'INT_BICOH_INDEP assumes independent cross-band Fourier noise.', $
        'Dead-time cross-band covariance is not included.', $
        'Cross-spectrum high-frequency Re subtraction is not used as a denominator correction.']

gh_write_bispec_fits,filename,frequency1,frequency2, $
                     breal,bimag,bmod,bphase,bicoh,nprod_used, $
                     raw_bsum_real=raw_bsum_real,raw_bsum_imag=raw_bsum_imag, $
                     raw_den1=raw_den1,raw_den2=raw_den2, $
                     poisson_method=poisson_method,poisson_level=poisson_level, $
                     poisson_freq_range=poisson_freq_range, $
                     corr_mode='INDEP_NOISE',diag_corr='FOURTH_MOMENT', $
                     dtcovar='NOT_INCLUDED',numbias='NONE_ASSUMED', $
                     band1=band1,band2=band2,band3=band3,bandovlp=bandovlp, $
                     xorder='X1*X2*CONJ(X3)', $
                     fconv='GHATS positive lag: hard lags soft', $
                     product=product,srcfile1=srcfile1,srcfile2=srcfile2, $
                     srcfile3=srcfile3,history=hist,leahy=leahy,rms=rms, $
                     irf1=irf1,irf2=irf2

if(n_elements(int_bicoh_indep) gt 0) then $
   gh_wcb_append_image,filename,int_bicoh_indep,'INT_BICOH_INDEP', $
                       frequency1,frequency2,corr_mode='INDEP_NOISE', $
                       diag_corr='FOURTH_MOMENT',dtcovar='NOT_INCLUDED', $
                       numbias='NONE_ASSUMED'
if(n_elements(den1_corr_indep) gt 0) then $
   gh_wcb_append_image,filename,den1_corr_indep,'RAW_DEN1_CORR_INDEP', $
                       frequency1,frequency2,corr_mode='INDEP_NOISE', $
                       diag_corr='FOURTH_MOMENT',dtcovar='NOT_INCLUDED', $
                       numbias='NONE_ASSUMED'
if(n_elements(den2_corr_indep) gt 0) then $
   gh_wcb_append_image,filename,den2_corr_indep,'RAW_DEN2_CORR_INDEP', $
                       frequency1,frequency2,corr_mode='INDEP_NOISE', $
                       diag_corr='FOURTH_MOMENT',dtcovar='NOT_INCLUDED', $
                       numbias='NONE_ASSUMED'
if(n_elements(intr_valid) gt 0) then $
   gh_wcb_append_image,filename,byte(intr_valid),'INTR_VALID',frequency1,frequency2, $
                       corr_mode='INDEP_NOISE',diag_corr='FOURTH_MOMENT', $
                       dtcovar='NOT_INCLUDED',numbias='NONE_ASSUMED'
if(n_elements(corr_flags) gt 0) then $
   gh_wcb_append_image,filename,long(corr_flags),'CORR_FLAGS',frequency1,frequency2, $
                       corr_mode='INDEP_NOISE',diag_corr='FOURTH_MOMENT', $
                       dtcovar='NOT_INCLUDED',numbias='NONE_ASSUMED'
if(n_elements(diag_flag) gt 0) then $
   gh_wcb_append_image,filename,byte(diag_flag),'DIAG_FLAG',frequency1,frequency2, $
                       corr_mode='INDEP_NOISE',diag_corr='FOURTH_MOMENT', $
                       dtcovar='NOT_INCLUDED',numbias='NONE_ASSUMED'

end
