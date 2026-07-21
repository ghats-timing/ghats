pro gh_write_bispec_fits,filename,frequency1,frequency2, $
                         breal,bimag,bmod,bphase,bicoh,nprod_used, $
                         raw_bsum_real=raw_bsum_real, $
	                         raw_bsum_imag=raw_bsum_imag, $
	                         raw_den1=raw_den1, $
	                         raw_den2=raw_den2, $
	                         int_bicoh=int_bicoh, $
	                         den1_corr=den1_corr_kw, $
	                         den2_corr=den2_corr_kw, $
	                         poisson_method=poisson_method, $
	                         poisson_level=poisson_level, $
	                         poisson_freq_range=poisson_freq_range, $
	                         corr_mode=corr_mode,diag_corr=diag_corr, $
	                         dtcovar=dtcovar,numbias=numbias, $
	                         band1=band1,band2=band2,band3=band3, $
	                         bandovlp=bandovlp,xorder=xorder,fconv=fconv, $
	                         product=product, $
                         srcfile1=srcfile1,srcfile2=srcfile2,srcfile3=srcfile3, $
                         history=history,command_history=command_history, $
                         rms=rms,leahy=leahy,irf1=irf1,irf2=irf2,help=help
;+
; NAME:
;      GH_WRITE_BISPEC_FITS
;
; PURPOSE:
;      Clear-name wrapper for GH_BISPEC_FITS.
;
; EXPLANATION:
;      GH_BISPEC_FITS is kept for backward compatibility. New code can call
;      GH_WRITE_BISPEC_FITS with the same arguments and optional raw arrays.
;      PRODUCT and SRCFILE1/2/3 are optional FITS provenance keywords.
;-

if(keyword_set(help)) then begin
   print,''
   print,'GH_WRITE_BISPEC_FITS'
   print,''
   print,'Clear-name wrapper for GH_BISPEC_FITS.'
   print,'It accepts the same arguments and keywords; detailed help follows.'
   gh_bispec_fits, /help
   return
endif

gh_bispec_fits,filename,frequency1,frequency2, $
               breal,bimag,bmod,bphase,bicoh,nprod_used, $
               raw_bsum_real=raw_bsum_real, $
	               raw_bsum_imag=raw_bsum_imag, $
	               raw_den1=raw_den1, $
	               raw_den2=raw_den2, $
	               int_bicoh=int_bicoh, $
	               den1_corr=den1_corr_kw, $
	               den2_corr=den2_corr_kw, $
	               poisson_method=poisson_method, $
	               poisson_level=poisson_level, $
	               poisson_freq_range=poisson_freq_range, $
	               corr_mode=corr_mode,diag_corr=diag_corr, $
	               dtcovar=dtcovar,numbias=numbias, $
	               band1=band1,band2=band2,band3=band3, $
	               bandovlp=bandovlp,xorder=xorder,fconv=fconv, $
	               product=product, $
               srcfile1=srcfile1,srcfile2=srcfile2,srcfile3=srcfile3, $
               history=history,command_history=command_history, $
               rms=rms,leahy=leahy,irf1=irf1,irf2=irf2,help=help

end
