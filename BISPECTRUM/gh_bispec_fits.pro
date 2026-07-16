pro gh_bispec_fits,filename,frequency1,frequency2, $
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
;      GH_BISPEC_FITS
;
; PURPOSE:
;      Writes 2D bispectrum products to a multi-extension FITS image file.
;
; EXPLANATION:
;      This routine writes the output of GH_BISPEC2D to a FITS file with
;      one image extension per quantity:
;
;          BREAL
;          BIMAG
;          BMOD
;          BPHASE
;          BICOH
;          NPROD_USED
;          FREQ1
;          FREQ2
;
;      If supplied, it also appends the raw accumulated arrays needed to
;      recompute bicoherence and exactly rebin or estimate errors:
;
;          RAW_BSUM_REAL
;          RAW_BSUM_IMAG
;          RAW_DEN1
;          RAW_DEN2
;
;      If supplied, it can also append the Poisson-denominator-corrected
;      intrinsic observed bicoherence and its auditable corrected denominator
;      sums:
;
;          INT_BICOH
;          RAW_DEN1_CORR
;          RAW_DEN2_CORR
;
;      The first axis is f1 and the second axis is f2, following the
;      array convention used by GH_BISPEC2D:
;
;          image[f1_index,f2_index]
;
; CALLING SEQUENCE:
;      GH_BISPEC_FITS,FILENAME,FREQ1,FREQ2, $
;                     BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD_USED, $
;                     [,RMS=RMS] $
;                     [,/LEAHY] $
;                     [,PRODUCT=PRODUCT] $
;                     [,SRCFILE1=SRCFILE1] $
;                     [,SRCFILE2=SRCFILE2] $
;                     [,SRCFILE3=SRCFILE3] $
;                     [,HISTORY=HISTORY] $
;                     [,COMMAND_HISTORY=COMMAND_HISTORY] $
;                     [,RAW_BSUM_REAL=RAW_BSUM_REAL] $
;                     [,RAW_BSUM_IMAG=RAW_BSUM_IMAG] $
;                     [,RAW_DEN1=RAW_DEN1] $
;                     [,RAW_DEN2=RAW_DEN2] $
;                     [,IRF1=IRF1][,IRF2=IRF2] $
;                     [,/HELP]
;
; INPUTS:
;      FILENAME    = output FITS filename
;      FREQUENCY1  = frequency array for first bispectral axis
;      FREQUENCY2  = frequency array for second bispectral axis
;      BREAL       = real part of bispectrum
;      BIMAG       = imaginary part of bispectrum
;      BMOD        = modulus of bispectrum
;      BPHASE      = biphase, in radians
;      BICOH       = measured squared bicoherence
;      NPROD_USED  = number of products contributing to each pixel
;
; KEYWORDS:
;      RMS         = if set, BREAL/BIMAG/BMOD are labelled as rms^3.
;                    Otherwise they are labelled as native FFT units.
;
;      LEAHY       = if set, BREAL/BIMAG/BMOD are labelled as Leahy^(3/2).
;                    RMS and LEAHY are mutually exclusive.
;
;      RAW_BSUM_REAL = real part of accumulated complex bispectrum sums.
;
;      RAW_BSUM_IMAG = imaginary part of accumulated complex bispectrum sums.
;
;      RAW_DEN1      = accumulated sum |X(f1) X(f2)|^2.
;
;      RAW_DEN2      = accumulated sum |X(f1+f2)|^2.
;
;      INT_BICOH     = Poisson-denominator-corrected intrinsic observed
;                      squared bicoherence. The bispectrum numerator is
;                      unchanged from BICOH.
;
;      RAW_DEN1_CORR = corrected accumulated denominator corresponding to
;                      RAW_DEN1.
;
;      RAW_DEN2_CORR = corrected accumulated denominator corresponding to
;                      RAW_DEN2.
;
;      PRODUCT       = optional product label, e.g. 'BISPEC' or
;                      'CROSS_BISPEC'.
;
;      SRCFILE1      = optional first input data filename.
;
;      SRCFILE2      = optional second input data filename.
;
;      SRCFILE3      = optional third input data filename.
;
;      HISTORY       = optional string or string array written as FITS
;                      HISTORY records in the primary header.
;
;      COMMAND_HISTORY = optional command string or string array written as
;                        FITS HISTORY records. Long strings are split across
;                        multiple HISTORY cards.
;
;      IRF1        = rebin factor used on f1 axis.
;                    Positive = linear rebinning.
;                    Negative = logarithmic rebinning.
;
;      IRF2        = rebin factor used on f2 axis.
;                    Positive = linear rebinning.
;                    Negative = logarithmic rebinning.
;
;      HELP        = print short usage message and return.
;
; EXAMPLE:
;      IDL> gh_bispec2d,'file.fft',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $
;             f1range=[0.1,20.0],f2range=[0.1,20.0],rms=0.001
;
;      IDL> gh_bispec_fits,'out_bispec2d.fits',f1,f2, $
;             breal,bimag,bmod,bphase,bicoh,nprod, $
;             raw_bsum_real=bsr,raw_bsum_imag=bsi, $
;             raw_den1=den1,raw_den2=den2, $
;             rms=0.001,irf1=-100,irf2=-100
;
; COMMON BLOCKS:
;      None
;
; ROUTINES USED:
;      FXHMAKE
;      FXADDPAR
;      MWRFITS
;
; NOTES:
;      This routine assumes that the IDL Astronomy User's Library FITS
;      routines are available.
;
;      Invalid pixels should already be NaN in BREAL/BIMAG/BMOD/BPHASE/BICOH.
;      NPROD_USED should be zero for invalid pixels.
;
; MODIFICATION HISTORY:
;      M. Mendez  21 May 2026  first version, developed from the GHATS
;                              cross-spectrum routines, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC_FITS'
   print,'Write 2D bispectrum products to a multi-extension FITS image.'
   print,' '
   print,'Usage:'
   print,"  gh_bispec_fits,'out.fits',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod"
   print,"  gh_bispec_fits,'out.fits',f1,f2,breal,bimag,bmod,bphase,bicoh,nprod, $"
	   print,'                  raw_bsum_real=bsr,raw_bsum_imag=bsi, $'
	   print,'                  raw_den1=den1,raw_den2=den2, $'
	   print,'                  int_bicoh=bint,den1_corr=d1c,den2_corr=d2c, $'
	   print,"                  product='BISPEC',srcfile1='file.fft', $"
	   print,'                  command_history=cmd, $'
   print,"                  rms=0.001,irf1=-100,irf2=-100"
   print,' '
   return
endif

;--------------------------------------------------------------------------
; Basic checks
;--------------------------------------------------------------------------
if(n_params() lt 9) then begin
   massage,'Usage: gh_bispec_fits,filename,f1,f2,breal,bimag,bmod,bphase,bicoh,nprod'
   retall
endif

outfile = strtrim(string(filename),2)

nf1 = n_elements(frequency1)
nf2 = n_elements(frequency2)

if(nf1 le 0 or nf2 le 0) then begin
   massage,'Frequency arrays are undefined or empty'
   retall
endif

s = size(breal)
if(s[0] ne 2) then begin
   massage,'BREAL must be a 2D array'
   retall
endif

if((s[1] ne nf1) or (s[2] ne nf2)) then begin
   massage,'Array dimensions do not match frequency axes'
   retall
endif

;--------------------------------------------------------------------------
; Check that all image arrays have the same dimensions.
;--------------------------------------------------------------------------
arrays_ok = 1

sb = size(bimag)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(bmod)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(bphase)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(bicoh)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0
sb = size(nprod_used)
if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then arrays_ok = 0

if(arrays_ok eq 0) then begin
   massage,'One or more image arrays do not match frequency axes'
   retall
endif

	have_raw = (n_elements(raw_bsum_real) gt 0) or $
	           (n_elements(raw_bsum_imag) gt 0) or $
	           (n_elements(raw_den1) gt 0) or $
	           (n_elements(raw_den2) gt 0)

	have_intobs = (n_elements(int_bicoh) gt 0) or $
	              (n_elements(den1_corr_kw) gt 0) or $
	              (n_elements(den2_corr_kw) gt 0)

if(have_raw) then begin
   if((n_elements(raw_bsum_real) eq 0) or $
      (n_elements(raw_bsum_imag) eq 0) or $
      (n_elements(raw_den1) eq 0) or $
      (n_elements(raw_den2) eq 0)) then begin
      massage,'Give all raw arrays: RAW_BSUM_REAL, RAW_BSUM_IMAG, RAW_DEN1, RAW_DEN2'
      retall
	endif

if(have_intobs) then begin
   if((n_elements(int_bicoh) eq 0) or $
      (n_elements(den1_corr_kw) eq 0) or $
      (n_elements(den2_corr_kw) eq 0)) then begin
      massage,'Give all intrinsic bicoherence arrays: INT_BICOH, RAW_DEN1_CORR, RAW_DEN2_CORR'
      retall
   endif

   int_ok = 1
   sb = size(int_bicoh)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then int_ok = 0
   sb = size(den1_corr_kw)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then int_ok = 0
   sb = size(den2_corr_kw)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then int_ok = 0

   if(int_ok eq 0) then begin
      massage,'One or more intrinsic-observed arrays do not match frequency axes'
      retall
   endif
endif

   raw_ok = 1
   sb = size(raw_bsum_real)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then raw_ok = 0
   sb = size(raw_bsum_imag)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then raw_ok = 0
   sb = size(raw_den1)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then raw_ok = 0
   sb = size(raw_den2)
   if((sb[0] ne 2) or (sb[1] ne nf1) or (sb[2] ne nf2)) then raw_ok = 0

   if(raw_ok eq 0) then begin
      massage,'One or more raw arrays do not match frequency axes'
      retall
   endif
endif

	if(have_raw) then begin
	   print,'Raw accumulated bispectrum sums will be written.'
	endif else begin
	   print,'Warning: raw accumulated bispectrum sums were not supplied; writing derived products only.'
	endelse

if(have_intobs) then begin
   print,'Poisson-denominator-corrected intrinsic observed bicoherence will be written.'
endif

;--------------------------------------------------------------------------
; Determine frequency WCS-like keywords.
; This assumes native, regularly spaced Fourier frequencies.
;--------------------------------------------------------------------------
df1 = 0.0d0
df2 = 0.0d0

if(nf1 gt 1) then df1 = double(frequency1[1]-frequency1[0])
if(nf2 gt 1) then df2 = double(frequency2[1]-frequency2[0])

if(keyword_set(rms) and keyword_set(leahy)) then begin
   massage,'Use either RMS or LEAHY, not both'
   retall
endif

if(keyword_set(rms)) then begin
   bunit_bispec = 'rms^3'
endif else begin
   if(keyword_set(leahy)) then begin
      bunit_bispec = 'Leahy^(3/2)'
   endif else begin
      bunit_bispec = 'native'
   endelse
endelse

if(n_elements(product) gt 0) then begin
   product_label = strtrim(string(product),2)
endif else begin
   product_label = 'BISPEC'
endelse

;--------------------------------------------------------------------------
; Rebinning metadata for plotting.
;--------------------------------------------------------------------------
if(keyword_set(irf1)) then begin
   rebin1 = irf1
endif else begin
   rebin1 = 1
endelse

if(keyword_set(irf2)) then begin
   rebin2 = irf2
endif else begin
   rebin2 = 1
endelse

axlog1 = byte(rebin1 lt 0)
axlog2 = byte(rebin2 lt 0)

if(n_elements(poisson_method) gt 0) then begin
   poi_method_label = strtrim(string(poisson_method),2)
endif else begin
   poi_method_label = 'UNKNOWN'
endelse

if(n_elements(poisson_level) gt 0) then begin
   poi_level_value = double(poisson_level)
endif else begin
   poi_level_value = !VALUES.D_NAN
endelse

;--------------------------------------------------------------------------
; Primary HDU.
; MWRFITS does not accept a zero-length image in some installations,
; so write one dummy byte in the primary HDU.
;--------------------------------------------------------------------------
primary = bytarr(1)
primary[0] = 0B

mkhdr,hdr,primary,/extend

fxaddpar,hdr,'CREATOR','GH_BISPEC_FITS','Created by GHATS bispectrum writer'
fxaddpar,hdr,'CONTENT','2D bispectrum products'
fxaddpar,hdr,'PRODUCT',product_label,'Stored Fourier product'
if(n_elements(srcfile1) gt 0) then fxaddpar,hdr,'SRCFILE1',strtrim(string(srcfile1),2),'Input file 1'
if(n_elements(srcfile2) gt 0) then fxaddpar,hdr,'SRCFILE2',strtrim(string(srcfile2),2),'Input file 2'
if(n_elements(srcfile3) gt 0) then fxaddpar,hdr,'SRCFILE3',strtrim(string(srcfile3),2),'Input file 3'
fxaddpar,hdr,'NFREQ1',nf1,'Number of f1 bins'
fxaddpar,hdr,'NFREQ2',nf2,'Number of f2 bins'
fxaddpar,hdr,'F1MIN',double(frequency1[0]),'Minimum f1 frequency'
fxaddpar,hdr,'F1MAX',double(frequency1[nf1-1]),'Maximum f1 frequency'
fxaddpar,hdr,'F2MIN',double(frequency2[0]),'Minimum f2 frequency'
fxaddpar,hdr,'F2MAX',double(frequency2[nf2-1]),'Maximum f2 frequency'
fxaddpar,hdr,'DF1',df1,'Frequency spacing on f1 axis'
fxaddpar,hdr,'DF2',df2,'Frequency spacing on f2 axis'
fxaddpar,hdr,'BUNITB',bunit_bispec,'Units of BREAL/BIMAG/BMOD'
fxaddpar,hdr,'RAWSUMS',long(have_raw),'1 if raw accumulated sums are present'
fxaddpar,hdr,'INTBICOH',long(have_intobs),'1 if INT_BICOH is present'
if(have_intobs) then begin
   fxaddpar,hdr,'POIMETH',poi_method_label,'Poisson method for INT_BICOH'
   fxaddpar,hdr,'POILEVEL',poi_level_value,'Scalar/mean Poisson level used'
   if(n_elements(poisson_freq_range) eq 2) then begin
      fxaddpar,hdr,'POIFMIN',double(poisson_freq_range[0]),'Poisson range fmin'
      fxaddpar,hdr,'POIFMAX',double(poisson_freq_range[1]),'Poisson range fmax'
   endif
   fxaddpar,hdr,'BICNUM','UNCHANGED','INT_BICOH numerator correction'
   fxaddpar,hdr,'BICDEN','POISSON_CORR','INT_BICOH denominator correction'
endif
if(n_elements(corr_mode) gt 0) then fxaddpar,hdr,'CORRMODE',strtrim(string(corr_mode),2),'Correction mode'
if(n_elements(diag_corr) gt 0) then fxaddpar,hdr,'DIAGCORR',strtrim(string(diag_corr),2),'Diagonal correction'
if(n_elements(dtcovar) gt 0) then fxaddpar,hdr,'DTCOVAR',strtrim(string(dtcovar),2),'Dead-time covariance status'
if(n_elements(numbias) gt 0) then fxaddpar,hdr,'NUMBIAS',strtrim(string(numbias),2),'Numerator bias correction'
if(n_elements(band1) gt 0) then fxaddpar,hdr,'BAND1',strtrim(string(band1),2),'Band in X1(f1)'
if(n_elements(band2) gt 0) then fxaddpar,hdr,'BAND2',strtrim(string(band2),2),'Band in X2(f2)'
if(n_elements(band3) gt 0) then fxaddpar,hdr,'BAND3',strtrim(string(band3),2),'Band in conjugated X3(F)'
if(n_elements(bandovlp) gt 0) then fxaddpar,hdr,'BANDOVLP',strtrim(string(bandovlp),2),'Whether input bands overlap'
if(n_elements(xorder) gt 0) then fxaddpar,hdr,'XORDER',strtrim(string(xorder),2),'Fourier product ordering'
if(n_elements(fconv) gt 0) then fxaddpar,hdr,'FCONV',strtrim(string(fconv),2),'Fourier convention note'
nhist = 0L
if(n_elements(history) gt 0) then nhist = nhist + n_elements(history)
if(n_elements(command_history) gt 0) then nhist = nhist + n_elements(command_history)
if(nhist gt 0L) then begin
   hist_all = strarr(nhist)
   ihist = 0L
   if(n_elements(history) gt 0) then begin
      for ih=0L,n_elements(history)-1L do begin
         hist_all[ihist] = strtrim(string(history[ih]),2)
         ihist = ihist + 1L
      endfor
   endif
   if(n_elements(command_history) gt 0) then begin
      for ih=0L,n_elements(command_history)-1L do begin
         hist_all[ihist] = strtrim(string(command_history[ih]),2)
         ihist = ihist + 1L
      endfor
   endif
   for ih=0L,nhist-1L do begin
      htxt = strtrim(hist_all[ih],2)
      if(strlen(htxt) eq 0) then htxt = ' '
      p0 = 0L
      while p0 lt strlen(htxt) do begin
         fxaddpar,hdr,'HISTORY',strmid(htxt,p0,68)
         p0 = p0 + 68L
      endwhile
   endfor
endif
fxaddpar,hdr,'NEXTEND',8L + 4L*long(have_raw) + 3L*long(have_intobs),'Number of image extensions'
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,primary,outfile,hdr,/create

;--------------------------------------------------------------------------
; BREAL
;--------------------------------------------------------------------------
mkhdr,hdr,breal,/image
fxaddpar,hdr,'EXTNAME','BREAL'
fxaddpar,hdr,'BUNIT',bunit_bispec
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,breal,outfile,hdr

;--------------------------------------------------------------------------
; BIMAG
;--------------------------------------------------------------------------
mkhdr,hdr,bimag,/image
fxaddpar,hdr,'EXTNAME','BIMAG'
fxaddpar,hdr,'BUNIT',bunit_bispec
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,bimag,outfile,hdr

;--------------------------------------------------------------------------
; BMOD
;--------------------------------------------------------------------------
mkhdr,hdr,bmod,/image
fxaddpar,hdr,'EXTNAME','BMOD'
fxaddpar,hdr,'BUNIT',bunit_bispec
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,bmod,outfile,hdr

;--------------------------------------------------------------------------
; BPHASE
;--------------------------------------------------------------------------
mkhdr,hdr,bphase,/image
fxaddpar,hdr,'EXTNAME','BPHASE'
fxaddpar,hdr,'BUNIT','rad'
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,bphase,outfile,hdr

;--------------------------------------------------------------------------
; BICOH
;--------------------------------------------------------------------------
mkhdr,hdr,bicoh,/image
fxaddpar,hdr,'EXTNAME','BICOH'
fxaddpar,hdr,'BUNIT','dimensionless'
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,bicoh,outfile,hdr

;--------------------------------------------------------------------------
; NPROD_USED
; Use a count image. Invalid pixels should have value 0.
;--------------------------------------------------------------------------
nprod_long = long(nprod_used)

mkhdr,hdr,nprod_long,/image
fxaddpar,hdr,'EXTNAME','NPROD_USED'
fxaddpar,hdr,'BUNIT','count'
fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
fxaddpar,hdr,'CUNIT1','Hz'
fxaddpar,hdr,'CRPIX1',1.0d0
fxaddpar,hdr,'CRVAL1',double(frequency1[0])
fxaddpar,hdr,'CDELT1',df1
fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
fxaddpar,hdr,'CUNIT2','Hz'
fxaddpar,hdr,'CRPIX2',1.0d0
fxaddpar,hdr,'CRVAL2',double(frequency2[0])
fxaddpar,hdr,'CDELT2',df2
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,nprod_long,outfile,hdr

;--------------------------------------------------------------------------
; Frequency axes.
; These are needed because rebinned/logarithmic axes are not fully described
; by CRVAL/CDELT.
;--------------------------------------------------------------------------
freq1_double = double(frequency1)
mkhdr,hdr,freq1_double,/image
fxaddpar,hdr,'EXTNAME','FREQ1'
fxaddpar,hdr,'BUNIT','Hz'
fxaddpar,hdr,'REBIN1',long(rebin1),'Rebin factor on f1 axis'
fxaddpar,hdr,'AXLOG1',long(axlog1),'1 if f1 axis is logarithmically rebinned'
mwrfits,freq1_double,outfile,hdr

freq2_double = double(frequency2)
mkhdr,hdr,freq2_double,/image
fxaddpar,hdr,'EXTNAME','FREQ2'
fxaddpar,hdr,'BUNIT','Hz'
fxaddpar,hdr,'REBIN2',long(rebin2),'Rebin factor on f2 axis'
fxaddpar,hdr,'AXLOG2',long(axlog2),'1 if f2 axis is logarithmically rebinned'
mwrfits,freq2_double,outfile,hdr

if(have_raw) then begin

   ;-----------------------------------------------------------------------
   ; RAW_BSUM_REAL
   ;-----------------------------------------------------------------------
   mkhdr,hdr,raw_bsum_real,/image
   fxaddpar,hdr,'EXTNAME','RAW_BSUM_REAL'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,raw_bsum_real,outfile,hdr

   ;-----------------------------------------------------------------------
   ; RAW_BSUM_IMAG
   ;-----------------------------------------------------------------------
   mkhdr,hdr,raw_bsum_imag,/image
   fxaddpar,hdr,'EXTNAME','RAW_BSUM_IMAG'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,raw_bsum_imag,outfile,hdr

   ;-----------------------------------------------------------------------
   ; RAW_DEN1
   ;-----------------------------------------------------------------------
   mkhdr,hdr,raw_den1,/image
   fxaddpar,hdr,'EXTNAME','RAW_DEN1'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,raw_den1,outfile,hdr

   ;-----------------------------------------------------------------------
   ; RAW_DEN2
   ;-----------------------------------------------------------------------
   mkhdr,hdr,raw_den2,/image
   fxaddpar,hdr,'EXTNAME','RAW_DEN2'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,raw_den2,outfile,hdr

	endif

if(have_intobs) then begin

   ;-----------------------------------------------------------------------
   ; INT_BICOH
   ;-----------------------------------------------------------------------
   mkhdr,hdr,int_bicoh,/image
   fxaddpar,hdr,'EXTNAME','INT_BICOH'
   fxaddpar,hdr,'BUNIT','dimensionless'
   fxaddpar,hdr,'POIMETH',poi_method_label,'Poisson method'
   fxaddpar,hdr,'POILEVEL',poi_level_value,'Scalar/mean Poisson level used'
   fxaddpar,hdr,'BICNUM','UNCHANGED','Bispectrum numerator correction'
   fxaddpar,hdr,'BICDEN','POISSON_CORR','Bicoherence denominator correction'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,int_bicoh,outfile,hdr

   ;-----------------------------------------------------------------------
   ; RAW_DEN1_CORR
   ;-----------------------------------------------------------------------
   mkhdr,hdr,den1_corr_kw,/image
   fxaddpar,hdr,'EXTNAME','RAW_DEN1_CORR'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'POIMETH',poi_method_label,'Poisson method'
   fxaddpar,hdr,'BICNUM','UNCHANGED','Bispectrum numerator correction'
   fxaddpar,hdr,'BICDEN','POISSON_CORR','Bicoherence denominator correction'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,den1_corr_kw,outfile,hdr

   ;-----------------------------------------------------------------------
   ; RAW_DEN2_CORR
   ;-----------------------------------------------------------------------
   mkhdr,hdr,den2_corr_kw,/image
   fxaddpar,hdr,'EXTNAME','RAW_DEN2_CORR'
   fxaddpar,hdr,'BUNIT','sum'
   fxaddpar,hdr,'POIMETH',poi_method_label,'Poisson method'
   fxaddpar,hdr,'BICNUM','UNCHANGED','Bispectrum numerator correction'
   fxaddpar,hdr,'BICDEN','POISSON_CORR','Bicoherence denominator correction'
   fxaddpar,hdr,'CTYPE1','FREQ1','Bispectrum f1 axis'
   fxaddpar,hdr,'CUNIT1','Hz'
   fxaddpar,hdr,'CRPIX1',1.0d0
   fxaddpar,hdr,'CRVAL1',double(frequency1[0])
   fxaddpar,hdr,'CDELT1',df1
   fxaddpar,hdr,'CTYPE2','FREQ2','Bispectrum f2 axis'
   fxaddpar,hdr,'CUNIT2','Hz'
   fxaddpar,hdr,'CRPIX2',1.0d0
   fxaddpar,hdr,'CRVAL2',double(frequency2[0])
   fxaddpar,hdr,'CDELT2',df2
   mwrfits,den2_corr_kw,outfile,hdr

endif

	print,'Wrote ',outfile

end
