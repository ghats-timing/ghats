pro gh_cross_bispec2d,filename1,filename2,filename3,frequency1,frequency2, $
                       breal,bimag,bmod,bphase,bicoh,nprod_used,nseg, $
                       f1range=f1range,f2range=f2range, $
                       index=index,time=time,sel=sel,rate=rate, $
                       rms1=back1,rms2=back2,rms3=back3,leahy=leahy, $
                       raw_bsum_real=raw_bsum_real, $
                       raw_bsum_imag=raw_bsum_imag, $
                       raw_den1=raw_den1, $
                       raw_den2=raw_den2, $
                       help=help
;+
; NAME:
;      GH_CROSS_BISPEC2D
;
; PURPOSE:
;      Computes the 2D cross-bispectrum and squared cross-bicoherence
;      from three compatible GHATS FFT files.
;
; EXPLANATION:
;      This routine reads three GHATS FFT files and computes
;
;          B_123(f1,f2) = < X1(f1) X2(f2) X3*(f1+f2) >
;
;      over a selected rectangular region in the positive-frequency
;      (f1,f2) plane.
;
;      This V1 is strict: the three FFT files must have the same time grid,
;      same number of FFTs, same FFT length, same start time, and matching
;      RMJD values when read in lockstep. If not, the routine stops.
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
;      Pixels for which f1+f2 is above the Nyquist frequency are invalid.
;
;      Raw accumulated sums can be returned:
;
;          RAW_BSUM_REAL = real part of sum X1 X2 X3*
;          RAW_BSUM_IMAG = imaginary part of sum X1 X2 X3*
;          RAW_DEN1      = sum |X1 X2|^2
;          RAW_DEN2      = sum |X3|^2
;
; CALLING SEQUENCE:
;      GH_CROSS_BISPEC2D,FILE1,FILE2,FILE3,FREQ1,FREQ2, $
;                         BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD,NSEG, $
;                         F1RANGE=F1RANGE,F2RANGE=F2RANGE, [/LEAHY]
;
; INPUTS:
;      FILENAME1 = file providing X1(f1)
;      FILENAME2 = file providing X2(f2)
;      FILENAME3 = file providing X3(f1+f2)
;
; OUTPUTS:
;      FREQUENCY1  = frequency array for first cross-bispectral axis
;      FREQUENCY2  = frequency array for second cross-bispectral axis
;      BREAL       = real part of averaged cross-bispectrum
;      BIMAG       = imaginary part of averaged cross-bispectrum
;      BMOD        = modulus of averaged cross-bispectrum
;      BPHASE      = cross-biphase in radians
;      BICOH       = measured squared cross-bicoherence
;      NPROD_USED  = number of products contributing to each pixel
;      NSEG        = number of selected FFT segments
;
; KEYWORDS:
;      F1RANGE, F2RANGE = mandatory two-element frequency ranges.
;      INDEX, TIME, SEL, RATE = selection keywords, as in GH_BISPEC2D.
;      LEAHY = convert each FFT amplitude using sqrt(2/cnts_i).
;      RMS1,RMS2,RMS3 = optional rms^3-like conversion of BREAL/BIMAG/BMOD.
;                       Use only if you know the background rates for all bands.
;                       RMS keywords and LEAHY are mutually exclusive.
;      RAW_* = optional raw accumulated sums.
;
; NOTES:
;      Cross-bicoherence denominator used here is
;
;          b^2 = |sum X1 X2 X3*|^2 / [sum |X1 X2|^2 sum |X3|^2]
;
;      For exact rebinning, rebin RAW_BSUM_REAL, RAW_BSUM_IMAG, RAW_DEN1,
;      RAW_DEN2 and NPROD_USED, then recompute derived products.
;
; MODIFICATION HISTORY:
;      M. Mendez  25 May 2026  first V1, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_CROSS_BISPEC2D'
   print,'Compute 2D cross-bispectrum B_123(f1,f2)=<X1(f1) X2(f2) X3*(f1+f2)>.'
   print,' '
   print,'Usage:'
   print,"  gh_cross_bispec2d,file1,file2,file3,f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $"
   print,'                      f1range=[0.1,20.0],f2range=[0.1,20.0],/leahy'
   print,' '
   print,'Optional raw sums:'
   print,'  raw_bsum_real=bsr, raw_bsum_imag=bsi, raw_den1=den1, raw_den2=den2'
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

if(keyword_set(leahy) and (keyword_set(back1) or keyword_set(back2) or keyword_set(back3))) then begin
   massage,'Use either LEAHY or RMS1/RMS2/RMS3, not both'
   retall
endif

if((keyword_set(back1) or keyword_set(back2) or keyword_set(back3)) and $
   ~(keyword_set(back1) and keyword_set(back2) and keyword_set(back3))) then begin
   massage,'If using RMS conversion, provide RMS1, RMS2 and RMS3'
   retall
endif

;--------------------------------------------------------------------------
; Open FFT files and read headers.
;--------------------------------------------------------------------------
ghats_openfft,filename1,unit1,/dialog
ghats_openfft,filename2,unit2,/dialog
ghats_openfft,filename3,unit3,/dialog

ntrafos1 = 0L
ntrafos2 = 0L
ntrafos3 = 0L

dummy1 = bytarr(100)
dummy2 = bytarr(100)
dummy3 = bytarr(100)

gh_version_string1 = '                '
gh_version_string2 = '                '
gh_version_string3 = '                '
observatory1 = '                '
observatory2 = '                '
observatory3 = '                '
instrument1  = '                '
instrument2  = '                '
instrument3  = '                '
target1      = '                '
target2      = '                '
target3      = '                '
rmjd01 = 0.0D0
rmjd02 = 0.0D0
rmjd03 = 0.0D0

ghats_getheader,unit1,gh_version_string1,observatory1,instrument1,target1,rmjd01, $
                     nft_header1,T1,ntrafos1,e1,proliferation1,baryflag1, $
                     n_spectral_bins1,background_flag1,dummy1
muflag1 = dummy1[0]

ghats_getheader,unit2,gh_version_string2,observatory2,instrument2,target2,rmjd02, $
                     nft_header2,T2,ntrafos2,e2,proliferation2,baryflag2, $
                     n_spectral_bins2,background_flag2,dummy2
muflag2 = dummy2[0]

ghats_getheader,unit3,gh_version_string3,observatory3,instrument3,target3,rmjd03, $
                     nft_header3,T3,ntrafos3,e3,proliferation3,baryflag3, $
                     n_spectral_bins3,background_flag3,dummy3
muflag3 = dummy3[0]

;--------------------------------------------------------------------------
; Strict compatibility checks.
;--------------------------------------------------------------------------
if(nft_header1 ne nft_header2 or nft_header1 ne nft_header3) then begin
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'FFT lengths are not compatible'
   retall
endif

tol_t = 1.0d-6

if((abs(T1-T2) gt tol_t) or (abs(T1-T3) gt tol_t)) then begin
   print,T1,T2,T3
   massage,'Time resolution not compatible'
   retall
endif

tol_rmjd0 = 1.0d-8

if((abs(rmjd01-rmjd02) gt tol_rmjd0) or (abs(rmjd01-rmjd03) gt tol_rmjd0)) then begin
   print,'WARNING: Start times differ:'
   print,format='(3F20.10)',rmjd01,rmjd02,rmjd03
   print,'Proceeding, but segment RMJD values will be checked line by line.'
endif

if(ntrafos1 ne ntrafos2 or ntrafos1 ne ntrafos3) then begin
   print,ntrafos1,ntrafos2,ntrafos3
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'Number of FFTs not compatible'
   retall
endif

nft = nft_header1/2
T   = T1
df  = T
ntrafos = ntrafos1

allfreq = (findgen(nft)+1.0) * df
fnyq = allfreq[nft-1]

w1 = where(allfreq ge f1range[0] and allfreq le f1range[1], nf1)
w2 = where(allfreq ge f2range[0] and allfreq le f2range[1], nf2)

if(nf1 le 0) then begin
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'No Fourier bins found in F1RANGE'
   retall
endif

if(nf2 le 0) then begin
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'No Fourier bins found in F2RANGE'
   retall
endif

frequency1 = allfreq[w1]
frequency2 = allfreq[w2]

;--------------------------------------------------------------------------
; Allocate arrays.
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

wgood = where((goodarray ge firstfft) and (goodarray le lastfft), ngood)
if(ngood le 0) then begin
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'No FFTs found in selected INDEX/SEL range'
   retall
endif

goodarray = goodarray[wgood]
lastindex = max(goodarray)

;--------------------------------------------------------------------------
; Read FFTs and accumulate cross-bispectral sums.
;--------------------------------------------------------------------------
rdata1 = complexarr(nft)
rdata2 = complexarr(nft)
rdata3 = complexarr(nft)

nseg = 0L
flux1 = 0.0d0
flux2 = 0.0d0
flux3 = 0.0d0

time_tol = 1.0d-4

print,'Computing 2D cross-bispectrum'
print,'FILE1: ',filename1
print,'FILE2: ',filename2
print,'FILE3: ',filename3
print,'F1RANGE: ',frequency1[0],' - ',frequency1[nf1-1],' Hz (',nf1,' bins)'
print,'F2RANGE: ',frequency2[0],' - ',frequency2[nf2-1],' Hz (',nf2,' bins)'

for itrafos=0L,ntrafos-1L do begin

   read_fft_line,unit1,muflag1,rmjd1,cnts1,poisson1,current_vle_rate1,fndet1,rdata1
   read_fft_line,unit2,muflag2,rmjd2,cnts2,poisson2,current_vle_rate2,fndet2,rdata2
   read_fft_line,unit3,muflag3,rmjd3,cnts3,poisson3,current_vle_rate3,fndet3,rdata3

   if((abs(rmjd1-rmjd2) gt time_tol) or (abs(rmjd1-rmjd3) gt time_tol)) then begin
      print,rmjd1,rmjd2,rmjd3
      free_lun,unit1 & free_lun,unit2 & free_lun,unit3
      massage,'Times of FFTs not compatible'
      retall
   endif

   t_1 = ((rmjd1-rmjd01)*86400.0d0)
   t_2 = t_1 + 1.0d0/df

   ; For cross products, use total rate for RATE selection.
   cratem = (cnts1 + cnts2 + cnts3) * df

   gotcha = where(goodarray eq itrafos)

   if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (cratem ge firstrate) and (cratem le lastrate) and $
      (gotcha[0] ge 0)) then begin

      if(keyword_set(leahy)) then begin
         if((cnts1 le 0.0) or (cnts2 le 0.0) or (cnts3 le 0.0)) then begin
            free_lun,unit1 & free_lun,unit2 & free_lun,unit3
            massage,'Cannot apply LEAHY normalization: segment has cnts <= 0'
            retall
         endif
         ampnorm1 = sqrt(2.0d0/double(cnts1))
         ampnorm2 = sqrt(2.0d0/double(cnts2))
         ampnorm3 = sqrt(2.0d0/double(cnts3))
      endif else begin
         ampnorm1 = 1.0d0
         ampnorm2 = 1.0d0
         ampnorm3 = 1.0d0
      endelse

      nseg = nseg + 1L
      flux1 = flux1 + cnts1 * df
      flux2 = flux2 + cnts2 * df
      flux3 = flux3 + cnts3 * df

      for ii1=0L,nf1-1L do begin

         i1 = w1[ii1]

         for ii2=0L,nf2-1L do begin

            i2 = w2[ii2]
            isum = i1 + i2 + 1L

            if(isum lt nft) then begin

               x1 = rdata1[i1]    * ampnorm1
               x2 = rdata2[i2]    * ampnorm2
               x3 = rdata3[isum] * ampnorm3

               z = x1 * x2 * conj(x3)

               bsum[ii1,ii2] = bsum[ii1,ii2] + z
               den1[ii1,ii2] = den1[ii1,ii2] + double(abs(x1*x2)^2)
               den2[ii1,ii2] = den2[ii1,ii2] + double(abs(x3)^2)
               nprod_used[ii1,ii2] = nprod_used[ii1,ii2] + 1L

            endif

         endfor
      endfor

   endif

   if(itrafos gt (lastindex-1)) then goto,finished_reading

endfor

finished_reading:

free_lun,unit1
free_lun,unit2
free_lun,unit3

if(nseg le 0) then begin
   massage,'No FFTs retrieved for input selection'
   retall
endif

flux1 = flux1 / nseg
flux2 = flux2 / nseg
flux3 = flux3 / nseg

print,'  ',strtrim(string(nseg),1),' FFTs selected'

;--------------------------------------------------------------------------
; Average cross-bispectrum and compute derived products.
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

      endif

   endfor
endfor

;--------------------------------------------------------------------------
; Optional rms^3-like conversion.
; Raw sums remain in internal normalization.
;--------------------------------------------------------------------------
if(keyword_set(back1)) then begin

   if(back1 ge flux1) then begin
      massage,'Error: background flux higher than source+bkg flux in band 1'
      retall
   endif
   if(back2 ge flux2) then begin
      massage,'Error: background flux higher than source+bkg flux in band 2'
      retall
   endif
   if(back3 ge flux3) then begin
      massage,'Error: background flux higher than source+bkg flux in band 3'
      retall
   endif

   ; Approximate cross-bispectrum rms^3 normalization.
   norm_bispec = (2.0d0^1.5d0) * (df^1.5d0) / $
                 ((flux1-back1)*(flux2-back2)*(flux3-back3))

   breal = breal * norm_bispec
   bimag = bimag * norm_bispec
   bmod  = bmod  * norm_bispec

   print,'Cross-bispectrum units: rms^3-like'

endif else begin

   if(keyword_set(leahy)) then begin
      print,'Cross-bispectrum units: Leahy^(3/2)'
   endif else begin
      print,'Cross-bispectrum units: native FFT units'
   endelse

endelse

raw_bsum_real = float(bsum)
raw_bsum_imag = imaginary(bsum)
raw_den1      = den1
raw_den2      = den2

end
