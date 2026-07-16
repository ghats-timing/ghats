pro gh_cross_bispec2d_match,filename1,filename2,filename3,frequency1,frequency2, $
                           breal,bimag,bmod,bphase,bicoh,nprod_used,nseg, $
                           f1range=f1range,f2range=f2range, $
                           index=index,time=time,sel=sel,rate=rate, $
                           rms1=back1,rms2=back2,rms3=back3,leahy=leahy, $
                           raw_bsum_real=raw_bsum_real, $
                           raw_bsum_imag=raw_bsum_imag, $
                           raw_den1=raw_den1, $
                           raw_den2=raw_den2, $
                           poisson=poisson_key, $
                           manual_poisson=manual_poisson, $
                           zhang=zhang, $
                           int_bicoh_indep=int_bicoh_indep, $
                           den1_corr_indep=den1_corr_indep, $
                           den2_corr_indep=den2_corr_indep, $
                           intr_valid=intr_valid, $
                           corr_flags=corr_flags, $
                           diag_flag=diag_flag, $
                           poi_method_used=poisson_method_used, $
                           poi_level_used=poisson_level_used, $
                           time_tol=time_tol,freq_tol=freq_tol, $
                           help=help
;+
; NAME:
;      GH_CROSS_BISPEC2D_MATCH
;
; PURPOSE:
;      Computes a 2D cross-bispectrum from three GHATS FFT files, matching
;      FFT segments by RMJD instead of assuming identical segment ordering.
;
; EXPLANATION:
;      This routine computes
;
;          B_123(f1,f2) = < X1(f1) X2(f2) X3*(f1+f2) >
;
;      where X1 is read from FILENAME1, X2 from FILENAME2, and X3 from
;      FILENAME3. The three input FFT files do not need to have the same
;      number of FFT segments or the same RMJD0. Segments are matched by
;      their per-segment RMJD values within TIME_TOL.
;
;      The frequency grids may differ slightly. The first axis uses the
;      frequency grid of FILENAME1 and the second axis uses the grid of
;      FILENAME2. For the sum frequency f1+f2, the nearest Fourier bin in
;      FILENAME3 is used if it lies within FREQ_TOL.
;
;      The GHATS FFT files do not include the DC bin:
;
;          rdata[0] = X(df)
;          rdata[1] = X(2df)
;
; CALLING SEQUENCE:
;      GH_CROSS_BISPEC2D_MATCH,FILE1,FILE2,FILE3,F1,F2, $
;                              BREAL,BIMAG,BMOD,BPHASE,BICOH,NPROD,NSEG, $
;                              F1RANGE=F1RANGE,F2RANGE=F2RANGE, [/LEAHY]
;
; KEYWORDS:
;      F1RANGE/F2RANGE  required frequency ranges.
;      INDEX/TIME/SEL/RATE selection is applied to the matched file1 index,
;                         file1-relative time, and summed count rate.
;      LEAHY             apply per-file Leahy amplitude normalisation.
;      RMS1/RMS2/RMS3    optional rms^3-like conversion of BREAL/BIMAG/BMOD.
;      TIME_TOL          RMJD matching tolerance in days. Default 1d-7.
;      FREQ_TOL          frequency matching tolerance in Hz for f1+f2 in
;                         file3. Default 0.25*minimum(df1,df2,df3).
;      RAW_*             return accumulated numerator and denominator sums.
;      POISSON           two-value frequency range used to estimate a constant
;                        noise level separately for each input FFT file.
;      MANUAL_POISSON    scalar noise level applied to all three factors.
;      ZHANG             use GHATS POISSON_ESTIMATE for frequency-dependent
;                        per-factor noise powers.
;      *_INDEP           independent-noise corrected denominator products and
;                        intrinsic cross-bicoherence.  Dead-time cross-band
;                        covariance is not included.
;
; NOTES:
;      The bicoherence denominator is
;
;          den1 = sum |X1(f1) X2(f2)|^2
;          den2 = sum |X3(f1+f2)|^2
;
;      so
;
;          b^2 = |sum X1 X2 X3*|^2 / (den1 den2).
;
;      If INT_BICOH_INDEP is requested, the corrected denominator assumes
;      independent Fourier noise between factors.  Same-band f1=f2 input
;      pixels use the fourth-moment diagonal correction for RAW_DEN1.  This is
;      not a full dead-time covariance-aware correction.
;
; MODIFICATION HISTORY:
;      M. Mendez  25 May 2026  first matched-segment cross-bispectrum
;                              version, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_CROSS_BISPEC2D_MATCH'
   print,'Compute matched-segment 2D cross-bispectrum from three FFT files.'
   print,' '
   print,'Usage:'
   print,"  gh_cross_bispec2d_match,file1,file2,file3,f1,f2,breal,bimag,bmod,bphase,bicoh,nprod,nseg, $"
   print,'                           f1range=[0.061,50],f2range=[0.061,50],/leahy'
   print,' '
   return
endif

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

have_rms = keyword_set(back1) or keyword_set(back2) or keyword_set(back3)
if(have_rms and keyword_set(leahy)) then begin
   massage,'Use either RMS1/RMS2/RMS3 or LEAHY, not both'
   retall
endif
if(have_rms) then begin
   if((~keyword_set(back1)) or (~keyword_set(back2)) or (~keyword_set(back3))) then begin
      massage,'For RMS conversion, give all three: RMS1=, RMS2=, RMS3='
      retall
   endif
endif

;--------------------------------------------------------------------------
; Optional source/noise denominator correction mode.
;--------------------------------------------------------------------------
poimode = 0L
poisson_freq_range = dblarr(2)
if(n_elements(manual_poisson) gt 0) then begin
   poimode = 1L
endif
if(n_elements(poisson_key) eq 2) then begin
   poimode = 2L
   poisson_freq_range = double(poisson_key)
endif
if(keyword_set(zhang)) then poimode = 3L
if(keyword_set(poisson_key) and n_elements(poisson_key) ne 2) then begin
   massage,'POISSON must be a two-element frequency range'
   retall
endif

if(keyword_set(time_tol)) then begin
   tol_time = double(time_tol)
endif else begin
   tol_time = 1.0d-7
endelse

;--------------------------------------------------------------------------
; Open files and read headers.
;--------------------------------------------------------------------------
ghats_openfft,filename1,unit1,/dialog
ghats_openfft,filename2,unit2,/dialog
ghats_openfft,filename3,unit3,/dialog

ntrafos1 = 0L & ntrafos2 = 0L & ntrafos3 = 0L
dummy1 = bytarr(100) & dummy2 = bytarr(100) & dummy3 = bytarr(100)
rmjd01 = 0.0d0 & rmjd02 = 0.0d0 & rmjd03 = 0.0d0
gh_version_string1 = '                ' & gh_version_string2 = '                ' & gh_version_string3 = '                '
observatory1 = '                ' & observatory2 = '                ' & observatory3 = '                '
instrument1  = '                ' & instrument2  = '                ' & instrument3  = '                '
target1      = '                ' & target2      = '                ' & target3      = '                '

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

i_vle1 = fix([dummy1[17:18]],0)+1
i_vle2 = fix([dummy2[17:18]],0)+1
i_vle3 = fix([dummy3[17:18]],0)+1
tdead = 1.0d-5

if(nft_header1 ne nft_header2 or nft_header1 ne nft_header3) then begin
   print,nft_header1,nft_header2,nft_header3
   free_lun,unit1 & free_lun,unit2 & free_lun,unit3
   massage,'FFT lengths are not compatible'
   retall
endif

nft1 = nft_header1/2
nft2 = nft_header2/2
nft3 = nft_header3/2

df1 = double(T1)
df2 = double(T2)
df3 = double(T3)

if(keyword_set(freq_tol)) then begin
   tol_freq = double(freq_tol)
endif else begin
   tol_freq = 0.25d0 * min([df1,df2,df3])
endelse

allfreq1 = (dindgen(nft1)+1.0d0) * df1
allfreq2 = (dindgen(nft2)+1.0d0) * df2
allfreq3 = (dindgen(nft3)+1.0d0) * df3

w1 = where(allfreq1 ge f1range[0] and allfreq1 le f1range[1], nf1)
w2 = where(allfreq2 ge f2range[0] and allfreq2 le f2range[1], nf2)

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

frequency1 = allfreq1[w1]
frequency2 = allfreq2[w2]

;--------------------------------------------------------------------------
; Read all FFTs into memory. This avoids unsafe lockstep assumptions.
;--------------------------------------------------------------------------
print,'Reading FFT files for RMJD matching'
print,'  File1 FFTs: ',ntrafos1,'  df: ',df1,'  rmjd0: ',rmjd01
print,'  File2 FFTs: ',ntrafos2,'  df: ',df2,'  rmjd0: ',rmjd02
print,'  File3 FFTs: ',ntrafos3,'  df: ',df3,'  rmjd0: ',rmjd03
print,'  TIME_TOL days: ',tol_time
print,'  FREQ_TOL Hz  : ',tol_freq

data1 = complexarr(nft1,ntrafos1)
data2 = complexarr(nft2,ntrafos2)
data3 = complexarr(nft3,ntrafos3)
rmjd_arr1 = dblarr(ntrafos1) & rmjd_arr2 = dblarr(ntrafos2) & rmjd_arr3 = dblarr(ntrafos3)
cnts_arr1 = dblarr(ntrafos1) & cnts_arr2 = dblarr(ntrafos2) & cnts_arr3 = dblarr(ntrafos3)
poi_hdr_sum1 = 0.0d0 & poi_hdr_sum2 = 0.0d0 & poi_hdr_sum3 = 0.0d0
npoi_hdr1 = 0L & npoi_hdr2 = 0L & npoi_hdr3 = 0L
power_sum1 = dblarr(nft1) & power_sum2 = dblarr(nft2) & power_sum3 = dblarr(nft3)
poi_zhang_sum1 = dblarr(nft1) & poi_zhang_sum2 = dblarr(nft2) & poi_zhang_sum3 = dblarr(nft3)

rdata1 = complexarr(nft1)
rdata2 = complexarr(nft2)
rdata3 = complexarr(nft3)
poitmp = fltarr(nft1)

for it=0L,ntrafos1-1L do begin
   read_fft_line,unit1,muflag1,rmjd,cnts,poisson,current_vle_rate,fndet,rdata1
   rmjd_arr1[it] = double(rmjd)
   cnts_arr1[it] = double(cnts)
   data1[*,it] = rdata1
   if(keyword_set(leahy)) then power_norm = 2.0d0/double(cnts) else power_norm = 1.0d0
   power_sum1 = power_sum1 + double(abs(rdata1)^2)*power_norm
   if(poimode eq 0L and finite(poisson) and double(poisson) gt 0.0d0) then begin
      poi_hdr_sum1 = poi_hdr_sum1 + double(poisson)*power_norm
      npoi_hdr1 = npoi_hdr1 + 1L
   endif
   if(poimode eq 3L) then begin
      poisson_estimate,poitmp,cnts,1.0d0/df1,current_vle_rate,nft1,fndet, $
                       i_vle1,tdead,/differential
      leahy_power_norm = 2.0d0/double(cnts)
      poi_zhang_sum1 = poi_zhang_sum1 + double(poitmp)*(power_norm/leahy_power_norm)
   endif
endfor

poitmp = fltarr(nft2)
for it=0L,ntrafos2-1L do begin
   read_fft_line,unit2,muflag2,rmjd,cnts,poisson,current_vle_rate,fndet,rdata2
   rmjd_arr2[it] = double(rmjd)
   cnts_arr2[it] = double(cnts)
   data2[*,it] = rdata2
   if(keyword_set(leahy)) then power_norm = 2.0d0/double(cnts) else power_norm = 1.0d0
   power_sum2 = power_sum2 + double(abs(rdata2)^2)*power_norm
   if(poimode eq 0L and finite(poisson) and double(poisson) gt 0.0d0) then begin
      poi_hdr_sum2 = poi_hdr_sum2 + double(poisson)*power_norm
      npoi_hdr2 = npoi_hdr2 + 1L
   endif
   if(poimode eq 3L) then begin
      poisson_estimate,poitmp,cnts,1.0d0/df2,current_vle_rate,nft2,fndet, $
                       i_vle2,tdead,/differential
      leahy_power_norm = 2.0d0/double(cnts)
      poi_zhang_sum2 = poi_zhang_sum2 + double(poitmp)*(power_norm/leahy_power_norm)
   endif
endfor

poitmp = fltarr(nft3)
for it=0L,ntrafos3-1L do begin
   read_fft_line,unit3,muflag3,rmjd,cnts,poisson,current_vle_rate,fndet,rdata3
   rmjd_arr3[it] = double(rmjd)
   cnts_arr3[it] = double(cnts)
   data3[*,it] = rdata3
   if(keyword_set(leahy)) then power_norm = 2.0d0/double(cnts) else power_norm = 1.0d0
   power_sum3 = power_sum3 + double(abs(rdata3)^2)*power_norm
   if(poimode eq 0L and finite(poisson) and double(poisson) gt 0.0d0) then begin
      poi_hdr_sum3 = poi_hdr_sum3 + double(poisson)*power_norm
      npoi_hdr3 = npoi_hdr3 + 1L
   endif
   if(poimode eq 3L) then begin
      poisson_estimate,poitmp,cnts,1.0d0/df3,current_vle_rate,nft3,fndet, $
                       i_vle3,tdead,/differential
      leahy_power_norm = 2.0d0/double(cnts)
      poi_zhang_sum3 = poi_zhang_sum3 + double(poitmp)*(power_norm/leahy_power_norm)
   endif
endfor

free_lun,unit1 & free_lun,unit2 & free_lun,unit3

power_avg1 = power_sum1 / double(ntrafos1)
power_avg2 = power_sum2 / double(ntrafos2)
power_avg3 = power_sum3 / double(ntrafos3)
poisson_power1 = dblarr(nft1) + !values.d_nan
poisson_power2 = dblarr(nft2) + !values.d_nan
poisson_power3 = dblarr(nft3) + !values.d_nan
poisson_method_used = 'NONE'
poisson_level_used = !values.d_nan
compute_indep = (arg_present(int_bicoh_indep) or arg_present(den1_corr_indep) or $
                 arg_present(den2_corr_indep) or arg_present(intr_valid) or $
                 arg_present(corr_flags) or arg_present(diag_flag))

case poimode of
   1L: begin
      poisson_power1 = double(manual_poisson) + dblarr(nft1)
      poisson_power2 = double(manual_poisson) + dblarr(nft2)
      poisson_power3 = double(manual_poisson) + dblarr(nft3)
      poisson_method_used = 'MANUAL_POISSON'
      poisson_level_used = double(manual_poisson)
   end
   2L: begin
      wpoi1 = where(allfreq1 ge poisson_freq_range[0] and $
                    allfreq1 le poisson_freq_range[1],npoi1)
      wpoi2 = where(allfreq2 ge poisson_freq_range[0] and $
                    allfreq2 le poisson_freq_range[1],npoi2)
      wpoi3 = where(allfreq3 ge poisson_freq_range[0] and $
                    allfreq3 le poisson_freq_range[1],npoi3)
      if(npoi1 le 0 or npoi2 le 0 or npoi3 le 0) then begin
         massage,'POISSON range contains no bins in at least one cross-bispectrum factor'
         retall
      endif
      plev1 = total(power_avg1[wpoi1],/double)/double(npoi1)
      plev2 = total(power_avg2[wpoi2],/double)/double(npoi2)
      plev3 = total(power_avg3[wpoi3],/double)/double(npoi3)
      poisson_power1 = plev1 + dblarr(nft1)
      poisson_power2 = plev2 + dblarr(nft2)
      poisson_power3 = plev3 + dblarr(nft3)
      poisson_method_used = 'POISSON_RANGE'
      poisson_level_used = (plev1 + plev2 + plev3)/3.0d0
   end
   3L: begin
      poisson_power1 = poi_zhang_sum1 / double(ntrafos1)
      poisson_power2 = poi_zhang_sum2 / double(ntrafos2)
      poisson_power3 = poi_zhang_sum3 / double(ntrafos3)
      poisson_method_used = 'ZHANG'
      wpn = where(finite(poisson_power1),npn)
      if(npn gt 0) then poisson_level_used = total(poisson_power1[wpn],/double)/double(npn)
   end
   else: begin
      if(npoi_hdr1 gt 0L and npoi_hdr2 gt 0L and npoi_hdr3 gt 0L) then begin
         plev1 = poi_hdr_sum1/double(npoi_hdr1)
         plev2 = poi_hdr_sum2/double(npoi_hdr2)
         plev3 = poi_hdr_sum3/double(npoi_hdr3)
         poisson_power1 = plev1 + dblarr(nft1)
         poisson_power2 = plev2 + dblarr(nft2)
         poisson_power3 = plev3 + dblarr(nft3)
         poisson_method_used = 'HEADER_POISSON'
         poisson_level_used = (plev1 + plev2 + plev3)/3.0d0
      endif
   end
endcase

source_power1 = power_avg1 - poisson_power1
source_power2 = power_avg2 - poisson_power2
source_power3 = power_avg3 - poisson_power3

if(compute_indep and poisson_method_used eq 'NONE') then begin
   massage,'Independent-noise cross-bicoherence requested but no Poisson/noise correction is available'
   retall
endif

;print,'DEBUG filenames: ',filename1,' | ',filename2,' | ',filename3
;print,'DEBUG data diff seg0 1-2, 1-3, 2-3: ',total(abs(data1[*,0]-data2[*,0])),total(abs(data1[*,0]-data3[*,0])),total(abs(data2[*,0]-data3[*,0]))
;print,'DEBUG cnts first: ',cnts_arr1[0],cnts_arr2[0],cnts_arr3[0]
;print,'DEBUG rmjd first: '
;print,format='(3F20.10)',rmjd_arr1[0],rmjd_arr2[0],rmjd_arr3[0]

;--------------------------------------------------------------------------
; Build matched triplets using file1 as the reference.
;--------------------------------------------------------------------------
match1 = lonarr(ntrafos1)
match2 = lonarr(ntrafos1)
match3 = lonarr(ntrafos1)
match1[*] = -1L & match2[*] = -1L & match3[*] = -1L
nmatch = 0L

for i=0L,ntrafos1-1L do begin
   dt2 = abs(rmjd_arr2 - rmjd_arr1[i])
   dt3 = abs(rmjd_arr3 - rmjd_arr1[i])
   mindt2 = min(dt2,j2)
   mindt3 = min(dt3,j3)
   if((mindt2 le tol_time) and (mindt3 le tol_time)) then begin
      match1[nmatch] = i
      match2[nmatch] = j2
      match3[nmatch] = j3
      nmatch = nmatch + 1L
   endif
endfor

if(nmatch le 0L) then begin
   print,'No matched FFT triplets found.'
   print,'Try a larger TIME_TOL, but do not exceed a small fraction of the segment length.'
   return
endif

match1 = match1[0:nmatch-1L]
match2 = match2[0:nmatch-1L]
match3 = match3[0:nmatch-1L]

;print,'Matched FFT triplets: ',nmatch
;print,'DEBUG first matched indices: ',match1[0],match2[0],match3[0]
;print,'DEBUG first matched rmjd: '
;print,format='(3F20.10)',rmjd_arr1[match1[0]],rmjd_arr2[match2[0]],rmjd_arr3[match3[0]]
;print,'DEBUG data diff first matched 1-2, 1-3, 2-3: ',total(abs(data1[*,match1[0]]-data2[*,match2[0]])),total(abs(data1[*,match1[0]]-data3[*,match3[0]])),total(abs(data2[*,match2[0]]-data3[*,match3[0]]))
;
;debug_sum12=0d0 & debug_sum13=0d0 & debug_sum23=0d0
;for kk=0L,nmatch-1L do begin
;   debug_sum12=debug_sum12+total(abs(data1[*,match1[kk]]-data2[*,match2[kk]]))
;   debug_sum13=debug_sum13+total(abs(data1[*,match1[kk]]-data3[*,match3[kk]]))
;   debug_sum23=debug_sum23+total(abs(data2[*,match2[kk]]-data3[*,match3[kk]]))
;endfor
;print,'DEBUG total matched data diff 1-2,1-3,2-3: ',debug_sum12,debug_sum13,debug_sum23


;--------------------------------------------------------------------------
; Data selection on matched file1 index/time and summed count rate.
;--------------------------------------------------------------------------
if(keyword_set(index)) then begin
   firstfft = index[0]
   lastfft  = index[1]
endif else begin
   firstfft = 0L
   lastfft  = 10000000L
endelse

if(keyword_set(time)) then begin
   firsttime = double(time[0])
   lasttime  = double(time[1])
endif else begin
   firsttime = 0.0d0
   lasttime  = 1.0d10
endelse

if(keyword_set(rate)) then begin
   firstrate = double(rate[0])
   lastrate  = double(rate[1])
endif else begin
   firstrate = 0.0d0
   lastrate  = 2.0d9
endelse

if(keyword_set(sel)) then begin
   goodarray = long(sel)
endif else begin
   goodarray = lindgen(ntrafos1)
endelse

wgood = where((goodarray ge firstfft) and (goodarray le lastfft), ngood)
if(ngood le 0) then begin
   massage,'No FFTs found in selected INDEX/SEL range'
   retall
endif
goodarray = goodarray[wgood]

;--------------------------------------------------------------------------
; Allocate outputs and accumulators.
;--------------------------------------------------------------------------
nan = !values.f_nan
bsum = complexarr(nf1,nf2)
den1 = dblarr(nf1,nf2)
den2 = dblarr(nf1,nf2)
nprod_used = lonarr(nf1,nf2)

breal  = fltarr(nf1,nf2) + nan
bimag  = fltarr(nf1,nf2) + nan
bmod   = fltarr(nf1,nf2) + nan
bphase = fltarr(nf1,nf2) + nan
bicoh  = fltarr(nf1,nf2) + nan

if(compute_indep) then begin
   den1_corr_indep = dblarr(nf1,nf2) + !values.d_nan
   den2_corr_indep = dblarr(nf1,nf2) + !values.d_nan
   int_bicoh_indep = fltarr(nf1,nf2) + !values.f_nan
   intr_valid = bytarr(nf1,nf2)
   corr_flags = lonarr(nf1,nf2)
   diag_flag = bytarr(nf1,nf2)
endif

same12 = strcmp(strtrim(string(filename1),2),strtrim(string(filename2),2))

; Precompute sum-frequency index in file3 for each pixel.
isum3_arr = lonarr(nf1,nf2) - 1L
for ii1=0L,nf1-1L do begin
   f1v = frequency1[ii1]
   for ii2=0L,nf2-1L do begin
      fsum = f1v + frequency2[ii2]
      k3 = long(round(fsum/df3)) - 1L
      if(k3 ge 0L and k3 lt nft3) then begin
         if(abs(allfreq3[k3] - fsum) le tol_freq) then isum3_arr[ii1,ii2] = k3
      endif
   endfor
endfor

nvalid_sum = total(isum3_arr ge 0L)
if(nvalid_sum le 0L) then begin
   massage,'No valid f1+f2 pixels found within FREQ_TOL on file3 grid'
   retall
endif

;--------------------------------------------------------------------------
; Accumulate cross-bispectral sums.
;--------------------------------------------------------------------------
print,'Computing matched 2D cross-bispectrum'
print,'F1RANGE file1: ',frequency1[0],' - ',frequency1[nf1-1],' Hz (',nf1,' bins)'
print,'F2RANGE file2: ',frequency2[0],' - ',frequency2[nf2-1],' Hz (',nf2,' bins)'

nseg = 0L
flux1 = 0.0d0 & flux2 = 0.0d0 & flux3 = 0.0d0

for im=0L,nmatch-1L do begin

   i1seg = match1[im]
   i2seg = match2[im]
   i3seg = match3[im]

   gotcha = where(goodarray eq i1seg, ngot)
   if(ngot le 0) then goto,next_match

   t_1 = (rmjd_arr1[i1seg] - rmjd01) * 86400.0d0
   t_2 = t_1 + 1.0d0/df1

   cratem = cnts_arr1[i1seg]*df1 + cnts_arr2[i2seg]*df2 + cnts_arr3[i3seg]*df3

   if((t_1 lt firsttime) or (t_2 gt lasttime)) then goto,next_match
   if((cratem lt firstrate) or (cratem gt lastrate)) then goto,next_match

   if(keyword_set(leahy)) then begin
      if((cnts_arr1[i1seg] le 0.0d0) or (cnts_arr2[i2seg] le 0.0d0) or (cnts_arr3[i3seg] le 0.0d0)) then begin
         massage,'Cannot apply LEAHY normalization: segment has cnts <= 0'
         retall
      endif
      ampnorm1 = sqrt(2.0d0/cnts_arr1[i1seg])
      ampnorm2 = sqrt(2.0d0/cnts_arr2[i2seg])
      ampnorm3 = sqrt(2.0d0/cnts_arr3[i3seg])
   endif else begin
      ampnorm1 = 1.0d0
      ampnorm2 = 1.0d0
      ampnorm3 = 1.0d0
   endelse

   nseg = nseg + 1L
   flux1 = flux1 + cnts_arr1[i1seg]*df1
   flux2 = flux2 + cnts_arr2[i2seg]*df2
   flux3 = flux3 + cnts_arr3[i3seg]*df3

   for ii1=0L,nf1-1L do begin
      i1freq = w1[ii1]
      x1 = data1[i1freq,i1seg] * ampnorm1

      for ii2=0L,nf2-1L do begin
         k3 = isum3_arr[ii1,ii2]
         if(k3 ge 0L) then begin
            i2freq = w2[ii2]

            x2 = data2[i2freq,i2seg] * ampnorm2
            x3 = data3[k3,i3seg]     * ampnorm3

            z = x1 * x2 * conj(x3)
            d1this = double(abs(x1*x2)^2)
            d2this = double(abs(x3)^2)

            bsum[ii1,ii2] = bsum[ii1,ii2] + z
            den1[ii1,ii2] = den1[ii1,ii2] + d1this
            den2[ii1,ii2] = den2[ii1,ii2] + d2this
            nprod_used[ii1,ii2] = nprod_used[ii1,ii2] + 1L
         endif
      endfor
   endfor

   next_match:
endfor

if(nseg le 0L) then begin
   massage,'No matched FFT triplets passed the selection'
   retall
endif

flux1 = flux1/double(nseg)
flux2 = flux2/double(nseg)
flux3 = flux3/double(nseg)

print,'  ',strtrim(string(nseg),1),' matched FFT triplets selected'

;--------------------------------------------------------------------------
; Derived products.
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
            bicoh[ii1,ii2] = float(abs(bsum[ii1,ii2])^2/(den1[ii1,ii2]*den2[ii1,ii2]))
         endif
         if(compute_indep) then begin
            i1freq = w1[ii1]
            i2freq = w2[ii2]
            k3 = isum3_arr[ii1,ii2]
            if(k3 ge 0L) then begin
               diagonal = same12 and (abs(allfreq1[i1freq]-allfreq2[i2freq]) le tol_freq)
               d1c = gh_bispec_den1_corr_indep(den1[ii1,ii2],nprod_used[ii1,ii2], $
                                                source_power1[i1freq],poisson_power1[i1freq], $
                                                source_power2[i2freq],poisson_power2[i2freq], $
                                                same12,diagonal,flag=flag1)
               d2c = gh_bispec_den2_corr_indep(nprod_used[ii1,ii2], $
                                                source_power3[k3],flag=flag2)
               den1_corr_indep[ii1,ii2] = d1c
               den2_corr_indep[ii1,ii2] = d2c
               cflag = flag1
               if((flag2 and 4L) ne 0L) then cflag = cflag + 4L
               if(((flag2 and 8L) ne 0L) and ((cflag and 8L) eq 0L)) then cflag = cflag + 8L
               corr_flags[ii1,ii2] = cflag
               if(diagonal) then diag_flag[ii1,ii2] = 1B
               if((d1c gt 0.0d0) and (d2c gt 0.0d0)) then begin
                  int_bicoh_indep[ii1,ii2] = float(abs(bsum[ii1,ii2])^2/(d1c*d2c))
                  intr_valid[ii1,ii2] = 1B
               endif
            endif
         endif
      endif
   endfor
endfor

;--------------------------------------------------------------------------
; Optional rms^3 conversion of averaged bispectrum products.
; Raw sums remain in internal units.
;--------------------------------------------------------------------------
if(have_rms) then begin
   if(back1 ge flux1) then begin
      massage,'Error: background flux higher than source+bkg flux in file1'
      retall
   endif
   if(back2 ge flux2) then begin
      massage,'Error: background flux higher than source+bkg flux in file2'
      retall
   endif
   if(back3 ge flux3) then begin
      massage,'Error: background flux higher than source+bkg flux in file3'
      retall
   endif

   norm_bispec = (2.0d0^1.5d0) * ((df1*df2*df3)^0.5d0) / $
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
raw_den1 = den1
raw_den2 = den2

end
