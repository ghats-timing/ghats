pro gh_cross_bispec_target_segments,filename1,filename2,filename3, $
                                   band_edges,target_c,target_w, $
                                   zseg,d1seg,d2seg,npix_target,seg_rmjds,nseg, $
                                   f1=f1,f2=f2, $
                                   f1range=f1range,f2range=f2range, $
                                   index=index,time=time,sel=sel,rate=rate, $
                                   rms1=back1,rms2=back2,rms3=back3,leahy=leahy, $
                                   matched_i1=matched_i1, $
                                   matched_i2=matched_i2, $
                                   matched_i3=matched_i3, $
                                   time_tol=time_tol,freq_tol=freq_tol, $
                                   help=help
;+
; NAME:
;      GH_CROSS_BISPEC_TARGET_SEGMENTS
;
; PURPOSE:
;      Computes target-level segment sums for matched-segment cross-bispectra.
;
; EXPLANATION:
;      For each selected matched FFT triplet this routine computes
;
;          z_s = sum X1_s(f1) X2_s(f2) X3_s*(f1+f2)
;
;      inside each requested target region, together with the denominator
;      terms used for bicoherence. It stores only [ntarget,nseg] arrays,
;      not [nf1,nf2,nseg] cubes.
;
; CALLING SEQUENCE:
;      GH_CROSS_BISPEC_TARGET_SEGMENTS,FILE1,FILE2,FILE3, $
;                                      BAND_EDGES,TARGET_C,TARGET_W, $
;                                      ZSEG,D1SEG,D2SEG,NPIX,SEG_RMJDS,NSEG, $
;                                      F1RANGE=F1RANGE,F2RANGE=F2RANGE, [/LEAHY]
;
; KEYWORDS:
;      F1RANGE/F2RANGE  required frequency ranges.
;      INDEX/TIME/SEL/RATE selection is applied to the matched file1 index,
;                         file1-relative time, and summed count rate.
;      LEAHY             apply per-file Leahy amplitude normalisation.
;      RMS1/RMS2/RMS3    accepted for interface consistency. The target
;                         denominator-normalized outputs are unchanged.
;      TIME_TOL          RMJD matching tolerance in days. Default 1d-7.
;      FREQ_TOL          frequency matching tolerance in Hz for f1+f2 in
;                         file3. Default 0.25*minimum(df1,df2,df3).
;      F1/F2             return selected file1/file2 frequency grids.
;      MATCHED_I1/I2/I3  original segment indices in file1/file2/file3.
;
; OUTPUTS:
;      ZSEG              complex target sums, [ntarget,nseg].
;      D1SEG             sum |X1(f1) X2(f2)|^2, [ntarget,nseg].
;      D2SEG             sum |X3(f1+f2)|^2, [ntarget,nseg].
;      NPIX_TARGET       number of selected pixels per target region.
;      SEG_RMJDS         RMJD of selected matched file1 segments.
;      NSEG              number of selected matched FFT triplets.
;
; MODIFICATION HISTORY:
;      M. Mendez  31 May 2026  target-level segment sums for bootstrap,
;                              with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_CROSS_BISPEC_TARGET_SEGMENTS'
   print,'Compute target-level matched-segment cross-bispectral sums.'
   print,' '
   print,'Usage:'
   print,"  gh_cross_bispec_target_segments,file1,file2,file3,band_edges,target_c,target_w, $"
   print,'                                  zseg,d1seg,d2seg,npix,seg_rmjds,nseg, $'
   print,'                                  f1range=[0.061,50],f2range=[0.061,50],/leahy'
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
if(n_elements(band_edges) lt 2) then begin
   massage,'BAND_EDGES must contain at least two elements'
   retall
endif
if(n_elements(target_c) ne n_elements(target_w)) then begin
   massage,'TARGET_C and TARGET_W must have the same number of elements'
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

f1 = allfreq1[w1]
f2 = allfreq2[w2]

;--------------------------------------------------------------------------
; Read all FFTs into memory. This follows GH_CROSS_BISPEC2D_MATCH.
;--------------------------------------------------------------------------
print,'Reading FFT files for target-segment RMJD matching'
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

rdata1 = complexarr(nft1)
rdata2 = complexarr(nft2)
rdata3 = complexarr(nft3)

for it=0L,ntrafos1-1L do begin
   read_fft_line,unit1,muflag1,rmjd,cnts,poisson,current_vle_rate,fndet,rdata1
   rmjd_arr1[it] = double(rmjd)
   cnts_arr1[it] = double(cnts)
   data1[*,it] = rdata1
endfor

for it=0L,ntrafos2-1L do begin
   read_fft_line,unit2,muflag2,rmjd,cnts,poisson,current_vle_rate,fndet,rdata2
   rmjd_arr2[it] = double(rmjd)
   cnts_arr2[it] = double(cnts)
   data2[*,it] = rdata2
endfor

for it=0L,ntrafos3-1L do begin
   read_fft_line,unit3,muflag3,rmjd,cnts,poisson,current_vle_rate,fndet,rdata3
   rmjd_arr3[it] = double(rmjd)
   cnts_arr3[it] = double(cnts)
   data3[*,it] = rdata3
endfor

free_lun,unit1 & free_lun,unit2 & free_lun,unit3

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
; Precompute valid target pixels.
;--------------------------------------------------------------------------
nt = n_elements(target_c)
nb = n_elements(band_edges) - 1L
npix_target = lonarr(nt)
target_mask = bytarr(nf1,nf2,nt)
isum3_arr = lonarr(nf1,nf2) - 1L

for ii1=0L,nf1-1L do begin
   ff1 = double(f1[ii1])
   if(ff1 lt band_edges[0] or ff1 ge band_edges[nb]) then continue

   for ii2=0L,nf2-1L do begin
      ff2 = double(f2[ii2])
      if(ff2 lt band_edges[0] or ff2 ge band_edges[nb]) then continue

      fsum = ff1 + ff2
      k3 = long(round(fsum/df3)) - 1L
      if(k3 lt 0L or k3 ge nft3) then continue
      if(abs(allfreq3[k3] - fsum) gt tol_freq) then continue

      for it=0L,nt-1L do begin
         if(abs(fsum - double(target_c[it])) le 0.5d0*double(target_w[it])) then begin
            target_mask[ii1,ii2,it] = 1B
            isum3_arr[ii1,ii2] = k3
            npix_target[it] = npix_target[it] + 1L
         endif
      endfor
   endfor
endfor

if(total(npix_target) le 0L) then begin
   massage,'No target pixels found within requested target windows'
   retall
endif

;--------------------------------------------------------------------------
; Accumulate target-level segment sums.
;--------------------------------------------------------------------------
print,'Computing target-level matched cross-bispectrum segment sums'
print,'F1RANGE file1: ',f1[0],' - ',f1[nf1-1],' Hz (',nf1,' bins)'
print,'F2RANGE file2: ',f2[0],' - ',f2[nf2-1],' Hz (',nf2,' bins)'
print,'NPIX_TARGET:'
print,npix_target

zseg = dcomplexarr(nt,nmatch)
d1seg = dblarr(nt,nmatch)
d2seg = dblarr(nt,nmatch)
seg_rmjds = dblarr(nmatch) + !values.d_nan
matched_i1 = lonarr(nmatch) - 1L
matched_i2 = lonarr(nmatch) - 1L
matched_i3 = lonarr(nmatch) - 1L

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

   iseg_out = nseg
   nseg = nseg + 1L
   flux1 = flux1 + cnts_arr1[i1seg]*df1
   flux2 = flux2 + cnts_arr2[i2seg]*df2
   flux3 = flux3 + cnts_arr3[i3seg]*df3

   seg_rmjds[iseg_out]  = rmjd_arr1[i1seg]
   matched_i1[iseg_out] = i1seg
   matched_i2[iseg_out] = i2seg
   matched_i3[iseg_out] = i3seg

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

            for itarg=0L,nt-1L do begin
               if(target_mask[ii1,ii2,itarg] ne 0B) then begin
                  zseg[itarg,iseg_out] = zseg[itarg,iseg_out] + dcomplex(float(z),imaginary(z))
                  d1seg[itarg,iseg_out] = d1seg[itarg,iseg_out] + d1this
                  d2seg[itarg,iseg_out] = d2seg[itarg,iseg_out] + d2this
               endif
            endfor
         endif
      endfor
   endfor

   next_match:
endfor

if(nseg le 0L) then begin
   massage,'No matched FFT triplets passed the selection'
   retall
endif

zseg = zseg[*,0:nseg-1L]
d1seg = d1seg[*,0:nseg-1L]
d2seg = d2seg[*,0:nseg-1L]
seg_rmjds = seg_rmjds[0:nseg-1L]
matched_i1 = matched_i1[0:nseg-1L]
matched_i2 = matched_i2[0:nseg-1L]
matched_i3 = matched_i3[0:nseg-1L]

print,'  ',strtrim(string(nseg),1),' matched FFT triplets selected'

if(have_rms) then begin
   flux1 = flux1/double(nseg)
   flux2 = flux2/double(nseg)
   flux3 = flux3/double(nseg)
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
   print,'Target segment sums kept in denominator-normalized bootstrap units'
endif else begin
   if(keyword_set(leahy)) then begin
      print,'Target segment sums use Leahy-normalized FFT amplitudes'
   endif else begin
      print,'Target segment sums use native FFT amplitudes'
   endelse
endelse

end
