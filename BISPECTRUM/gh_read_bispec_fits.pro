pro gh_read_bispec_fits,filename,f1,f2,bsr,bsi,den1,den2,nprod, $
                         breal=breal,bimag=bimag,bmod=bmod, $
                         bphase=bphase,bicoh=bicoh, $
                         int_bicoh=int_bicoh, $
                         den1_corr=den1_corr_kw, $
                         den2_corr=den2_corr_kw, $
                         product=product, $
                         srcfile1=srcfile1,srcfile2=srcfile2,srcfile3=srcfile3, $
                         help=help

if keyword_set(help) then begin
   print,'GH_READ_BISPEC_FITS'
   print,'Reads saved bispectrum FITS arrays into memory.'
   print,'Usage: gh_read_bispec_fits,file,f1,f2,bsr,bsi,den1,den2,nprod'
   print,'Optional keyword outputs: breal=, bimag=, bmod=, bphase=, bicoh='
   print,'                          int_bicoh=, den1_corr=, den2_corr='
   print,'Optional provenance outputs: product=, srcfile1=, srcfile2=, srcfile3='
   return
endif

primary_header = headfits(filename)
product = fxpar(primary_header,'PRODUCT')
srcfile1 = fxpar(primary_header,'SRCFILE1')
srcfile2 = fxpar(primary_header,'SRCFILE2')
srcfile3 = fxpar(primary_header,'SRCFILE3')

nextend = long(fxpar(primary_header,'NEXTEND'))
if(nextend le 0L) then begin
   rawsum_flag = long(fxpar(primary_header,'RAWSUMS'))
   if(rawsum_flag ne 0L) then begin
      nextend = 12L
   endif else begin
      nextend = 8L
   endelse
endif

;--------------------------------------------------------------------------
; Preferred image-extension format. Read by EXTNAME so adding new extensions
; does not change the meaning of older plotting extension numbers.
;--------------------------------------------------------------------------
got_f1 = 0B
got_f2 = 0B
got_bsr = 0B
got_bsi = 0B
got_den1 = 0B
got_den2 = 0B
got_nprod = 0B
got_named = 0B

for iext=1L,nextend do begin
   data = mrdfits(filename,iext,h,/silent,status=status)
   if(status ne 0) then goto,done_named_extensions

   extname = strupcase(strtrim(sxpar(h,'EXTNAME'),2))

   case extname of
      'FREQ1': begin
         f1 = data
         got_f1 = 1B
         got_named = 1B
      end
      'FREQ2': begin
         f2 = data
         got_f2 = 1B
         got_named = 1B
      end
      'RAW_BSUM_REAL': begin
         bsr = data
         got_bsr = 1B
         got_named = 1B
      end
      'RAW_BSUM_IMAG': begin
         bsi = data
         got_bsi = 1B
         got_named = 1B
      end
      'RAW_DEN1': begin
         den1 = data
         got_den1 = 1B
         got_named = 1B
      end
      'RAW_DEN2': begin
         den2 = data
         got_den2 = 1B
         got_named = 1B
      end
      'NPROD_USED': begin
         nprod = data
         got_nprod = 1B
         got_named = 1B
      end
      'BREAL': begin
         breal = data
         got_named = 1B
      end
      'BIMAG': begin
         bimag = data
         got_named = 1B
      end
      'BMOD': begin
         bmod = data
         got_named = 1B
      end
      'BPHASE': begin
         bphase = data
         got_named = 1B
      end
      'BICOH': begin
         bicoh = data
         got_named = 1B
      end
      'INT_BICOH': begin
         int_bicoh = data
         got_named = 1B
      end
      'BICOH_INT_OBS': begin
         int_bicoh = data
         got_named = 1B
      end
      'RAW_DEN1_CORR': begin
         den1_corr_kw = data
         got_named = 1B
      end
      'RAW_DEN2_CORR': begin
         den2_corr_kw = data
         got_named = 1B
      end
      else: begin
      end
   endcase
endfor

done_named_extensions:

if(got_named and got_f1 and got_f2) then begin
   if((~got_bsr) or (~got_bsi) or (~got_den1) or (~got_den2)) then begin
      bsr = !NULL
      bsi = !NULL
      den1 = !NULL
      den2 = !NULL
      print,'Warning: raw bispectrum sums are not present in this FITS file.'
   endif

   if(~got_nprod) then nprod = !NULL

   print,'Read bispectrum from image-extension FITS format.'
   return
endif

;--------------------------------------------------------------------------
; Older binary-table compatibility path.
;--------------------------------------------------------------------------
data=mrdfits(filename,1,h,/silent,status=status)

if status eq 0 then begin
   if(size(data,/type) eq 8) then begin
      tags=strupcase(tag_names(data))

      if total(tags eq 'FREQUENCY1') gt 0 then f1=data.frequency1
      if total(tags eq 'F1') gt 0 then f1=data.f1

      if total(tags eq 'FREQUENCY2') gt 0 then f2=data.frequency2
      if total(tags eq 'F2') gt 0 then f2=data.f2

      if total(tags eq 'RAW_BSUM_REAL') gt 0 then bsr=data.raw_bsum_real
      if total(tags eq 'BSR') gt 0 then bsr=data.bsr

      if total(tags eq 'RAW_BSUM_IMAG') gt 0 then bsi=data.raw_bsum_imag
      if total(tags eq 'BSI') gt 0 then bsi=data.bsi

      if total(tags eq 'RAW_DEN1') gt 0 then den1=data.raw_den1
      if total(tags eq 'DEN1') gt 0 then den1=data.den1

      if total(tags eq 'RAW_DEN2') gt 0 then den2=data.raw_den2
      if total(tags eq 'DEN2') gt 0 then den2=data.den2

      if total(tags eq 'NPROD') gt 0 then nprod=data.nprod
      if total(tags eq 'NPROD_USED') gt 0 then nprod=data.nprod_used

      if total(tags eq 'BREAL') gt 0 then breal=data.breal
      if total(tags eq 'BIMAG') gt 0 then bimag=data.bimag
      if total(tags eq 'BMOD') gt 0 then bmod=data.bmod
      if total(tags eq 'BPHASE') gt 0 then bphase=data.bphase
      if total(tags eq 'BICOH') gt 0 then bicoh=data.bicoh

      print,'Read bispectrum from binary-table FITS extension.'
      return
   endif
endif

;--------------------------------------------------------------------------
; Last fallback: old image-extension layout assumed by the original reader.
; ext 1=f1, 2=f2, 3=bsr, 4=bsi, 5=den1, 6=den2, 7=nprod.
;--------------------------------------------------------------------------
f1=mrdfits(filename,1,h1,/silent)
f2=mrdfits(filename,2,h2,/silent)
bsr=mrdfits(filename,3,h3,/silent)
bsi=mrdfits(filename,4,h4,/silent)
den1=mrdfits(filename,5,h5,/silent)
den2=mrdfits(filename,6,h6,/silent)
nprod=mrdfits(filename,7,h7,/silent)

print,'Read bispectrum from legacy numbered image-extension FITS format.'

end
