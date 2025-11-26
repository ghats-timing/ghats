pro mu_ask,filename,patte
;
; Procedure for prompting for a file name
; MU 3.0 version, spawn-free and with a new flavor
; MU 5.0 version, adapted for GDL
;
;-------------------------------------------------------------------------
; Parameters
;
; filename                 O: selected filename
;-------------------------------------------------------------------------
;
;       T. Belloni  13 Aug 2008  title added and check for empty string
;       T. Belloni  03 Feb 2009  adepted for GDL
;
common sis, sistema   ; common block with system variable
;patte  = 'FS37*;FS3b*;FS3f*;FS4f*;@*;e*'
gpatte = '*'

if(sistema eq 'IDL') then begin
   filename=dialog_pickfile(filter=patte,title='Select input file for FFT', /read)
 endif else begin
   listona = file_search(gpatte)
   nfiles  = n_elements(listona)
   numeri  = indgen(nfiles)+1
   for i=0,nfiles-1 do begin
	   print,numeri(i),': ',listona(i)
   endfor
   print,'--------------------------------------------'
   print,format='($,"(INPUT FILE) -> ")'
   read,chout
   scelto = fix(chout)-1
   filename=listona(scelto)
endelse

if(filename eq '') then begin
   filename = ' '
endif

return
end
