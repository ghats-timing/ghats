pro mu_get_time_resolution_xmm_new,filenames,nfiles,nmetafiles,tress
;
; Procedure to extract time resolution info
;
;------------------------------------------------------------------------
; Parameters
;
; filenames                 I: array of input filenames
; nfiles                    I: X dimension of filenames
; nmetafiles                I: Y dimension of filenames
; tress                     O: output time resolution
;------------------------------------------------------------------------
;
type=''
unit  = 10
unit2 = 11
for j=0,nmetafiles-1 do begin
   fxbopen,unit,filenames(0,j),1,header,errmsg=errmsg
   instrument=strtrim(fxpar(header,'INSTRUME'))
   extname   =strtrim(fxpar(header,'EXTNAME'))
   datamode  =strtrim(fxpar(header,'DATAMODE'))
fxbclose,unit	; close FITS file (can give problems?)

   type=''
   if((strmid(instrument,0,3) eq 'EPN')) then begin
;     Complicated extraction for pn
	  type = 'epn'
	  ;timedel   = fxpar(header,'FRMTIME')/1000.0D0
	
	  rdfits_struct,filenames(0,0),struttura,/SILENT,/HEADER_ONLY
	  ntags = n_tags(struttura)
	  exposure_extension = -1
;     Now look for the EXPOSU* keyword
	  for k = 0,ntags-1 do begin
		 testa   = struttura(k)    ; testa is a string array with the header
		 trovato = total(strmatch(testa,'*EXPOSU*')) ; if 1, found
		 IF(trovato EQ 1) THEN exposure_extension = k
	  endfor
	  if(exposure_extension eq -1) then begin
		massage,'Time resolution of pn file not found!'
		retall
	  endif
;     Found extension
      fxbopen,unit,filenames(0,0),k,header,errmsg=errmsg
      timedel   = fxpar(header,'EXPOSU*')/1000.0D0  ; does it work?
   endif
  if(strmid(instrument,1,3) eq 'MOS') then begin
;  if((instrument eq 'EMOS1') or (instrument eq 'EMOS2')) then begin
	type = 'mos'
	fxbopen,unit2,filenames(0,j),3,header2,errmsg=errmsg
	timedel   = fxpar(header2,'FRMTIME')/1000.0D0
	fxbclose,unit2
   endif

   if(type eq '') then begin
      print,type
      massage,'Unrecognized data file!'
      retall
   endif

   tress(0,j) = timedel
endfor
end
