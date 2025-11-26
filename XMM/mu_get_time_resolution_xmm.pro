pro mu_get_time_resolution_xmm,filenames,nfiles,nmetafiles,tress
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
   if((instrument eq 'EPN')) then begin
	type = 'epn'
	timedel   = fxpar(header,'FRMTIME')/1000.0D0
	;
	; ****** WARNING: TO REMOVE!
	;
	;timedel=5.965e-3   for Timing mode
   endif
  if((instrument eq 'EMOS1') or (instrument eq 'EMOS2')) then begin
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
