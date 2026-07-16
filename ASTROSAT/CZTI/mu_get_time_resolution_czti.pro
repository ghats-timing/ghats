pro mu_get_time_resolution_czti,filenames,nfiles,nmetafiles,tress
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
unit = 10
for j=0,nmetafiles-1 do begin
   fxbopen,unit,filenames(0,j),1,header,errmsg=errmsg
   instrument=strtrim(fxpar(header,'INSTRUME'))
;   extname   =strtrim(fxpar(header,'EXTNAME'))
;   datamode  =strtrim(fxpar(header,'DATAMODE'))

   type=''
   if((instrument eq 'CZTI')) then begin
	type = 'czti'
   endif


   if(type eq '') then begin
      print,type
      massage,'Unrecognized data file!'
      retall
   endif

   timedel   =fxpar(header,'TIMEDEL')

   tress(0,j) = timedel
   fxbclose,unit	; close FITS file (can give problems?)
endfor
end
