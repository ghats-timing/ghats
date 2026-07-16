pro mu_get_time_resolution_laxpc,filenames,nfiles,nmetafiles,tress
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
   ;fxbopen,unit,filenames(0,j),1,header,errmsg=errmsg
   header=headfits(filenames[0,j])     ; reads in primary header now
   instrument=strtrim(fxpar(header,'INSTRUME'))
   extname   =strtrim(fxpar(header,'EXTNAME'))
   datamode  =strtrim(fxpar(header,'DATAMODE'))

   type=''
   ;if(((instrument eq 'LAXPC1') or (instrument eq 'LAXPC2') or (instrument eq 'LAXPC3'))$
   ;    and (datamode eq 'EA')) then begin
	;type = 'ea'
   ;endif
   
   if((instrument eq 'LAXPC1') or (instrument eq 'LAXPC2') or (instrument eq 'LAXPC3') or $
   (instrument eq 'LAXPC')) then begin
	type = 'ea'
   endif


   if(type eq '') then begin
      print,type
      massage,'Unrecognized data file!'
      retall
   endif
   fxbopen,unit,filenames(0,j),1,header2,errmsg=errmsg
   timedel   =fxpar(header2,'TIMEDEL')
   fxbclose,unit

   tress(0,j) = timedel
   ;fxbclose,unit	; close FITS file (can give problems?)
endfor
end
