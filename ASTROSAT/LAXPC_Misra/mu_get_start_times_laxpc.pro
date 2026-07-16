pro mu_get_start_times_laxpc,      $
               filenames,nfiles,nmetafiles,    $
	       tstarts
;
;  Obtain start times (TSTART keyword) from input files
;
;-------------------------------------------------------------------------
; Parameters
;
; filenames                    I: array with input filenames
; nfiles                       I: X dimension of filenames
; nmetafiles                   I: Y dimension of filenames
;
; tstarts                      O: array of start times from keywords
;-------------------------------------------------------------------------
for j=0,nmetafiles-1 do begin
   for i=0,nfiles-1 do begin
      unit = 10
      errmsg=''
      ;fxbopen,unit,filenames(i,j),1,header,errmsg=errmsg
	  header=headfits(filenames[i,j])     ; reads in primary header now
      if (errmsg ne '') then begin
	       massage,'Cannot open FITS file for reading'
	       retall
      endif
      tstarts(i,j) = fxpar(header,'TSTART')
      ;fxbclose,unit
   endfor
endfor
end
