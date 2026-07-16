pro mu_get_time_resolution_huiyan,filenames,nfiles,nmetafiles,tress
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
   extname   =strtrim(fxpar(header,'EXTNAME'))
   datamode  =strtrim(fxpar(header,'DATAMODE'))
   type=''
   SWITCH instrument OF
      'LE': BEGIN
		      type = 'Huiyan_LE'
			  tress(0,j) = 1.0e-6
			  BREAK
		    END
      'ME': BEGIN
		      type = 'Huiyan_ME'
			  tress(0,j) = 6.0e-6
			  BREAK
		    END
      'HE': BEGIN
		      type = 'Huiyan_HE'
			  tress(0,j) = 2.0e-6
			  BREAK
		    END
	  ELSE: BEGIN
	      massage,'Unrecognized data file!'
	      retall
	  END
   ENDSWITCH
   fxbclose,unit	; close FITS file (can give problems?)
endfor
end
