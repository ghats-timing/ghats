PRO mu_read_user_gti_file, filename,gt1,gt2
;+
; NAME:
;      MU_READ_USER_GTI_FILE
; PURPOSE:
;      Obtain the GTIs from a user supplied ASCII file
; EXPLANATION:
;      The ASCII file must be in a format from TIMETRANS (for RXTE)
;
; CALLING SEQUENCE:
;       MU_READ_USER_GTI_FILE,filename,gt1.gt2
; INPUTS:
;       FILENAME = ASCII file to read
; OUTPUTS:
;       GT1: array with start times
;       GT2: array with end times
;
; EXAMPLE:
;       NONE
; COMMON BLOCKS:
;       None
; ROUTINES USED:
;
; NOTES
;       None
; MODIFICATION HISTORY:
;       T. Belloni  11 Mar 2009
;		T. Belloni  09 Jun 2009 possible FITS input
;		T. Belloni  21 Apr 2022 added TIMEZERO reading for OGIP compatibility 
;
; Determine how many lines are in the file
;
; Find out whether it is a FITS file by trying to open it
;-			
   errmsg = ''
   fxbopen,lun,filename,1,header,errmsg=errmsg
   errore=errmsg
   fxbclose,lun
   if(strlen(errore) gt 0) then begin
	; It's not a FITS file, it must be ASCII
    ngti = mu_file_lines(filename)
    OPENR, unit, filename, /GET_LUN  
    gt1 = dblarr(ngti)
    gt2 = gt1
	g1 = 0d0
	g2 = g1
    for i=0,ngti-1 do begin
        readf,unit,g1,g2
        gt1(i) = g1
        gt2(i) = g2
    endfor 
    FREE_LUN, unit
   endif else begin
	; It's a FITS file.
	tzero=fxpar(header,'TIMEZERO')   ; added 21 Apr 2022
	fxbopen,lun,filename,1,header,errmsg=errmsg
	fxbreadm,lun,[1,2],gt1,gt2
	fxbclose,lun
	; Now add TIMEZERO to the GTI values (added 21 Apr 2022)
	gt1 = gt1 + tzero
	gt2 = gt2 + tzero
		
   endelse
   
END