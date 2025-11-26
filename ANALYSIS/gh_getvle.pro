pro gh_getvle,filename,i_vle,nodialog=nodialog
;+
; NAME:
;      MU_GETVLE
; PURPOSE:
;      Reads is the VLE time parameter from a HK file (FH5a)
; EXPLANATION:
;      This procedure reads in the VLE time parameter from a HK
;      file. 
;
; CALLING SEQUENCE:
;       GH_GETVLE,FILENAME,I_VLE[,/NODIALOG]
; INPUTS:
;       FILENAME = name of the input FH5a file
;
; OUTPUTS:
;       I_VLE    = Output VLE parameter
;
; KEYWORDS:
;       NODIALOG = If set, no dialog window is displayed for file opening
;
; EXAMPLE:
;       None
;
; COMMON BLOCKS:
;       None
; ROUTINES USED:
;       None
; NOTES:
;       None
; MODIFICATION HISTORY:
;       T. Belloni  17 Dev 2002  implementation
;-
;--------------------------------------------------------------------------
;
;  Reads ivle information from FH5a
;
if(not keyword_set(nodialog)) then begin
   filename=dialog_pickfile(filter='FH5a*',/read)
endif
unit = 11
fxbopen,unit,filename,1,header,errmsg=errmsg
;
; Read in first row only
;
fxbread,unit,i_vle,'dsVle   ',1
fxbclose,unit

end
