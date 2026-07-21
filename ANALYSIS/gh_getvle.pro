pro gh_getvle,filename,i_vle,nodialog=nodialog,help=help
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
;       GH_GETVLE,FILENAME,I_VLE[,/NODIALOG][,/HELP]
; INPUTS:
;       FILENAME = name of the input FH5a file
;
; OUTPUTS:
;       I_VLE    = Output VLE parameter
;
; KEYWORDS:
;       NODIALOG = If set, no dialog window is displayed for file opening
;       HELP     = If set, print usage information and return
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
if(keyword_set(help)) then begin
   print,''
   print,'GH_GETVLE'
   print,''
   print,'Read the RXTE/PCA VLE time parameter from an FH5a housekeeping file.'
   print,''
   print,'Usage:'
   print,"  GH_GETVLE, 'FH5a_file', i_vle, /NODIALOG"
   print,'  GH_GETVLE, filename, i_vle'
   print,''
   print,'Output:'
   print,'  i_vle  VLE parameter read from the dsVle column'
   print,''
   print,'Keywords: /NODIALOG suppresses file picker; /HELP prints this message.'
   print,''
   return
endif
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
