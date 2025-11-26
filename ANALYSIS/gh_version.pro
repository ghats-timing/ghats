pro gh_version
;+
; NAME:
;      MU_VERSION
; PURPOSE:
;      Prints to terminal the current MU version
; EXPLANATION:
;      This procedure prints to the terminal information about the
;      MU version.
;
; CALLING SEQUENCE:
;       MU_VERSION
; INPUTS:
;       NONE
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       NONE
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS:
;       None
; ROUTINES USED:
;       NONE
; NOTES:
;       NONE
; MODIFICATION HISTORY:
;       T. Belloni  12 Nov 2001  implementation
;		T. Belloni  01 Dec 2010  from mu6. Version from common block
;-
;--------------------------------------------------------------------------
common vers, versione, data_versione
;versione = 'GH Version 0.0.3'
;day  = '2010 Dec 01'
stringa = [versione, data_versione]

result=dialog_message(stringa,title='Welcome to GH!',/information,/center)

end
