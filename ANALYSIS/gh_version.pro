pro gh_version,help=help
;+
; NAME:
;      GH_VERSION
; PURPOSE:
;      Print the current GHATS version.
; EXPLANATION:
;      This procedure reports the GHATS version stored in the version
;      common block.
;
; CALLING SEQUENCE:
;       GH_VERSION
;       GH_VERSION,/HELP
; INPUTS:
;       NONE
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       HELP = If set, print usage information and return.
;
; EXAMPLE:
;       GH_VERSION
;       GH_VERSION,/HELP
;
; COMMON BLOCKS:
;       vers
; ROUTINES USED:
;       NONE
; NOTES:
;       NONE
; MODIFICATION HISTORY:
;       T. Belloni  12 Nov 2001  implementation
;		T. Belloni  01 Dec 2010  from mu6. Version from common block
;       M. Mendez/Codex  17 Jul 2026  added /HELP
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_VERSION'
   print,''
   print,'Print the current GHATS version string.'
   print,''
   print,'Usage:'
   print,'  GH_VERSION'
   print,'  GH_VERSION,/HELP'
   print,''
   print,'Output:'
   print,'  Shows the version and release date stored in the GHATS version common block.'
   print,''
   return
endif

common vers, versione, data_versione
;versione = 'GH Version 0.0.3'
;day  = '2010 Dec 01'
stringa = [versione, data_versione]

result=dialog_message(stringa,title='Welcome to GH!',/information,/center)

end
