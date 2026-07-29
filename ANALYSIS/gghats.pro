;pro ghats
;+
; NAME: 
;      GHATS
; PURPOSE: 
;      Starts the GHATS environment
; EXPLANATION:
;      This procedure changes the prompt to Ghats> and
;      prints the Ghats welcome message
;
; CALLING SEQUENCE: 
;       GHATS
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
;       Start the current version of Ghats:
;
;       IDL> ghats
;
;
;-----------------------------------------------------------------
;    Welcome to GHATS 3.3.1      28 July 2026
;-----------------------------------------------------------------
;
;
; COMMON BLOCKS: 
;       sis   : for saving the IDL/GDL flag
; ROUTINES USED: 
;       None
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  11 Nov 2009 first testing
;       T. Belloni  03 Jan 2012 release version
;       M. Mendez   29 Aug 2023 release version
;       M. Mendez   28 Jul 2026 release version 3.3.1
;-

common sis, sistema
common vers, versione, data_versione
common versione,version_id
sistema = 'GDL'
versione = 'GHATS  V3.3.1'
data_versione = '28 July 2026'
version_id = 'GHATSIB0331     '


;--------------------------------------------------------------------------
!PROMPT = 'Ghats> '
print,' '
print,'-----------------------------------------------------------------'
print,'           Welcome to ',versione,'      ',data_versione
print,'-----------------------------------------------------------------'
print,' '

;device,true=24,retain=2
