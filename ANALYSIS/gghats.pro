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
;    Welcome to GHATS 3.4.0      29 August 2026
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
;       M. Mendez/Codex 26 Jul 2026 update startup version banner to 3.4.0
;-

common sis, sistema
common vers, versione, data_versione
common versione,version_id
sistema = 'GDL'
versione = 'GHATS  V3.4.0'
data_versione = '29 August 2026'
version_id = 'GHATSIB0340     '


;--------------------------------------------------------------------------
!PROMPT = 'Ghats> '
print,' '
print,'-----------------------------------------------------------------'
print,'           Welcome to ',versione,'      ',data_versione
print,'-----------------------------------------------------------------'
print,' '

;device,true=24,retain=2
