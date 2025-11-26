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
;    Welcome to GHATS 2.0.0      2017 Jan 20
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
;-

common sis, sistema
common vers, versione, data_versione
common versione,version_id
sistema = 'IDL'
versione = 'GHATS  V3.1.0'
data_versione = '2021 Dec 01'
version_id = 'GHATSIB0300     '


;--------------------------------------------------------------------------
!PROMPT = 'Ghats> '
print,' '
print,'-----------------------------------------------------------------'
print,'           Welcome to ',versione,'      ',data_versione
print,'-----------------------------------------------------------------'
print,' '

device,true=24,retain=2
