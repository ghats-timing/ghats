pro ghats_openpds,filename,unit,dialog=dd,help=help
;+
; NAME: 
;      GHATS_OPENPDS
; PURPOSE: 
;      Open a PDS file  and returns file unit
; EXPLANATION:
;      This procedure opens a PDS file and returns the file unit
;      assigned by the opening process. If the file does not exist
;      it aborts the program, unless the pickup option is specified,
;      in which case it uses dialog_pickup() to prompt for a file
;
; CALLING SEQUENCE: 
;       GHATS_OPENPDS,FILENAME,UNIT,[/DIALOG]
; INPUTS:
;       FILENAME (I) = name of the input PDS file
;
; OUTPUTS:
;       UNIT     (O) = name of the output unit for opened file
;
; KEYWORDS:
;       DIALOG    K  = optional keyword, If set and file is not found,
;                      opens a dialog window
;       HELP         = Print help and return
;
; EXAMPLE:
;       Open a PDS file
;
;       MU> filename='data/gx339.pds'
;       MU> ghats_openpds,filename,unit
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       NONE
; NOTES:
;       This is a helper routine normally called by higher-level PDS readers.
; MODIFICATION HISTORY: 
;       T. Belloni   9 Nov 2001  implementation
;       T. Belloni  12 Mag 2002  adapted for MUFFT format
;       T. Belloni  13 Apr 2003  check for compression
;       T. Belloni  14 May 2005  endian independence (for mac)
;	T. Belloni  08 Dec 2008  MU5: no more endian dychotomy: only Intel. Check for compression fixed
;       T. Belloni  13 Feb 2009  free_lun
;       T. Belloni  11 Nov 2009  from muxana_openpds to ghats_openpds
;		T. Belloni  11 Oct 2018  added filename existence check for GDL
;       2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,'GHATS_OPENPDS'
   print,'Helper routine: normally called by GHATS PDS readers, not run standalone.'
   print,'Calling sequence: ghats_openpds, filename, unit [, /dialog]'
   print,'Opens a GHATS .pds file and returns an IDL logical unit.'
   return
endif
;
;--------------------------------------------------------------------------
;
; Open pds file
;
IF(not(isa(filename))) THEN filename=''

openr,unit,filename,/get_lun,error=err
if(err ne 0) then begin
  if(keyword_set(dd)) then begin
     filename=dialog_pickfile(filter='*.pds',/read)
     openr,unit,filename,/get_lun,error=err
     if(err ne 0) then begin
        print,'Problems opening PDS file ',filename
        retall
     endif
    endif else begin
     print,'Problems opening PDS file ',filename
     retall
  endelse
endif
end
