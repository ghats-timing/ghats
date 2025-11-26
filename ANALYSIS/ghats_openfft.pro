pro ghats_openfft,filename,unit,dialog=dd
;+
; NAME: 
;      MUXANA_OPENFFT
; PURPOSE: 
;      Open a FFT file  and returns file unit
; EXPLANATION:
;      This procedure opens a FFT file and returns the file unit
;      assigned by the opening process. If the file does not exist
;      it aborts the program, unless the pickup option is specified,
;      in which case it uses dialog_pickup() to prompt for a file
;
; CALLING SEQUENCE: 
;       MUXANA_OPENFFT,FILENAME,UNIT,[/DIALOG]
; INPUTS:
;       FILENAME (I) = name of the input FFT file
;
; OUTPUTS:
;       UNIT     (O) = name of the output unit for opened file
;
; KEYWORDS:
;       DIALOG    K  = optional keyword, If set and file is not found,
;                      opens a dialog window
;
; EXAMPLE:
;       Open a FFT file
;
;       MU> filename='data/gx339.fft'
;       MU> muxana_openpds,filename,unit
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       NONE
; NOTES:
;       
; MODIFICATION HISTORY: 
;       T. Belloni  24 Sep 2002  from muxana_openpds
;		T. Belloni  03 Dec 2010  from Mu6
;-
;--------------------------------------------------------------------------
;
; Open tra file
;
openr,unit,filename,/get_lun,error=err
;;openr,unit,filename,/compress,/get_lun,error=err
if(err ne 0) then begin
  if(keyword_set(dd)) then begin
     filename=dialog_pickfile(filter='*.fft',/read)
     ;;openr,unit,filename,/compress,/get_lun,error=err
     openr,unit,filename,/get_lun,error=err
     if(err ne 0) then begin
        print,'Problems opening FFT file ',filename
        retall
     endif
    endif else begin
     print,'Problems opening FFT file ',filename
     retall
  endelse
endif
;
end
