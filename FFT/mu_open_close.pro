pro mu_open_close,unit,filename,task
;
; Routine to open and close input files
; I hope it works in IDL as it does in FORTRAN
;
;-----------------------------------------------------------------
; Parameters
;
; unit                  I: unit to open or close
; filename              I: name of file to open or close
; task                  I: 'open' or 'close'
;-----------------------------------------------------------------
;
if(task eq 'open') then begin
   errmsg=''
   fxbopen,unit,filename,1,header,errmsg=errmsg
   if(errmsg ne '') then begin
      print,'Error opening file ',filename
      retall
   endif
endif
if (task eq 'close') then begin
   fxbclose,unit
endif

end
