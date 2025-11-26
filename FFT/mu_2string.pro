pro mu_2string,pot,chdpot
;
; Procedure to produce string with power of two from number (!)
;
;-----------------------------------------------------------------------
; Parameters
;
; pot                   I: input power of two
; chdpot                O: output string "1/pot"
;-----------------------------------------------------------------------
;
if(pot le 0.0) then begin
   massage,'Error in chdpot, negative argument'
   retall
  endif ELSE begin
     if(pot ge 1.0) then begin
        ipot=long(nint(pot))
        if(abs(ipot-pot) gt 1.e-5) then begin
           massage,'Error in chdpot, non-integer argument'
           retall
        endif
        if(mu_power_of_two(ipot) eq 0) then begin
           print,'Warning! number is not a power of two!'
        endif
        chdpot=strtrim(string(ipot),1)
       endif ELSE begin
        ipot=nint(1.0/pot)
;       if(abs(ipot-1.0/pot) gt 1.e-5) then begin  ; HEXTE CHANGE
        if(abs(ipot-1.0/pot) gt 1.e-1) then begin
           massage,'Error in chdpot, non-integer argument'
           retall
        endif
        if(mu_power_of_two(ipot) eq 0) then begin
           print,'Warning! number is not a power of two!'
        endif
        chdpot='1/'+strtrim(string(ipot),1)
     endelse
endelse

end
