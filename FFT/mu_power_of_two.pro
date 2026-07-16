function mu_power_of_two,j
;
; Function to calculate power of two
;
;-----------------------------------------------------------------------
; Parameters
;
; j                              I: input number
; power_of_two                   O: output power of two
;-----------------------------------------------------------------------
;
a=1l

if(j lt 0l) then begin
   return,0
  endif ELSE begin
   for i=0l,30l do begin
     a=2l^i
     if(a eq j) then begin
       return,1
     endif
   endfor
endelse

return,0

end
