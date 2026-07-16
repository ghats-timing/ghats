pro mu_parse_event,event_byte,nbytes,ichan,ioffset,ibit_chan
;
; Extract channel number from PCA event word
;
j    = fix((ioffset -1)/16 + 1)
ioff = ioffset - (j-1)*16
ichan = 0
;
;  Here converts event_byte to the short integer event
;
;event=fix(event_byte,0)
;event=fix([event_byte(1),event_byte(0)],0)  ; need to swap here
; New for GX data 18-jan-2004
event = intarr(2)
event(0)=fix([event_byte(1),event_byte(0)],0)  ; need to swap here
if(j ge 2) then begin
;    event(1)=fix([event_byte(2),0b],0)
    event(1)=fix([0b,event_byte(2)],0)
endif


for ibit=0,ibit_chan-1 do begin
   if(mu_bit_extractor(event(j-1),ibit+ioff) eq 1) then begin
      ichan = ichan+2^(ibit_chan-1-ibit)
   endif
endfor 

ichan = ichan + 1

end
