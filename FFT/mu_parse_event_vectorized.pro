pro mu_parse_event_vectorized,event_byte,nbytes,ichan,ioffset,ibit_chan;
;************* STILL UNDER DEVELOPMENT!!!
;
;      Author: T. Belloni      Date: 17-APR-2009  vectorized version
;
; Extract channel number from PCA event word
; Input is event_byte, the event to parse, ioffset/ibit_chan: single numbers
;
j    = fix((ioffset -1)/16 + 1)
ioff = ioffset - (j-1)*16

n = size(event_byte)
n = n(2)
ichan = intarr(n)
;
;  Here converts event_byte to the short integer event
;
event = intarr(n,2)
event(*,0)=fix([event_byte(1),event_byte(0)],0,1,n)  ; need to swap here
if(j ge 2) then begin
    event(*,1)=fix([0b,event_byte(2)],0,1,n)
endif

; How do I do it here?
ibit  = indgen(ibit_chan)
bitto = indgen(n,ibit_chan)

bitto(*,????) = mu_bit_extractor(event(*,j-1),ibit+ioff)

total(2^(ibit_chan-1-ibit)*bitto)


for ibit=0,ibit_chan-1 do begin
   if(mu_bit_extractor(event(j-1),ibit+ioff) eq 1) then begin
      ichan = ichan+2^(ibit_chan-1-ibit)
   endif
endfor 

ichan = ichan + 1

end
