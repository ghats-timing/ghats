pro mu_parse_channel_descriptor_sa,                                 $
             descript,nchan,ichans,           $
             ichane,ioffset,ibit_chan
;
; Procedure for the parsing of the FITSIO channel descriptor (science SA data)
;
;-------------------------------------------------------------------------
; Parameters
;
; descript                   I: input descriptor string
; nchan                      I: supposed number of bins

; ichans                     O: starting channels of energy bins
; ichane                     O: ending channels of energy bins
; ioffset                    O: ioffset for event data
; ibit_chan                  O: bit channel info for event data
;--------------------------------------------------------------------------
;
; Parsing was based on spawning a Linux pipeline, for simplicity, but
; this was changed in Mu 3.0 and a full internal portable parsing was implemented
;
ipos1 = strpos(descript,'C[')

ipos2 = strpos(strmid(descript,ipos1+2),']')

temp = ','+strmid(descript,ipos1+2,ipos2)+',' ; temp contains only relevant channel info, must be parsed

goodx = strpos(temp,':')
if (goodx lt 0) then begin
;   data is NOT goodxenon

; event data do not get number of channels from input.
; They are computed here from the number of commas
; in the string.
nchan = 0
for i=1,strlen(temp)-1 do begin
   su = strmid(temp,i,1)
   if((strcmp(su,',') eq 1)) then begin
      nchan = nchan+1
   endif
endfor

;
; Loop to replace tildes
;
for i=0,strlen(temp)-1 do begin
   su = strmid(temp,i,1)
   if((strcmp(su,'~') eq 1)) then begin
      strput,temp,' ',i
   endif
endfor

chanstring = strarr(nchan)
su         = ' '
j          = 0
for i=0,strlen(temp)-2 do begin
   su = strmid(temp,i,1)            ; find comma
   if((strcmp(su,',') eq 1)) then begin
      prox = strpos(temp,',',i+1)   ; find next comma
      temp2 = strmid(temp,i+1,prox-i-1)
      ; Check whether there is a space. If not, it's a single channel and must be duplicated
      if(strpos(temp2,' ') eq -1) then begin
         temp2 = temp2+' '+temp2
      endif
      chanstring(j) = temp2
      j=j+1
   endif
endfor


;   awkstring='echo '+temp+' | tr '','' ''\n'' | tr ''~'' '' '' | awk ''{if(NF == 2) {print $0} else {print $0," ",$0}}'' '

;   spawn,awkstring,temp2

;   n=n_elements(temp2)

   for i=0,nchan-1 do begin
      reads,chanstring(i),dum1,dum2
      ichans(i)=dum1
      ichane(i)=dum2
   endfor

   reads,strmid(descript,ipos2+ipos1+4,1),ibit_chan
   nchan = 2^ibit_chan
;   if(nchan ne n) then begin
;      massage,'Inconsistency in channel reading!'
;      retall
;   endif
  endif ELSE begin
;  data IS goodxenon
   nchan = 256
   ichans=indgen(256)
   ichane=ichans
   ibit_chan=8
endelse

;
; now get ioffset
;
ioffset = 1
iposl   = ipos1

iposr=strpos(strmid(descript,0,iposl-1),'}',/REVERSE_SEARCH)
iposl=strpos(strmid(descript,0,iposl-1),'{',/REVERSE_SEARCH)
i=0

while (iposl ge 0) do begin
   reads,strmid(descript,iposl+1,iposl+iposr-1),i
   ioffset = ioffset + i
   iposr=strpos(strmid(descript,0,iposl-1),'}',/REVERSE_SEARCH)
   iposl=strpos(strmid(descript,0,iposl-1),'{',/REVERSE_SEARCH)
endwhile

end
