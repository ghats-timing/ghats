pro mu_parse_channel_descriptor_se,                          $
                descript,nchan,           $
                ichans,ichane
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
;--------------------------------------------------------------------------
;
; Parsing is based on spawning a Linux pipeline, for simplicity
;
ipos1 = strpos(descript,'C[')
ipos2 = strpos(strmid(descript,ipos1+2),']')

temp = ','+strmid(descript,ipos1+2,ipos2)+','   ; temp contains only relevant channel info, must be parsed

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

;awkstring='echo '+temp+' | tr '','' ''\n'' | tr ''~'' '' '' | awk ''{if(NF == 2) {print $0} else {print $0," ",$0}}'' '

;spawn,awkstring,temp2

;n=n_elements(temp2)

for i=0,nchan-1 do begin
   reads,chanstring(i),dum1,dum2
   ichans(i)=dum1
   ichane(i)=dum2
endfor

;if(nchan ne n) then begin
;   massage,'Inconsistency in channel reading!'
;   retall
;endif

end
