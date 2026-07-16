pro mu_channel_selection_ep,metafiles,filenames,channels,channels1,channels2, $
                    nchannels,nfiles,nmetafiles,istart,iend
;
; Procedure for choosing the desired channels for accumlation
;
;------------------------------------------------------------------------
; Parameters
;
; metafiles                I: input array of metafiles
; filenames                I: input array of filenames
; channels                 O: array of selected channels (1 = selected)
; channels1                I: array of arrays of start channels
; channels2                I: array of arrays of end channels
; nchannels                I: how many input channels
; nfiles                   I: number of input files
; nmetafiles               I: number if input metafiles
; istart                   O: final start channel
; iend                     O: final end channel
;------------------------------------------------------------------------
;
;     Author: T. Belloni      Date: 5-FEB-2009 GMU revision
;
common sis,sistema
cpar=''
print,''
print,'Overview of available channels: ',channels1(0,0,0),'-', $
      channels2(nchannels(0,nmetafiles-1)-1,0,nmetafiles-1)
;for j=0,nmetafiles-1 do begin
;   print,''
;   if(metafiles(0) ne 'NONE') then begin
;      print,'Metafile ',metafiles(j)
;     endif ELSE begin
;      print,'File ',filenames(0,j)
;   endelse
;   canali = strarr(249)   ; maximum numbers of RXTE channels (goodXenon)
;   for k=0,nchannels(0,j)-1 do begin
;	if(sistema eq 'IDL') then begin
;      print,channels1(k,0,j),channels2(k,0,j),format='($,i3," -",i3," ")'
;          endif else begin
;      canali(k) = string(channels1(k,0,j),channels2(k,0,j),format='(i3," -",i3," ")')
;	  lista_canali = strjoin(canali,/single)
;	  ;print,channels1(k,0,j),channels2(k,0,j),format='(i3," -",i3," ")'
;	endelse
;   endfor 
;   if(sistema eq 'GDL') then begin
;	  	  print,lista_canali
;	endif
;endfor

andato=0
while (andato eq 0) do begin
   print,''
   print,''
   print,'Choose channels for accumulation'
   print,'ALL: All channels'
   print,'RANGE: Continuous range of channels'
   print,''
   print,format='($,"(ALL) -> ")'
   read,cpar
   if(cpar eq '') then cpar='ALL'
   cpar=strlowcase(strmid(cpar,0,1))

   case cpar of

   'a': begin
;
;  ALL channels
;
           channels = channels*0
           for j=0,nmetafiles-1 do begin
              for i=0,nfiles-1 do begin
                 channels(0:nchannels(i,j)-1,i,j) = 1
              endfor
           endfor
           andato = 1
           istart = channels1(0,0,0)
           iend   = channels2(nchannels(0,nmetafiles-1)-1,0,nmetafiles-1)
        end

   'r': begin
;
;  RANGE of channels
;
           istart=0
           iend=0
           while(istart ge iend) do begin
              print,format='($,"Start channel:  ")'
              read,istart
              print,format='($,"End   channel:  ")'
              read,iend
           endwhile
           for ichan=istart,iend do begin
              for j=0,nmetafiles-1 do begin
                 for i=0,nfiles-1 do begin
                    for k=0,nchannels(i,j)-1 do begin
                       if((ichan ge channels1(k,i,j)) and $
                          (ichan le channels2(k,i,j))) then $
                             channels(k,i,j) = 1
                    endfor
                 endfor
              endfor
           endfor
           andato = 1
        end

   else: begin
            print,'Unrecognized case, try again'
         end

   endcase
endwhile

end
