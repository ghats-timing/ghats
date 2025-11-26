pro mu_channel_selection_laxpc,metafiles,filenames,channels,channels1,channels2, $
                    nchannels,nfiles,nmetafiles,istart,iend,multiple,resp_flag,respfilelist, $
						energies_from_event,energie
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
print,'Overview of available channels: all channels 0-1023'
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

multiple = 0    ; GIGUS
andato=0
while (andato eq 0) do begin
   print,''
   print,''
   print,'Choose channels for accumulation'
   print,'ALL: All channels'
   print,'RANGE: Continuous range of channels'
   print,'UNIT: Separate range for LAXPC units'
   ;IF(resp_flag eq 1) THEN BEGIN
	   print,'ENERGY: energy range in keV. If negative, channels will be computed from DRM. If positive from event file itself'
   ;ENDIF
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
           channels= channels*0+1
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
		'u': BEGIN
			FOR IU=0,2 DO BEGIN
				multiple = 1     ; GIGUS
	            istart=0
	            iend=0
	            while(istart ge iend) do begin
					print,'LAXPC Unit',(IU+1)*10
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
	                              channels(k,i,j,IU) = 1
	                     endfor
	                  endfor
	               endfor
	            endfor
			ENDFOR
			andato = 1
		END
		'e': BEGIN
            istart=0
            iend=0
            while(abs(istart) ge abs(iend)) do begin
               print,format='($,"Start energy (keV):  ")'
               read,istart
               print,format='($,"End   energy (keV):  ")'
               read,iend
            endwhile
			IF(istart le 0) THEN BEGIN
				istart = -istart   ; bring them back to positive
				iend   = -iend
				IF(resp_flag eq 0) THEN BEGIN
					print,'Option ENERGY not valid if DRM files are not specified with keyword RESP'
					stop
				ENDIF ELSE BEGIN
					; Read channel ranges from DRM files
					OPENR,lun,respfilelist,/GET_LUN
					   rmf = STRARR(3)
					   READF,lun,rmf
					FREE_LUN,lun
				ENDELSE
				FOR irmf=0,2 DO BEGIN
					fxbopen,unitr,rmf[irmf],1,hea,errmsg=errmsg
					   ttype       = fxpar(hea,'TTYPE*')
					   cha         = where(ttype eq 'CHANNEL ')+1
					   e1          = where(ttype eq 'E_MIN   ')+1
					   e2          = where(ttype eq 'E_MAX   ')+1
					   fxbreadm,unitr,[cha,e1,e2],can,emin,emax
					fxbclose,unitr
					; interpolation for channels 0-1023
					IF(irmf eq 1) THEN BEGIN
						nreb = 4
					ENDIF ELSE BEGIN
						nreb = 2
					ENDELSE
					nc       = n_elements(emin)
					canale   = nreb*findgen(nc)
					canalone = findgen(nreb*nc)
					emin_new = INTERPOL(emin, canale, canalone)
					emax_new = INTERPOL(emax, canale, canalone)
					c1 = min(where(emin_new ge istart))
					c2 = max(where(emax_new le iend))
					;print,'*** ',c1,c2
					for ichan=c1,c2 do begin
				       for j=0,nmetafiles-1 do begin
						  for i=0,nfiles-1 do begin
						     for k=0,nchannels(i,j)-1 do begin
						        if((ichan ge channels1(k,i,j)) and $
						           (ichan le channels2(k,i,j))) then $
						             channels(k,i,j,irmf) = 1
						     endfor
						  endfor
						endfor
					endfor	
				ENDFOR
			ENDIF ELSE BEGIN
				energies_from_event = 1
				energie = [istart,iend]
			ENDELSE
			andato = 1
		END

   else: begin
            print,'Unrecognized case, try again'
         end

   endcase
endwhile

end
