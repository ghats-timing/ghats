pro mu_set_all_channels_xmm,metafiles,filenames,channels,channels1,channels2, $
                    nchannels,nfiles,nmetafiles,canali
;
; Procedure for setting channels to be accumulated
;
;----------------------------------------------------------------------
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
; canali                   I: canali selezionati (start, end)
;----------------------------------------------------------------------
   for j=0,nmetafiles-1 do begin
      for i=0,nfiles-1 do begin
	     channels(canali[0]:canali[1],i,j) = 1
      endfor
   endfor
end
