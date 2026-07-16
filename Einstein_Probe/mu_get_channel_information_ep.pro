pro mu_get_channel_information_ep,filenames,             $
                                     channels1,channels2,nchannels,$
                                     nfiles,nmetafiles
;
; Procedure for obtaining channel information from input files
;
;----------------------------------------------------------------------------
; Parameters
;
; filenames                   I: array of input filenames
; channels1                   O: start channel list for each file
; channels2                   O: end channel list for each file
; nchannels                   O: array of number of channels
; nfiles                      I: input number of files
; nmetafiles                  I: input number of metafiles
;----------------------------------------------------------------------------
unit = 10
for i=0,nmetafiles-1 do begin
   for j=0,nfiles-1 do begin
      errmsg = ''
      fxbopen,unit,filenames(j,i),1,header,errmsg=errmsg
      if(errmsg ne '') then begin
         massage,'Cannot open EP FITS file for channel inspection!'
         retall
      endif

      mu_ep_event_columns,unit,time_col,channel_col,channel_name,bary_col=bary_col

      kmin = 'TLMIN'+strtrim(string(channel_col),2)
      kmax = 'TLMAX'+strtrim(string(channel_col),2)
      cmin = long(fxpar(header,kmin))
      cmax = long(fxpar(header,kmax))

      if(cmax le cmin) then begin
         cmin = 0L
         cmax = 4095L
      endif

      nchan = cmax+1L
      if(nchan gt n_elements(channels1[*,j,i])) then begin
         massage,'EP channel range exceeds MAX_N_CHANNELS in GH_EP!'
         retall
      endif

      channels1[0:nchan-1,j,i] = lindgen(nchan)
      channels2[0:nchan-1,j,i] = channels1[0:nchan-1,j,i]
      nchannels[j,i] = nchan
      fxbclose,unit
   endfor
endfor

end
