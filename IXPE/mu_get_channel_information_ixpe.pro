pro mu_get_channel_information_ixpe,                      $
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
;
FOR i=0,nmetafiles-1 DO BEGIN
	FOR j=0,nfiles-1 DO BEGIN
           channels1[*,j,i] = indgen(1501)
	   nchannels[j,i] = 1501
    ENDFOR
ENDFOR
channels2 = channels1

end
