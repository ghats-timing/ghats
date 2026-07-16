pro mu_check_time,files_units,types,valtimes1,valtimes2, $
               nvaltimes,nfiles,nmetafiles,           $
               tstart_fft,tend_fft,T,i,good,instrument
;
; Procedure to check on timing of current FFT
;
;---------------------------------------------------------------------
; Parameters
;
; files_units            I: array of units of opened input files
; types                  I: array of file types (for ALL files)
; valtimes1              I: array of valid start times for ALL files
; valtimes2              I: array of valid end times for ALL files
; nvaltimes              I: array for numbers of valid intervals for ALL files
; tstart_fft             I: start time for PDS interval
; tend_fft               O: end time for PDS interval
; T                      I: length of a PDS interval
; i                      I: index of current file which is worked on
; good
;---------------------------------------------------------------------	       

; -----------------------------------------------------------------------------
uno:                ; major return label (square one?)
tend_fft = tstart_fft + T
goodtime = 0
tbetter  = 0.0d0
for j=0,nmetafiles-1 do begin
;
;  determine if the FFT fits within the gti's of metafile j
;
   goodtime=0
   for k=0L,nvaltimes(i,j)-1 do begin
      if((tstart_fft ge valtimes1(k,i,j)) and $
         (tend_fft   le valtimes2(k,i,j)) ) then goodtime = 1
   endfor
   if(goodtime eq 0) then begin
;
;     FFT does not fit in metafile j; find new tstart_fft
;
      for k=0L,nvaltimes(i,j)-1 do begin
         if(tstart_fft lt valtimes1(k,i,j)) then begin
;
;           Update tstart_fft to the start of the next GTI after the
;           previous tstart_fft
;
            tstart_fft = valtimes1(k,i,j)
;
;           If metafile j is not an SA file, make sure that the new
;           tstart_fft corresponds to the start of a histogram
;           in an SA file (if there is one)
;
            if((strmid(types(i,j),4,2) ne 'sa') and (instrument eq 'PCA')) then begin
               for jj=0,nmetafiles-1 do begin
                  if(strmid(types(i,jj),4,2) eq 'sa') then begin
		     unita= files_units(jj)
                     mu_check_sa_time,unita,tstart_fft, $
                                      tbetter
                     tstart_fft = tbetter
                     goto,uno
                  endif
               endfor
            endif
;
;           tstart_fft not updated, as it already came from a SA file,
;           or there are no SA files
;
            goto,uno
         endif
      endfor
;
;     If you get out of the loop here, then we're out of GTIs
;
      good = -1
      return
   endif
endfor
;
;If you get our of the loop here, then FFT file fits in gti's of all metafiles
;
good = 0 

end
