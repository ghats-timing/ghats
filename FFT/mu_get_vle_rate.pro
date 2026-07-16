pro mu_get_vle_rate,tmid_fft,current_vle_rate,ndet,vle_times, $
                vle_counts,npcus,nvle
;
; Extract current VLE rate
;
;-------------------------------------------------------------------------
; Parameters
;
; tmid_fft             I: middle time for the input pds
; current_vle_rate     O: output VLE rate
; ndet                 O: output number of detectors
; vle_times            I: array with VLE times
; vle_counts           I: array with VLE rates
; npcus                I: array with number of PCUs
; nvle                 I: length of input arrays
;-------------------------------------------------------------------------
;
if(tmid_fft lt vle_times(0)) then begin
   ndet             = npcus(0)
   current_vle_rate = vle_counts(0)
   return
  endif else begin
   if(tmid_fft gt vle_times(nvle-1)) then begin
      ndet             = npcus(nvle-1)
      current_vle_rate = vle_counts(nvle-1)
      return
   endif
endelse
;
i_low = 0
i_up  = nvle-1
;
; Find row number of requested starting time
;
; OLD MU <= 6.0 WAY -------------------------------
;while((i_up-i_low) gt 1) do begin
;   i_mid = (i_up+i_low)/2
;   if(tmid_fft gt vle_times(i_mid)) then begin
;      i_low = i_mid
;   endif else begin
;      i_up = i_mid
;   endelse
;endwhile
;-------------------------------------------------
vv = vle_times(where(vle_times gt 0))
i_low = max(where(vv lt tmid_fft))
i_up  = min(where(vv ge tmid_fft))
if(i_low eq -1) then begin
    i_low = 0
    i_up  = 1
endif
if(i_up eq -1)  then begin
    i_up  = nvle-1
    i_low = i_up-1
endif
;
ndet = fix(float(npcus(i_up) - npcus(i_low)) /               $
       (vle_times(i_up) - vle_times(i_low))  *               $
       (tmid_fft - vle_times(i_low)) + npcus(i_low) )
       
current_vle_rate = (vle_counts(i_up) - vle_counts(i_low)) /  $
                   (vle_times(i_up) - vle_times(i_low))   *  $
		   (tmid_fft - vle_times(i_low))          +  $
		   vle_counts(i_low)
		   
end
