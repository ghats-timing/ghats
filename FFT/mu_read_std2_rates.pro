pro mu_read_std2_rates,tmid_fft,std21,std22,std23,vle_times, $
                rate1,rate2,rate3,nvle
;
; Extract current STD2/VLE rate
;
;-------------------------------------------------------------------------
; Parameters
;
; tmid_fft             I: mid   time for the input pds
; std21                O: output band 1 rate
; std22                O: output band 2 rate
; std23                O: output band 3 rate
; vle_times            I: array with times
; rate1                I: array with band 1 light curve
; rate2                I: array with band 2 light curve
; rate3                I: array with band 3 light curve
; nvle                 I: length of input arrays
;-------------------------------------------------------------------------
;
if(tmid_fft lt vle_times(0)) then begin
   std21            = rate1(0)
   std22            = rate2(0)
   std23            = rate3(0)
   return
  endif else begin
   if(tmid_fft gt vle_times(nvle-1)) then begin
      std21            = rate1(nvle-1)
      std22            = rate2(nvle-1)
      std23            = rate3(nvle-1)
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
;centodieci:

;if((i_up-i_low) gt 1) then begin
;   i_mid = (i_up+i_low)/2
;   if(tmid_fft gt vle_times(i_mid)) then begin
;      i_low = i_mid
;   endif else begin
;      i_up = i_mid
;   endelse
;   goto, centodieci
;endif
;-------------------------------------------------
; MU 6.0+
vv = vle_times(where(vle_times gt 0))  ; why this selection?
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
;==========================================================
;
; Computation of overlap  (not used yet)
;
tt1 = vle_times(i_up)-tmid_fft
tt2 = tmid_fft-vle_times(i_low)
;
std21  = float(mean(rate1(i_low:i_up)))
std22  = float(mean(rate2(i_low:i_up)))
std23  = float(mean(rate3(i_low:i_up)))

;======================================================================== 
;  TRIAL NEW VERSION (WITH PROBLEMS)
;vle_resolution = vle_times(1)-vle_times(0)

;i1 = where(vle_times(0:nvle-1) le tstart_fft)
;i1 = i1(n_elements(i1)-1)
;if (i1 eq -1) then begin
;	i1 = 0
;endif
;i2 = where((vle_times(0:nvle-1)+vle_resolution) ge tend_fft)
;i2 = i2(0)
;if (i2 eq -1) then begin
;	i2 = 0
;endif
;print,i1,i2,nvle
;std21 = float(total(rate1(i1:i2)))
;std22 = float(total(rate2(i1:i2)))
;std23 = float(total(rate3(i1:i2)))
;=========================================================================

end
