pro mu_ep_event_columns, unit, time_col, channel_col, channel_name, bary_col=bary_col
;+
; NAME:
;      MU_EP_EVENT_COLUMNS
; PURPOSE:
;      Locate Einstein Probe event-table columns used by GH_EP.
; EXPLANATION:
;      The EP reader calls this on the opened FITS event extension before
;      reading events, so GH_EP does not rely on a hard-coded column position.
;      TIME is required. PI is preferred for channel selection, with PHA as a
;      fallback. BARYTIME is optional and is returned when present.
;-

header = fxbheader(unit)
ttype  = fxpar(header,'TTYPE*')

time_col    = 0
channel_col = 0
bary_col    = 0
channel_name = ''

for icol=0,n_elements(ttype)-1 do begin
   name = strupcase(strtrim(ttype[icol],2))
   case name of
      'TIME':    time_col = icol+1
      'BARYTIME': bary_col = icol+1
      'PI': begin
         channel_col = icol+1
         channel_name = 'PI'
      end
      'PHA': begin
         if(channel_col eq 0) then begin
            channel_col = icol+1
            channel_name = 'PHA'
         endif
      end
      else:
   endcase
endfor

if(time_col eq 0) then begin
   massage,'EP event extension has no TIME column!'
   retall
endif

if(channel_col eq 0) then begin
   massage,'EP event extension has neither PI nor PHA column!'
   retall
endif

end
