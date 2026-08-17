pro mu_get_time_resolution_ep,filenames,nfiles,nmetafiles,tress
;
; Procedure to extract time resolution info
;
;------------------------------------------------------------------------
; Parameters
;
; filenames                 I: array of input filenames
; nfiles                    I: X dimension of filenames
; nmetafiles                I: Y dimension of filenames
; tress                     O: output time resolution
;------------------------------------------------------------------------
type=''
unit = 10
for j=0,nmetafiles-1 do begin
   FOR i=0, nfiles-1 DO BEGIN
      errmsg = ''
      fxbopen, unit, filenames(i,j), 1, header, errmsg=errmsg
      instrument = strtrim(fxpar(header, 'INSTRUME'), 2)
      datamode   = strtrim(fxpar(header, 'DATAMODE'), 2)

      IF (instrument NE 'FXT') THEN BEGIN
         massage, 'Unrecognized instrument: ' + instrument
         retall
      ENDIF

      timedel   =fxpar(header,'TIMEDEL')

      ; Fallback if TIMEDEL is not in the FITS HEADER
      IF (timedel LE 0 OR timedel EQ !VALUES.F_NAN) THEN BEGIN
          CASE datamode OF
              'FF': timedel = 0.050D0
              'PW': timedel = 0.0022D0
              'TM': timedel = 0.00002368D0
              ELSE: timedel = 0.050D0
          ENDCASE
      ENDIF

       tress(i,j) = timedel
       fxbclose,unit	; close FITS file (can give problems?)
   ENDFOR
endfor
end
