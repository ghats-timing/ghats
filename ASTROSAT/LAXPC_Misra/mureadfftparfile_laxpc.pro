pro mureadfftparfile_laxpc, $
                     infilename, $
		     outtype,channels,unita,layers,treb,npds,oufilename
;
; Reads parameter file for mufft
;
;--------------------------------------------------------------------
; Parameters
;
;   infilename      	I : input parameter file name
;
;   outtype             O : type of production (power or fft)
;   channels            O : array with start-end channels
;   treb                O : rebin factor
;   npds                O : number of points per pds (if int/long)
;                           length of pds (s) (if float/double)
;   oufilename          O : name for output pds file
;--------------------------------------------------------------------
;
; temporary variables for ascii file reading
;
dum1 = 0
dum2 = 0
;npds = 0l   ; forced to be long
snpds = '' 
stringona = ''
;
openr,1,infilename
   readf,1,infilename
   readf,1,outtype
   readf,1,stringona   ; read a string
   perio = strpos(stringona,'.')    ; finds whether there is a period
   schannels = strsplit(stringona)
   IF(perio eq -1) THEN BEGIN   ; integer values
	   channels = fix(schannels)
      ENDIF ELSE BEGIN   ; real
	   channels = float(schannels)
   ENDELSE
   readf,1,unita
   readf,1,layers
   readf,1,treb
   readf,1,snpds            ; string version
   readf,1,oufilename
close,1
;
; Deal with snpds now
;
perio = strpos(snpds,'.')    ; finds whether there is a period
if(perio eq -1) then begin   ; not there
   npds = 0l                 ; long version
  endif else begin
   npds = 0.0                ; float version
endelse
reads,snpds,npds

end
