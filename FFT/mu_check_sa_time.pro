pro mu_check_sa_time,unit,tgood,tbetter
;
; Procedure to check whether there is a SA file with a compatible 
; start time. Called by CHECK_TIME
;
;-------------------------------------------------------------------------
; Parameters
;
; unit                   I: unit number for file
; tgood                  I: start time to check
; tbetter                O: output new possible value for tstart
;-------------------------------------------------------------------------
type=''
header    = fxbheader(unit)
extname	  = strtrim(fxpar(header,'EXTNAME'))
nrows     = fxpar(header,'NAXIS*')
nrows     = nrows(1)

if(extname eq 'XTE_SA') then begin
   type = 'sa'
  endif ELSE begin
   massage,'Only know how to handle XTE_SA files'
endelse
;
; Determine time resolution of the observation
;
timedel   = double(0.0)
timedel	  = strtrim(fxpar(header,'TIMEDEL'))
tfields   = strtrim(fxpar(header,'TFIELDS'))
;
; Obtain info on data columns
;
;for icol=0,tfields-1 do begin
   fxbtform,header,tbcol,idltype,formato,numval,maxval
;endfor
;
; declaration of some double precision variables
;
timel = 0.0d0
timeu = timel
time  = timel
;
; Loop through the columns
;
for icol=0,tfields-1 do begin
   if(idltype(icol) eq 5) then begin
;
; Double precision values: it means that it's a time column. What a waste!
;
      irowl = 1
      irowu = nrows
;
; Read in first and last time values
;
      fxbread,unit,timel,icol+1,irowl
      fxbread,unit,timeu,icol+1,irowu
;
; Add time resolution to last time value
;
      timeu = timeu + timedel
;
; Find row number of requested starting time
;
      time_wanted = tgood
   
      if((time_wanted lt timel) or (time_wanted gt timeu)) then begin
         massage,'Requested time out of range'
         retall
      endif

      centodieci:
      if((irowu-irowl) gt 1 ) then begin
         irowm = (irowu+irowl)/2
         fxbread,unit,time,icol+1,irowm
         if(time_wanted gt time) then begin
            irowl = irowm
           endif ELSE begin
            irowu = irowm
         endelse
         goto,centodieci
      endif

      irow_tbetter = irowu
      fxbread,unit,tbetter,icol+1,irow_tbetter
      goto,seisei
   endif
endfor

seisei:

end
