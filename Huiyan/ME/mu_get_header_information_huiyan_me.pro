pro mu_get_header_information_huiyan_me,filenames,types,tstarts,tends,tress, $
         valtimes1,valtimes2, $
	 nfiles,nmetafiles,nvaltimes, $
         sourcename,obsid,startdate,telescope,instrument
;
; Procedure to obtain necessary information on input files (all of them)
;
;----------------------------------------------------------------------
; Parameters
;
; filenames                     I: input list of filenames
; types                         O: list of data types
; tstarts                       O: list of start times
; tends                         O: list of end times
; tress                         O: list of time resolutions
; valtimes1                     O: array of valid start times
; valtimes2                     O: array of valid end times
; nvaltimes                     O: dimension of valtimes1/2
; sourcename                    O: string with source name
; obsid                         O: observation ID
; startdate                     O: start data of observation
; telescope                     O: satellite name
; instrument                    O: instrument name
;----------------------------------------------------------------------
unit = 10
errmsg=''
for j=0,nmetafiles-1 do begin
   for i=0,nfiles-1 do begin

      fxbopen,unit,filenames(i,j),1,header,errmsg=errmsg
;
;  get additional reference information
;
      if((i+j) eq 0) then begin
         sourcename= fxpar(header,'OBJECT')
         obsid     = fxpar(header,'OBS_ID')
         startdate = fxpar(header,'DATE-OBS')
         telescope = fxpar(header,'TELESCOP')
		 ;telescope = '慧眼'
      endif
; If the date is in old format, cosntruct it in the new one
      if(strpos(startdate,'/') ge 0) then begin
	     starttime = fxpar(header,'TIME-OBS')
	     giorno    = strmid(startdate,0,2)
	     mese      = strmid(startdate,3,2)
	     anno      = strmid(startdate,6,2)
	     ora       = strmid(starttime,0,2)
	     minuto    = strmid(starttime,3,2)
	     secondo   = strmid(starttime,6,2)
	     startdate = '19'+anno+'-'+mese+'-'+giorno+'T'+ $
	                  ora+':'+minuto+':'+secondo
	  endif

      instrument= strtrim(fxpar(header,'INSTRUME'))
      datamode	= strtrim(fxpar(header,'DATAMODE'))
      extname	  = strtrim(fxpar(header,'EXTNAME'))

	   if(instrument eq 'ME') then begin
		types(i,j) = 'Huiyan_ME'
	   endif
	   if(instrument ne 'ME') then begin
		massage,'Unrecognized data type!'
		retall
	endif

      tstart = fxpar(header,'TSTART')
      tstop  = fxpar(header,'TSTOP')
      timedel= fxpar(header,'TIMEDEL')
	  ;timedel = 6.0e-6 ; ******** HARDCODED
      tdim	=fxpar(header,'TDIM*')
      ntdim	=n_elements(tdim)
      tstarts(i,j) = tstart
      tends(i,j)   = tstop
      tress(i,j)   = timedel
;
;     now moving to GTI extension
;     for speed reasons, the file is closed and opened again
;
      fxbclose,unit
      errmsg=''
      fxbopen,unit,filenames(i,j),2,header,errmsg=errmsg
;      status = fxmove(unit,2,/Silent)  ; too slow......_@_Y

      if (errmsg ne '') then begin
	 print,'GTI extension not found!'
	 print,'Using all times'
	 nvaltimes(i,j) = 1
	 valtimes1(0,i,j) = tstarts(i,j)
	 valtimes2(0,i,j) = tends(i,j)
      endif ELSE begin
;	    get number of rows in the GTI list
	 nrows=fxpar(header,'NAXIS*')
	 nvaltimes(i,j)=nrows(1)
	 tzero=fxpar(header,'TIMEZERO')
;	    reads in the GTIs
	 fxbread,unit,gti,1
	 gti = gti+tzero
	 valtimes1(0:nvaltimes(i,j)-1,i,j) = gti
	 fxbread,unit,gti,2
	 gti = gti+tzero
	 valtimes2(0:nvaltimes(i,j)-1,i,j) = gti
      endelse
      fxbclose,unit
   endfor
endfor
;
; trim arrays
;
   mass = max(nvaltimes)
   valtimes1 = valtimes1(0:mass-1,*,*)
   valtimes2 = valtimes2(0:mass-1,*,*)
;
;  consistency and sorting check
;
for j=0,nmetafiles-1 do begin
   for i=1,nfiles-1 do begin

      if(tstarts(i,j) lt tstarts(i-1,j)) then begin
	 massage,'Input files are not sorted in time!'
	 retall
	 
      endif

      if(types(i,j) ne types(i-1,j)) then begin
	 massage,'Inconsistent files within the same metafile!'
	 retall
      endif
   endfor
endfor

end
