pro mu_get_header_information_czti,filenames,types,tstarts,tends,tress, $
         valtimes1,valtimes2, $
	 nfiles,nmetafiles,nvaltimes, $
         sourcename,obsid,startdate,telescope,instrument,qua
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

	    IF(instrument eq 'CZTI') THEN BEGIN
              type = 'czti'
        ENDIF ELSE BEGIN
            massage,'Unrecognized data type!'
            retall
       ENDELSE

      tstart = fxpar(header,'TSTART')
      tstop  = fxpar(header,'TSTOP')
      timedel= fxpar(header,'TIMEDEL')
      tdim	=fxpar(header,'TDIM*')
      ntdim	=n_elements(tdim)
      tstarts(i,j) = tstart
      tends(i,j)   = tstop
      tress(i,j)   = timedel
;-----------------------------------------------------------
;     now moving to GTI extension
;     for speed reasons, the file is closed and opened again
;;-----------------------------------------------------------
      fxbclose,unit
      errmsg=''
	  
	  ; here we must open the column corresponding to the J quadrant
	  quadrant = string(qua(j),format='(i1)')
	  qext     = 'Q'+quadrant+'_GTI'
      fxbopen,unit,filenames(i,j),qext,header,errmsg=errmsg
      if (errmsg ne '') then begin
	 print,'GTI extension not found!'
	 print,'Using all times'
	 nvaltimes(i,j) = 1
	 valtimes1(0,i,j) = tstarts(i,j)
	 valtimes2(0,i,j) = tends(i,j)
      endif ELSE begin
		  ; GTI extension might be too distant for FXBREAD to be able to access it.
		  ; using MRDFITS instead
		  ; Modified by TMB, 20-01-2017 (IUCAA)
		  tzero=fxpar(header,'TIMEZERO')
		  fxbclose,unit
		  aaa = mrdfits(filenames[i,j],qext)
		  nvaltimes[i,j] = n_elements(aaa)
;	      get number of rows in the GTI list
	      ;nrows=fxpar(header,'NAXIS*')
	      ;nvaltimes(i,j)=nrows(1)
;	      reads in the GTIs
	      ;fxbread,unit,gti,1
	      valtimes1(0:nvaltimes(i,j)-1,i,j) = aaa.(0)+tzero
	      ;fxbread,unit,gti,2
	      valtimes2(0:nvaltimes(i,j)-1,i,j) = aaa.(1)+tzero
      endelse
      ;fxbclose,unit
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
