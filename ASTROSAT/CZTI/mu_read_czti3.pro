pro mu_read_czti3,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,rsour
;
;  Read in CZTI photon files
;  jmeta is the metafile number but ALSO the quadrant number!
;
;-------------------------------------------------------------------------
; Parameters
;
; unit                 I: unit for file reading (already opened)
; tag                  I: file tag identifier
; rdata                O: array with extracted light curve
; np                   I: number of points to accumulate
; channels             I: array of selected/non-selected channels
; selected             I: still a mystery to me
; jmeta                I: index of current metafile
; tres                 I: time resolution of current file
; tres_fft             I: time resolution needed for light curve
; tstart_fft           I: fft start time
; tend_fft             I: fft end time
; irowl                I: current row
;-------------------------------------------------------------------------
common dati,tempi0,tempi1,tempi2,tempi3,energy0,energy1,energy2,energy3,peso0,peso1,peso2,peso3,$
canale_corrente0,canale_corrente1,canale_corrente2,canale_corrente3, $
std21,std22,std23, $
nrows0,nrows1,nrows2,nrows3   ; common block for keeping the data
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)

;
; MU6 addition: at first call read in full file to memory using fxbreadm
if(tag ne canale_corrente3) then begin
		if(baryflag eq 1) then begin
		time_unit = fxbcolnum(unit,'barytime')
	endif
	
	;must find the columns (time and pi) and read in the whole file
	;it should be straightforward
	header      = fxbheader(unit)
	ttype       = fxpar(header,'TTYPE*')
	time_unit   = where(ttype eq 'Time    ')+1
	picol       = where(ttype eq 'PI      ')+1
	wcol        = where(ttype eq 'OPEN_FRACTION')+1
	ppcol       = where(ttype eq 'WEIGHT  ')

	nrows3     =fxpar(header,'NAXIS*')
	nrows3     =nrows3(1)

	if(wcol gt 0) then BEGIN
	    fxbreadm,unit,[time_unit,picol,wcol,ppcol],tempi3,energy3,weights,peso3,BUFFERSIZE=0   ; energy here is PI
	endif else BEGIN
		fxbreadm,unit,[time_unit,picol],tempi3,energy3,BUFFERSIZE=0
		weights = energy3 * 1.0  ; if no weigths are available, set them to 1
		peso3   = energy3 * 1.0  ; if no peso available, set source to 100%
	endelse
	canale_corrente3 = tag
	; Here thorw away all photons with fractions below threshold (and reset nrows)
	;;;weights = weights*0.0+randomu(12,n_elements(weights))   ; TEMPORARY RANDOMIZATION!!!!!!!!!
	buoni   = where(weights ge fraction)
	tempi3   = tempi3[buoni]
	energy3  = energy3[buoni]
	weights = weights[buoni]
	peso3    = peso3[buoni]
	nrows3   = n_elements(buoni)
endif
;
; Obtain info on data columns
;
;fxbtform,header,tbcol,idltype,format,numval,maxval
;
; Find row number of requested starting time
;
ibin   = 0l
irowu  = irowl + 1l
;
; To find the first data point, look for the row whose time covers the
; start time of the first bin of the output array plus half the time
; resolution of the data. This is correct, as any earlier data point
; has its mid-point BEFORE the begin of the output array.
;
time_wanted = tstart_fft
;
exit_flag = 0
;
; Loop to look for it
;

icol = 0

il_mio = where(time_wanted gt tempi3)+1
irowl  = max(il_mio)+1

icol=1   ; data column is #2
;
; Main loop for data reading
;
for irow=irowl,nrows3 do begin
;
; First read in the time info (time is assumed to be 1st column)
;
   ;fxbread,unit,time,1,irow
   time  = tempi3(irow-1) ; MU6

; Here channel energy must be found

ichan = energy3(irow-1)

;Swift   ; Check for > 15k, in which case go to 15k+1
;Swift   if(ichan gt 14999) then ichan = 15000


;
;     find the requested end time
;
      if(time ge tend_fft) then goto,trecento
;
;   The data bin is binned into the output array rdata based on where in
;   rdata its midtime falls. So tstart_fft is the start time of the first bin
;   in data
;
         tbin=time-tstart_fft
         if(tbin gt 0.0d0) then begin
            ibin = long(tbin/tres_fft); + 1l ; must start from 0
           if (ibin le (np-1)) then begin
            rdata(ibin) = rdata(ibin) + channels(ichan)
			rsour(ibin) = rsour(ibin) + channels(ichan)*peso3(ichan)
            std21       = std21 + channels_std21(ichan)  ; Band 1
            std22       = std22 + channels_std22(ichan)  ; Band 2
            std23       = std23 + channels_std23(ichan)  ; Band 3
           endif
         endif
endfor    ; Loop on irowl

trecento:

irowl = min([irow,nrows3])

end
