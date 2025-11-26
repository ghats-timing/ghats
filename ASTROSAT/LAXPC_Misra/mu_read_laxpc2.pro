pro mu_read_laxpc2,unit,laxpc_unit,laxpc_layer,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts
;
;  Read in LAXPC photon files
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
common dati,lay,lax,tempi1,energy1,tempi2,energy2,tempi3,energy3,canale_corrente1,canale_corrente2,canale_corrente3,std21,std22,std23 ; common block for keeping the data
;common dati,lay,lax,tempi1,energy1,canale_corrente1,std21,std22,std23 ; common block for keeping the data
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)

;must find the columns (time and pi) and read in the whole file
;it should be straightforward
header      = fxbheader(unit)
ttype       = fxpar(header,'TTYPE*')
time_unit   = where(ttype eq 'TIME    ')+1
picol       = where(ttype eq 'Channel ')+1
laxpc       = where(ttype eq 'LAXPC_No.')+1
layer       = where(ttype eq 'Layer   ')+1

nrows     =fxpar(header,'NAXIS*')
nrows     =nrows(1)
;
; MU6 addition: at first call read in full file to memory using fxbreadm
if(tag ne canale_corrente1) then begin
	if(baryflag eq 1) then begin
		time_unit = fxbcolnum(unit,'barytime')
	endif
	fxbreadm,unit,[time_unit,picol,laxpc,layer],tempi2,energy2,lax,lay
	canale_corrente1 = tag
endif
;
; Obtain info on data columns
;
fxbtform,header,tbcol,idltype,format,numval,maxval
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

il_mio = where(time_wanted gt tempi2)+1
irowl  = max(il_mio)+1

icol=1   ; data column is #2
;
; Main loop for data reading
;
for irow=irowl,nrows do begin
;
; First read in the time info (time is assumed to be 1st column)
;
   ;fxbread,unit,time,1,irow
   time  = tempi2(irow-1) ; MU6

; Here channel energy must be found


ichan   = energy2(irow-1)
quale_unita = lax(irow-1)
quale_layer = lay(irow-1)

;
;     find the requested end time 
;
      if(time ge tend_fft) then goto,trecento
;
;   The data bin is binned into the output array rdata based on where in
;   rdata its midtime falls. So tstart_fft is the start time of the first bin
;   in data
;
;   Here LAXPC unit selection is made
   IF((laxpc_unit[quale_unita-1] eq 1) and (laxpc_layer[quale_layer-1] eq 1)) THEN BEGIN
         tbin=time-tstart_fft
         if(tbin gt 0.0d0) then begin
            ibin = long(tbin/tres_fft); + 1l ; must start from 0
           if (ibin le (np-1)) then begin
            rdata(ibin) = rdata(ibin) + channels(ichan,0,0,quale_unita-1)  ; GIGUS
            std21       = std21 + channels_std21(ichan)  ; Band 1
            std22       = std22 + channels_std22(ichan)  ; Band 2
            std23       = std23 + channels_std23(ichan)  ; Band 3
           endif
         endif
	ENDIF
endfor    ; Loop on irowl

trecento:

irowl = min([irow,nrows])

end
