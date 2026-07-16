pro mu_read_ep,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy
;
;  Read in Einstein Probe photon event files.
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
common dati,tempi,eventi,canale_corrente,std21,std22,std23   ; common block for keeping the data
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)

; Find the columns and read the whole file into memory on first use.
header = fxbheader(unit)
mu_ep_event_columns,unit,time_unit,picol,channel_name,bary_col=bary_col

nrows     =fxpar(header,'NAXIS*')
nrows     =nrows(1)
;
; MU6 addition: at first call read in full file to memory using fxbreadm
if(tag ne canale_corrente) then begin
	if(baryflag eq 1) then begin
      if(bary_col eq 0) then begin
         massage,'/BARY requested, but EP event extension has no BARYTIME column!'
         retall
      endif
		time_unit = bary_col
	endif
   	fxbreadm,unit,[time_unit,picol],tempi,energy
	canale_corrente = tag
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

il_mio = where(time_wanted gt tempi)+1
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
   time  = tempi(irow-1) ; MU6
   ichan = energy(irow-1)
   if(ichan lt 0) then goto,duecento
   if(ichan gt n_elements(channels)-1) then goto,duecento
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
            std21       = std21 + channels_std21(ichan)  ; Band 1
            std22       = std22 + channels_std22(ichan)  ; Band 2
            std23       = std23 + channels_std23(ichan)  ; Band 3
           endif
         endif
duecento:
endfor    ; Loop on irowl

trecento:

irowl = min([irow,nrows])
;print,nffts,tstart_fft-543896000.0d0

end
