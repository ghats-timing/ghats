pro mu_read_pca_gx_buffer,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                tstart_fft,tend_fft,nchan,nhisto,                    $
		ioffset,ibit_chan,irowl
;
;  Read in GoodXenon SE files
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
; nchan                I: number of channels to accumulate
; nhisto               I: current number of histograms
; ioffset              I: current offset for bit reading
; ibit_chan            I: number of bits per channel
; irowl                I: current row
;-------------------------------------------------------------------------
common dati,tempi,eventi,canale_corrente   ; common block for keeping the data
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)
;
; Declaration of byte array for data reading
;
event_byte = bytarr(3)

header    =fxbheader(unit)
nrows     =fxpar(header,'NAXIS*')
nrows     =nrows(1)
tfields   =strtrim(fxpar(header,'TFIELDS'))
ttype     =strtrim(fxpar(header,'TTYPE*'))

time_unit = 1   ; normally the TIME column is the first one and the data are #2
;
; MU6 addition: at first call read in full file to memory using fxbreadm
if(tag ne canale_corrente) then begin
	if(baryflag eq 1) then begin
		time_unit = fxbcolnum(unit,'barytime')
	endif
	fxbreadm,unit,[time_unit,2],tempi,eventi
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
while (exit_flag eq 0) do begin
;
   if((irowl eq 0l) or (irowl eq nrows+1)) then begin
      massage,'We have a problem'
      retall
   endif
   if(irowl eq 1) then begin
      timel = tstart_fft
     endif else begin
      ;fxbread,unit,timel,icol+1,irowl
       timel = tempi(irowl-1) ; MU6
   endelse
   if(irowl eq nrows) then begin
      timeu = tend_fft
     endif ELSE begin
      ;fxbread,unit,timeu,icol+1,irowu
      timeu = tempi(irowu-1) ; MU6
   endelse

   if((time_wanted ge timel) and (time_wanted le timeu)) then begin
      exit_flag = 1   ; time to leave!
     endif ELSE begin
      if(time_wanted lt timel) then begin
         irowl = irowl - 1l
         irowu = irowl + 1l
         ;print,'Trouble with gx data'  ;; MU6 suppressed
;         trtc
         exit_flag = 0   ; try it again
        endif ELSE begin
         if(time_wanted gt timeu) then begin
            irowl = irowl + 1l
            irowu = irowl + 1l
            exit_flag = 0  ; try it again
         endif
      endelse
   endelse
endwhile

icol=1   ; data column is #2
;
; Main loop for data reading
;
for irow=irowl,nrows do begin
;
; First read in the time info (time is assumed to be 1st column)
;
   ;fxbread,unit,time,1,irow
   time = tempi(irow-1) ; MU6
;-------------------------------------------------------------------
;  Here the actual reading of data is taking place. It looks like in 
;  any possible input file, the data column to be read is ALWAYS
;  the second. Therefore, we are going to skip the crap and read
;  always that one.
;-------------------------------------------------------------------
; Type should be BYTE

   if((idltype(icol) eq 1) and (icol eq 1)) then begin
      ;fxbread,unit,event_byte,icol+1,irow
      event_byte = [eventi(0,irow-1),eventi(1,irow-1),eventi(2,irow-1)] ; MU6
;
;     Here the masked filtering must be assigned. It all depends on the
;     little/big endian convention
;
      event=fix(event_byte,0)
      mu_parse_event,event_byte,3,ichan,ioffset,ibit_chan
      ichan=ichan-1   ; start from zero!!
;     if(mu_bit_extractor(event,1) eq 0) then begin   ; powerPC
      if(mu_bit_extractor(event,9) eq 0) then begin   ; Intel
         time_marker = 0
	endif else begin
	 time_marker = 1
      endelse
     endif else begin
      massage,'Unexpected data format in GX data'
   endelse
;
;     find the requested end time
;
      if(time ge tend_fft) then goto,trecento
;
;   The data bin is binned into the output array rdata based on where in
;   rdata its midtime falls. So tstart_fft is the start time of the first bin
;   in data
;
      if(time_marker ne 0) then begin
         tbin=time-tstart_fft
         if(tbin gt 0.0d0) then begin
            ibin = long(tbin/tres_fft); + 1l ; must start from 0
           if (ibin le (np-1)) then begin
            rdata(ibin) = rdata(ibin) + channels(ichan)
           endif
         endif
      endif
endfor    ; Loop on irowl

trecento:

irowl = min([irow,nrows])
;  Deallocation of memory
bcounts = 0b
icounts = 0
rcounts = 0.0

end
