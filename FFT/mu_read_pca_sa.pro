pro mu_read_pca_sa,unit,rdata,np,channels,selected,jmeta, tres,tres_fft, $
                tstart_fft,tend_fft,nchan,nhisto,irowl,proliferation
;
;  Read in science array files
;
;-------------------------------------------------------------------------
; Parameters
;
; unit                 I: unit for file reading (already opened)
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
; irowl                I: current row
; proliferation        I: variable for sliding window
;-------------------------------------------------------------------------
;
;
; Declaration of byte array
;
header    =fxbheader(unit)
nrows     =fxpar(header,'NAXIS*')
nrows     =nrows(1)
tfields   =strtrim(fxpar(header,'TFIELDS'))
ttype     =strtrim(fxpar(header,'TTYPE*'))
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
time_wanted = tstart_fft + 0.50d0 * tres
;
exit_flag = 0
;
; Loop to look for it
;
icol = 0
while (exit_flag eq 0) do begin
;
   if((irowl eq 0l) or (irowl eq nrows+1)) then begin
      massage,'SA data file does not contain this row'
      retall
   endif
   fxbread,unit,timel,icol+1,irowl
   if(irowl eq nrows) then begin
      timeu = timel + nhisto*tres
     endif ELSE begin
      fxbread,unit,timeu,icol+1,irowu
   endelse

   if((time_wanted ge timel) and (time_wanted le timeu)) then begin
      exit_flag = 1   ; time to leave!
     endif ELSE begin
      if(time_wanted lt timel) then begin
         irowl = irowl - 1l
         irowu = irowl + 1l
         ;print,'Trouble with science array data	'   ; removed 2-2-2012
         ;trtc
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

rcounts=fltarr(nhisto,nchan)
icol=1   ; data column is #2
;
; Main loop for data reading
;
for irow=irowl,nrows do begin
;
; First read in the time info (time is assumed to be 1st column)
;
   fxbread,unit,time,1,irow
;-------------------------------------------------------------------
;  Here the actual reading of data is taking place. It looks like in 
;  any possible input file, the data column to be read is ALWAYS
;  the second. Therefore, we are going to skip the crap and read
;  always that one.
;-------------------------------------------------------------------
;  Integer type (no long integer??)
   if((idltype(icol) eq 2)) then begin
      fxbread,unit,icounts,icol+1,irow
      rcounts = float(icounts)
   endif

;  Float type
   if((idltype(icol) eq 4)) then begin
      fxbread,unit,rr,icol+1,irow
      rcounts = float(rr)
   endif

;  Byte type
   if((idltype(icol) eq 1) and (icol eq 1)) then begin
      fxbread,unit,bcounts,icol+1,irow
      rcounts = float(bcounts)
   endif
;
;  Remap counts to fft data
;
   for i=0,nhisto-1 do begin
;
;     find the requested end time
;
      tempo = time+(float(i+1)-0.5)*tres
      if(tempo ge tend_fft) then goto,trecento
;
;   The data bin is binned into the output array rdata based on where in
;   rdata its midtime falls. So tstart_fft is the start time of the first bin
;   in data
;
      tbin = time + (float(i+1)-0.5)*tres - tstart_fft

      if(tbin gt 0.0d0) then begin
         ibin = long(tbin/tres_fft); + 1l ; must start from 0
         if (ibin le (np-1)) then begin
            for j=0,nchan-1 do begin
;
;  This is where the photons in rcounts go into the array rdata that will be
;  fft'd channels(j) is either one or zero, depending on whether this channel
;  was selected or not
;
               rdata(ibin) = rdata(ibin) + rcounts(i,j)*channels(j)
            endfor
         endif
      endif
   endfor

endfor    ; Loop on irowl

trecento:

irowl = min([irow,nrows])
;
; Take care of sliding window by bringin the index back to allow re-reading the same photons
;

;iiii = irowl
;irowl = irowl - ((np - np/proliferation)/np)/nhisto -1
;print,'IROWL from ',iiii,' to ',irowl

;  Deallocation of memory
bcounts = 0b
icounts = 0
rcounts = 0.0

end
