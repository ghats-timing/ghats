pro gh_huiyan,infilename,outtype,canali,quali_grade,treb,npds,oufilename,ghx=ghx,bands=bands,gti=usergti,sliding=sliding,bary=bary, $
	       wind=wind,wpar=wpar
;+
; NAME:
;      GH_HUIYAN
; PURPOSE:
;      Production of Huiyan/LE-ME-HE PDS and FFT files for the GHATS project
; EXPLANATION:
;      This procedure, coming from MU, produces
;      PDS and FFT files from Huiyan/LE-ME-HE data
;
; CALLING SEQUENCE:
;       GH_HUIYAN
;       GH_HUIYAN,parfilename
;       GH_HUIYAN,infilename,outtype,channels,treb,npds,outfilename[,/muxana][,/turbo][,bands=B][,gti=GTI][,sliding=SLIDING][,/BARY]
; INPUTS:
;       PARFILENAME = name of parameter file containing the six
;                     parameters of the full command line version
;                     of this procedure
;
;       INFILENAME  = name of input file (or meta file or meta-meta file)
;       OUTTYPE     = POWER for power density spectra, FFT for full complex FFT
;       CHANNELS    = array with start,end of PHA channel range
;                     (example: [0,35])
;       TREB        = rebin factor
;       NPDS        = if int/long: number of light-curve points to be
;                                  FFTed in each interval
;                   = if float/doule: length of interval in seconds
;       OUTFILENAME = name of output .pds or .fft file
;		
; KEYWORDS:
;		GHX         = if set, GHATS_ALL is launched at the end of PDS production to generate
;					  a PDF with the output information. The PDS filename is 'gh.pds'
;		BANDS       = channel bands for the production of three standard2 light curves. It is an array of six
;			  		  elements ([0,10,11,20,21,40]) means three curves in the channel ranges 0-10, 11-20 and
;			 		  21-40 respectively. They can then be used to produce X-ray hardnesses. Default values
;					  are [3,6,7,14,15,20]
;		GTI         = optional user-supplied GTIs, to be intersected with the ones in the data files. 
;			          If set, it must be a valid filename. Two possible input formats are possible. Either
;					  an ASCII file with lines containing start end time of each GTI,
;					  or a standard GTI FITS file
;		SLIDING     = optional parameter for sliding the time window. If not set, the time intervals for 
;					  the FFT (data streth) do not overlap. If set to N, they overlap by an Nth of their
;					  total duration. As an example, if the data intervals are 16 seconds long, SLIDING=4
;					  will overlap two consecutive intervals for four seconds. If this keyword is
;					  set, the averaged PDS obtained with GHX does not obviously make statistical sense.
;		BARY        = flag for barycentered data. If set, the program will not read the TIME column,
;					  but BARYTIME. It is the user's responsibility to make sure that the BARYTIME 
;					  column exists.
;
;
; OUTPUTS:
;       NONE
;
; EXAMPLE:
;       GH_HUIYAN
;         Prompts user interactively for parameters
;       GH_HUIYAN,'#mysource.par'
;         Takes parameters from the specified parameter file
;       GH_HUIYAN,'@binned.lis','POWER',[0,35],1,4096,'mysource.pds'
;         Takes parameters from the command line
;
; COMMON BLOCKS:
;       dati         : common block to keep the actual data
;	    sis          : variable containing the system name (IDL/GDL)
;		barycentered : common block with the barycenter flag
;		versione     : common block with the software's version number
; ROUTINES USED:
;
; NOTES
;       None
; MODIFICATION HISTORY:
;		T. Belloni  28 Mar 2018  from GH_HUIYAN_LE (prototype)
;		T. Belloni  28 Feb 2019  general version
;		T. Belloni  09 Jul 2019  fixed color accumulation
;-
;-------------------------------------------------------------
common sis, sistema   ; common block with system variable
common barycentered,baryflag
common versione,version_id
common finestre,finestra,winn

; ================== VERSION ID ===============================
;version_id = 'GHATSIB0002     '
if(sistema eq 'GDL') then begin
	strput,version_id,'G',5
endif
LI = (byte(1,0,1))[0]
if(LI) then begin
	strput,version_id,'L',6
endif
;==============================================================
; Check whether windowing is requested, if not set default window
;
if(~keyword_set(wind)) then begin
	wind = 'Boxcar'
endif
;
; Check whether the bary flag is set
baryflag  = 0
if(keyword_set(bary)) then begin
	baryflag = 1
endif

common dati,tempi,energy,canale_corrente,std21,std22,std23,grade ; common block for keeping the data
canale_corrente     = -1 ; for fxbreadm call
std21            = 0.0
std22            = 0.0
std23            = 0.0
ndet             = 1
current_vle_rate = 0.0

;================================================================
if(keyword_set(sliding)) then begin
	proliferation = sliding   ; to estimate correct number of FFT
endif else begin
	proliferation = 1
endelse
if(keyword_set(usergti)) then begin
	gti_flag = 1
endif else begin
	gti_flag = 0
endelse
;----------
MAX_N_FILES =  200
MAX_N_METAFILES =  6
MAX_N_VALTIMES  =  1000
MAX_N_CHANNELS  =  1537
MAX_N_COL       =  256
MAX_VLE         = 6000    ; number for rates

turbo = 0
;================================================================
if(keyword_set(bands)) then begin
	BANDE           = bands
endif else begin
    BANDE           = [20,100,101,400,401,1023]
endelse

; Prepare arrays for accepted channels for hardness
channels_std21 = intarr(MAX_N_CHANNELS)*0
channels_std21[BANDE[0]:BANDE[1]] = 1
channels_std22 = intarr(MAX_N_CHANNELS)*0
channels_std22[BANDE[2]:BANDE[3]] = 1
channels_std23 = intarr(MAX_N_CHANNELS)*0
channels_std23[BANDE[4]:BANDE[5]] = 1

if(keyword_set(usergti)) then begin
	mu_read_user_gti_file,usergti,u_gti1,u_gti2 ; read in user-supplied GTI file
endif

metafiles = strarr(MAX_N_METAFILES)
filenames = strarr(MAX_N_FILES,MAX_N_METAFILES)
;
new_output_flag = 1   ; flag to open the output flag
output_unit     = 9   ; unit for output file writing
;
;
; Three fashions for calling this program:
;   0 pars: user will be prompted
;   1 pars: parameter file with input parameters
;   6 pars: all parameters are specified
case n_params() of

   0: begin
     	parameters_known = 0
     	infilename       = ''
     	outtype          = ''
     	canali           = intarr(2)
     	treb             = 0
     	npds             = 0l   ; irrelevant statement
     	oufilename       = ''
      end

   1: begin
;
; input parameter is a parameter file name
;
     	parameters_known = 1
        outtype          = ''
     	canali           = intarr(2)
     	treb            = 0
;       npds           = 0l
     	oufilename       = ''
     	mureadfftparfile,infilename,outtype,canali,treb,npds,oufilename
      end

   6: begin
     	parameters_known = 1
      end

   else: begin
     	print,'Usage: gh_huiyan_le,infile,outtype,channels,treb,npds,oufile'
     	stop
      end
endcase
;-------------------------------------------------------------------
;
;   get input files through a GUI window for IDL and a text list for GDL
;
if (parameters_known eq 0) then begin
   mu_ask,infilename,'*.fits;@*'
endif
;-------------------------------------------------------------------
;
; obtain input filenames
;
mu_fetch_filenames,infilename,MAX_N_FILES,MAX_N_METAFILES, $
          filenames,metafiles,filenameroot,nfiles,nmetafiles
;
filenames = filenames(0:nfiles-1,0:nmetafiles-1)  ; crop array
metafiles = metafiles(0:nmetafiles-1)       ; crop array
;-------------------------------------------------------------------
;
; obtain start times
;
tstarts   = dblarr(nfiles,nmetafiles)
mu_get_start_times,filenames,nfiles,nmetafiles,tstarts
;-------------------------------------------------------------------
;
; sort array of start times
;
mu_sort_start_times,filenames,tstarts,nfiles,nmetafiles
;-------------------------------------------------------------------
;
; get time resolutions of different files
;
tress     = dblarr(nfiles,nmetafiles)
; MM mu_get_time_resolution_huiyan_le,filenames,nfiles,nmetafiles,tress
; MM
; mu_get_time_resolution_huiyan_le does not exist
mu_get_time_resolution_huiyan,filenames,nfiles,nmetafiles,tress
; MM
;-------------------------------------------------------------------
;
; order files by time resolution
;
mu_sort_time_resolution,filenames,metafiles,tress,nfiles,nmetafiles
;-------------------------------------------------------------------
;
; obtain relevant info from the input files
;
types     = strarr(nfiles,nmetafiles)
tends     = dblarr(nfiles,nmetafiles)
valtimes1 = dblarr(MAX_N_VALTIMES,nfiles,nmetafiles)
valtimes2 = valtimes1
nvaltimes = intarr(nfiles,nmetafiles)
mu_get_header_information_huiyan_le,filenames,types,tstarts,tends,tress,valtimes1,valtimes2, $
     nfiles,nmetafiles,nvaltimes, $
         sourcename,obsid,startdate,telescope,instrument

;-------------------------------------------------------------------
;
; obtain channel information
;
channels1 = intarr(MAX_N_CHANNELS,nfiles,nmetafiles)
channels2 = channels1
nchannels = intarr(nfiles,nmetafiles)
mu_get_channel_information_huiyan_le,channels1,channels2,nchannels,nfiles,nmetafiles
;-------------------------------------------------------------------
channels = intarr(MAX_N_CHANNELS,nfiles,nmetafiles)*0
if (parameters_known eq 0) then begin
;
; select channels interactively
;
   mu_channel_selection_huiyan_le,metafiles,filenames,channels,channels1,channels2, $
                   nchannels,nfiles,nmetafiles,istart,iend
   canali(0) = istart
   canali(1) = iend
endif ELSE begin
;
; simply use the available range
;
   mu_set_all_channels,metafiles,filenames,channels,channels1,channels2, $
                   nchannels,nfiles,nmetafiles,canali
endelse

canali(0) = min(channels1(where(channels eq 1)))
canali(1) = max(channels2(where(channels eq 1)))
;-------------------------------------------------------------------
;
; choose type of output
;
chout=''
andato=0
while (andato eq 0) do begin
   if (parameters_known eq 0) then begin
      print,''
      print,''
      print,'Choose type of accumulation'
      print,'POWER: Power spectra (PDS file)'
      print,'FFT:   Fourier spectra (FFT file)'
      print,''
      print,format='($,"(POWER) -> ")'
      read,chout
      if(chout eq '') then chout='POWER'
     endif ELSE begin
      chout=outtype
   endelse
   chou2=strlowcase(strmid(chout,0,1))

   case chou2 of

   'p': begin
;
;  POWER spectra
;
           out_file = filenameroot+'.pds'
           andato=1
        end

   'f': begin
;
;  FOURIER spectra
;
           out_file = filenameroot+'.fft'
           andato=1
        end

   else: begin
            print,'Unrecognized output, try again'
         end

   endcase
endwhile
;-------------------------------------------------------------------
;  choose GRADE

;IF (parameters_known eq 0) THEN BEGIN
;   quali_grade=''
;   print,'Grades to accumulate    : 02'
;   print,format='($,"Enter units   -> ")'
;   read,quali_unita
;   IF(quali_grade eq '') THEN BEGIN
;      quali_grade = '02'
;   ENDIF
;ENDIF
;; Now parse input into array
;quanti = strlen(quali_grade)
;sel_grade=0
;FOR iu =0,quanti-1 DO BEGIN
;	IF(strmid(quali_grade,iu,1) eq '0') THEN sel_grade = sel_grade+1
;	IF(strmid(quali_grade,iu,1) eq '2') THEN sel_grade = sel_grade+3
;ENDFOR
sel_grade = 2
;-------------------------------------------------------------------
;
; choose output filename
;
if (parameters_known eq 0) then begin
   oufilename=''
   print,'Default output filename    : ',out_file
   print,format='($,"Enter output filename   -> ")'
   read,oufilename
   if(oufilename eq '') then begin
      oufilename=out_file
   endif
endif
;-------------------------------------------------------------------
; choose time resolution
;
chtres=''
if (parameters_known eq 0) then begin
;  interactive

   print,'Available time resolution: ',tress(0,0),' seconds'

   nbin    = 1L
   snbin   = ''
   print,format='($,"Enter rebinning factor in time (1) -> ")'
   read,snbin
   if(snbin ne '') then begin
      reads,snbin,nbin
   endif
  endif ELSE begin
   nbin=treb      ; from input line or parameter file
endelse

tres_fft = nbin * tress(0,0)     ; actual time resolution

;mu_2string,tres_fft,chtres
;
; choose number of points
;
if (parameters_known ne 0) then begin
;
; Check for period (real number input)
;
;   snp   = string(double(npds))  ; double for GDL compatibility ??? CRAZY!
   snp   = npds
   perio = strpos(string(npds),'.')
   if(perio ge 0) then begin    ; period found
        np = long(round(npds/tres_fft))
       endif else begin           ; period not found
        np = long(npds)
   endelse
endif else begin
   snp=''
   andato=0
   while (andato eq 0) do begin
      if(tress[0,0] le 1.0) then begin
	      ; Window timing mode
	      np = 8192L
	   endif else begin
	      ; Photon counting mode
	      np = 128L
	  endelse
	  lungh = np*tress(0,0)
      stringhina = ' ('+strtrim(lungh,2)+' s)'
      print,'Default # of points per FFT ',np,stringhina
      print,format='($,"Enter # of points (I) or time length (R) per FFT-> ")'
      read,snp
      if(snp eq '') then begin
          np=0L ; GHATS
    endif else begin
         ;
         ; Check for period (real number input)
         ;
         perio = strpos(snp,'.')
         if(perio ge 0) then begin    ; period found
            reads,snp,np2
            np = round(np2/tres_fft)
           endif else begin           ; period not found
            reads,snp,np
         endelse
      endelse

      if((np ge 1l) and (mu_power_of_two(np) eq 1)) then begin
         andato=1
        endif ELSE begin
            print,'Number of data points must be a power of two!'
      endelse
   endwhile
endelse

if((np lt 1l) or (mu_power_of_two(np) ne 1)) then begin
   massage,'Number of points must be a power of two!'
   retall
endif
;--------------------------------------------------------------------------
; Array allocation
;
rdata = fltarr(np)           ; light curve
pwr   = fltarr(np/2)         ; produced power spectrum
;--------------------------------------------------------------------------
; Additional variables
;
nffts       = 0l    ; must be long!
perc_old    = 0
ntotal_ffts = 0l
ntrafos_header_offset = 84L
T           = np*tres_fft
;--------------------------------------------------------------------------
; Compute final number of FFTs
;
for i=0,nfiles-1 do begin
   for j=0,nvaltimes(i,0)-1 do begin
      ntotal_ffts = ntotal_ffts + fix((valtimes2(j,i,0)-valtimes1(j,i,0))/T)
   endfor
endfor
ntotal_ffts = ntotal_ffts * proliferation    ; account for sliding window
;--------------------------------------------------------------------------
; Summary of input
;--------------------------------------------------------------------------
logfilename=filenameroot+'.log'
openw,1,logfilename
opsys = !version.os_family
if (opsys eq 'Windows') then begin
    dum = 'Win user'
   endif else begin
    spawn,'whoami',dum
endelse
prino,1,'----------------------------------------------------------------',''
prino,1,'Program gh_xte    ',systime()+'    User: '+dum
prino,1,'----------------------------------------------------------------',''
if (opsys eq 'unix') then begin
   spawn,'pwd',dum
   prino,1,'Working directory               : ',strtrim(dum,2)
endif
prino,1,'Input filename                  : ',strtrim(infilename,2)
prino,1,'Instrument                      : ',strtrim(telescope,2)+ $
                                         '/'+strtrim(instrument,2)
prino,1,'Observation ID                  : ',strtrim(obsid,2)
prino,1,'Source name                     : ',strtrim(sourcename,2)
prino,1,'Start date of observation       : ',strtrim(startdate,2)
prino,1,'',''
prino,1,'Channels for accumulation       : ',strtrim(string(canali(0)),2)$
                                            +'-'+strtrim(string(canali(1)),2)
prino,1,'Type of output                  : ',strtrim(chout,2)
prino,1,'Output filename                 : ',strtrim(oufilename,2)
prino,1,'Time resolution (s)             : ',strtrim(tress(0,0),2)
prino,1,'Nyquist frequency (Hz)          : ',strtrim(1.0/(2.0*tres_fft),2)
prino,1,'Number of points per FFT        : ',strtrim(string(np),2)
prino,1,'Corresponding to length (s)     : ',strtrim(string(T),2)
prino,1,'Total number of FFTs            : ',strtrim(string(ntotal_ffts),2)
prino,1,'Time step                       : ',strtrim(string(T/proliferation),2)
prino,1,'Window                          : ',strtrim(wind,2)
if(keyword_set(usergti)) then begin
   prino,1,'User GTI file                : ',strtrim(usergti,2)
endif
prino,1,'----------------------------------------------------------------',''
close,1
;--------------------------------------------------------------------------
; Main loop over input files
;=========================================================================

	if((wind eq 'Hamming') or (wind eq 'Triplet') or (wind eq 'Gauss') or (wind eq 'Kaiser')) then begin
		if(~keyword_set(wpar)) then begin
			print,'The selected window function requires an input parameter!'
		  endif else begin
			finestra = ghats_window(np,wind,par=wpar,winn=winn)
		endelse
	 endif else begin
		finestra = ghats_window(np,wind,winn=winn)
    endelse
;=======================end windowing=========================
;
; Definition and allocation of major arrays
;
time_offset     = double(0.0)
fields_selected = intarr(MAX_N_COL,nmetafiles)
mjreff          = double(0.0)

files_units     = intarr(nmetafiles)
files_tags      = files_units    ; separate bookkeeping for files
irowl           = lonarr(nmetafiles)
;  Here MU 3.0 -------------------
rate1 		= fltarr(MAX_VLE)
rate2 		= fltarr(MAX_VLE)
rate3 		= fltarr(MAX_VLE)
;
; MU6: tag assigned to each file to open
;
files_tags = indgen(nfiles,nmetafiles)


;-------------------------------------------------------------------
; If necessary, intersect user's input GTI files with the data GTIs
	if(keyword_set(usergti)) then begin
		mu_intersect_gtis,valtimes1,valtimes2,u_gti1,u_gti2,nvaltimes,  $
		                  MAX_N_VALTIMES,nfiles,nmetafiles
	endif
;---------------------------------
;
for i=0,nfiles-1 do begin    ;         ******** MAIN LOOP ********
;--------------------------------------------------------------------------
;    Obtain information for the header in output
;
   if (i eq 0) then begin
      for j=0,nmetafiles-1 do begin
         sw_first_file = 1
         mu_init,filenames(i,j),time_offset,sw_first_file, $
                    fields_selected,j,source,observatory, $
                    instrument,mjdrefi,mjdreff,timezero
      endfor
   endif ELSE begin
         sw_first_file = 0
         mu_init,filenames(i,0),time_offset,sw_first_file, $
                    fields_selected,0,source,observatory, $
                    instrument,mjdrefi,mjdreff,timezero
   endelse

;
; Obtain information about VLE counts
; (i.e. accumulate VLE rate light curve)
;
; New call for MU 3.0 with energy bands for colors
;
;   mu_get_housekeeping,filenames(i,0),tstarts(i,0),tends(i,0),housedata, $
;                 vle_times,vle_counts,npcus,MAX_VLE,nvle,         $
;		BANDE,rate1,rate2,rate3
;
; Open the relevant files in each metafile
;  and intialise its row number for use in the data_read routines
;
   for j=0,nmetafiles-1 do begin
      irowl(j)       = 1l
      mu_open_close,unita,filenames(i,j),'open'
; As a nice hidden feature, fxbopen gets its own logical units, so
; we need to take care of it. Our 12+j scheme turns out to be pointless
      files_units(j) = unita
   endfor
;
;  Set starting time for fft to zero
;
   tstart_fft       = 0.0d0
   tstart_fft_last  = -1.0d0
   tend_fft         = 0.0d0
   sw_first_segment = 1
;
;  Main loop for analysis within file(i)
;
   good             = 0
   quattrodieci:
;
;  Check the times (important routine!)
;
   mu_check_time,files_units,types,valtimes1,valtimes2, $
              nvaltimes,nfiles,nmetafiles,           $
              tstart_fft,tend_fft,T,i,good,instrument
;
   if(good eq -1) then goto,trecinquanta
;
;  Now we know the time interval to use, so get the corresponding data
;
   rdata = rdata*0.0

   for j=0,nmetafiles-1 do begin
	  unita = files_units(j)
      tag   = files_tags(i,j)
      cannali = channels(*,i,j)
      riga    = 0l
      mu_read_huiyan_le,unita,tag,rdata,np,cannali, $
                  fields_selected,j,tress(i,j),tres_fft,   $
                  tstart_fft,tend_fft,$
                  riga,channels_std21,channels_std22,channels_std23,nffts,energy,sel_grade
      irowl(j) = riga  ; subarrays are passed by value!!
   endfor     ; loop on metafiles (j)
;
;  Debug output
;
;   print,'T1, T2, Cnts: ',tstart_fft,tend_fft,   $
;          total(rdata)/n_elements(rdata)/tres_fft
;
;  time offset is obtained as follows from the FITS file header:
;  time_offset = dble(mjdrefi) + mjdreff + dble(timezero)/86400.0d0
;  tstart_fft is Mission Elapsed Time (MET) in raw spacecraft clock seconds,
;  to an accuracy of 4 s. By adding timezero to it, this becomes accurate to
;  100 usec. To make it even more accurate, to 5 usec, you have to add the
;  clock correction, which is NOT in the FITS file. start_time, obtained by
;  adding MJDREF to MET, is therefore MJD (JS-240000.5) on the TT time scale,
;  accurate to 100 usec.
;
   start_time = tstart_fft/86400.0d0 + time_offset
   if(tstart_fft le tstart_fft_last) then begin
      print,'WARNING: time reversal! Last: ',tstart_fft_last,' now: ',  $
            tstart_fft,' seconds'
      trtc
   endif
   tstart_fft_last = tstart_fft
;==========================================================================
; FFT section
;==========================================================================
;
;  ****  HERE THE GOING GETS TOUGH!
;  The FFT is actually done in here!
;
;pippo ;****************************************************************************************************
   crate = total(rdata)/n_elements(rdata)/tres_fft
   mu_compute_fft,rdata,pwr,np,T,start_time,source,observatory,           $
          instrument,oufilename,chout,1,0,                   $
          0.0,tstart_fft,new_output_flag,                    $
      output_unit,ntotal_ffts,canali,                        $
	  std21,std22,std23, $
	  BANDE,turbo,gti_flag,proliferation     ; MU6
   std21 = 0.0
   std22 = 0.0
   std23 = 0.0
;==========================================================================
   nffts = nffts + 1
   perc = fix(100.0*nffts/ntotal_ffts)
   print,'FFT# ',nffts,'/',ntotal_ffts,' Rate: ',crate
;  Set up new time interval
;
   if(keyword_set(sliding)) then begin ; here sliding is implemented
	  tstart_fft = tstart_fft + T /sliding
   endif else begin
	  tstart_fft = tend_fft   ; in case no sliding is requested
   endelse
   sw_first_segment = 0
   goto,quattrodieci
;
;  Second loop return
;
;  Close all relevant files
;
   trecinquanta:
;   close,/all      ; Ultimate choice: why bother?
   for j=0,nmetafiles-1 do begin
      unita=files_units(j)
      mu_open_close,unita,filenames(i,j),'close'
   endfor
endfor    ; loop on nfiles
;
; Here sanity check on nfft: did it reach the end? If not, something must
; be done to recovery header
;
if(nffts ne ntotal_ffts) then begin
   print,'Readjusting number of ffts to ',nffts
   point_lun,output_unit,ntrafos_header_offset
   writeu,output_unit,nffts
endif
;--------------------------------------------------------------------------
; End of program
;
close,/all
fine_dps:
print,'gh_huiyan_le: normal termination'
;beep
;
; Standard plot on a PS file
;
if(keyword_Set(gh)) then begin
   ghats_all,oufilename,/PS,/poisson
endif

end
