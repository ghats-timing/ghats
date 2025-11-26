pro gh_xmm,infilename,outtype,canali,treb,npds,oufilename,ghx=ghx,turbo=turbo,bands=bands,gti=usergti,sliding=sliding,bary=bary, $
	       wind=wind,wpar=wpar,mingap=mingap
;+
; NAME:
;      GH_XMM
; PURPOSE:
;      Production of XMM PDS and FFT files for the GHATS project
; EXPLANATION:
;      This procedure, coming from MU, produces
;      PDS and FFT files from Swift data
;
; CALLING SEQUENCE:
;       GH_XMM
;       GH_XMM,parfilename
;       GH_XMM,infilename,outtype,channels,treb,npds,outfilename[,/muxana][,/turbo][,bands=B][,gti=GTI][,sliding=SLIDING][,/BARY][,mingap=MINGAP]
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
;		TURBO       = if set, the PDS is produced through Craig Markwardt's RADPS routine. It is 
;					  faster, but it only works for a few missions
;		BANDS       = energy channels for the production of three standard light curves. It is an array of six
;			  		  elements ([0,10,11,20,21,40]) means three curves in the channel ranges 0-10, 11-20 and
;			 		  21-40 respectively. They can then be used to produce X-ray hardnesses. Default values
;					  are [70,500,501,4000,4001,10000].
;		GTI         = optional user-supplied GTIs, to be intersected with the ones in the data files. 
;			          If set, it must be a valid filename. Two possible input formats are possible. Either
;					  an ASCII file with lines containing start end time of each GTI (in Swift time units),
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
;       MINGAP      = maximum gap duration to fill in with local average 
;
;
; OUTPUTS:
;       NONE
;
; EXAMPLE:
;       GH_XMM
;         Prompts user interactively for parameters
;       GH_XMM,'#mysource.par'
;         Takes parameters from the specified parameter file
;       GH_XMM,'@binned.lis','POWER',[0,35],1,4096,'mysource.pds'
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
;       T. Belloni/S. Motta  9 May 2012  from GH_XTE
;		T. Belloni			 7 Sep 2012  GTI warning added
;       T. Belloni/S. Motta 11 Sep 2012  small gap filling implemented
;		T. Belloni  14 Dec 2013  fixed time rebinning from I to L
;		T. Belloni  26 Nov 2015  fixed npds in command line mode
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
;
; Check whether turbo flag is set.
;
if(keyword_Set(turbo)) then begin
;	if(sistema eq 'IDL') then begin
;   		turbo = 1
;    endif else begin
	    turbo = 0
	    ;massage,'Turbo mode is not enabled under GDL. Sorry.'
        massage,'Turbo mode is not enabled for XMM. Sorry.'
;	endelse
endif else begin
   turbo = 0
endelse
; MU6 turbo
common dati,tempi,eventi,canale_corrente,std21,std22,std23 ; common block for keeping the data
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
MAX_N_VALTIMES  =  10000
MAX_N_CHANNELS  =  32768   ; additional channel to collect all channels > 15k
MAX_N_COL       =  256
MAX_VLE         = 6000    ; number for rates

;================================================================
if(keyword_set(bands)) then begin
	BANDE           = bands
endif else begin
    BANDE           = [70,500,501,4000,4001,10000]
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
	endif else begin
	print,'WARNING! GH_XMM does not read the GTI information in the event files(s)'
	print,'         You should supply your own extrnal GTIs'
	;exit
endelse

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
        if((turbo eq 1) and (outtype eq 'FFT')) then begin
	        massage,'Turbo mode not enabled for FFT'
	        stop
	    endif
      end

   6: begin
     	parameters_known = 1
        if((turbo eq 1) and (outtype eq 'FFT')) then begin
	        massage,'Turbo mode not enabled for FFT'
	        stop
	    endif
      end

   else: begin
     	print,'Usage: gh_xmm,infile,outtype,channels,treb,npds,oufile'
     	stop
      end
endcase
;-------------------------------------------------------------------
;
;   get input files through a GUI window for IDL and a text list for GDL
;
if (parameters_known eq 0) then begin
   mu_ask,infilename,'*'
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
mu_get_time_resolution_xmm,filenames,nfiles,nmetafiles,tress
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
mu_get_header_information_xmm,filenames,types,tstarts,tends,tress,valtimes1,valtimes2, $
     nfiles,nmetafiles,nvaltimes, $
         sourcename,obsid,startdate,telescope,instrument
;-------------------------------------------------------------------
;
; obtain channel information
;
channels1 = intarr(MAX_N_CHANNELS,nfiles,nmetafiles)
channels2 = channels1
nchannels = intarr(nfiles,nmetafiles)
mu_get_channel_information_xmm,channels1,channels2,nchannels,nfiles,nmetafiles
;-------------------------------------------------------------------
channels = intarr(MAX_N_CHANNELS,nfiles,nmetafiles)*0
if (parameters_known eq 0) then begin
;
; select channels interactively
;
   mu_channel_selection_xmm,metafiles,filenames,channels,channels1,channels2, $
                   nchannels,nfiles,nmetafiles,istart,iend
   canali(0) = istart
   canali(1) = iend
endif ELSE begin
;
; simply use the available range
;
   canali(0) = 0
   canali(1) = 14999
   mu_set_all_channels_xmm,metafiles,filenames,channels,channels1,channels2, $
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
;      if(tress[0,0] le 1.0) then begin
;	      ; Window timing mode
;	      np = 8192L
;	   endif else begin
;	      ; Photon counting mode
;	      np = 128L
;	  endelse
      np = 65536L
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
;
; Here section where Craig's RADPS is called. It works, but
; only under very precise circumstances, which means I cannot use it. Moreover,
; it is not going to work for ASTROSAT data.
;
if (turbo eq 1 and chou2 eq 'p') then begin
	; TURBO NOT CHANGED FROM RXTE VERSION!  **** IT DOES NOT WORK
	print,'----------------------------------------------------------------'
	print,'                **** Turbo speed activated! ****'
	print,'----------------------------------------------------------------'
; Check whether a window file is required
    if((keyword_set(wind) and (wind ne 'Boxcar'))) then begin
	     print,'Windowing is not available in turbo mode!'
	     return
	end
; Call to Craig's routine'
   if(keyword_set(usergti)) then begin
	; assemble single array for user supplied GTIs
	turbo_gti = transpose([[u_gti1],[u_gti2]])
	radps,filenames,np*tres_fft,dps,steptime=np*tres_fft/proliferation, $
	      avgtime=np*tres_fft/proliferation,tbinsize=tres_fft,climits=canali, $
	      quiet=1,p0=licu,time=times,freqavg=nu,nspecsum=numerofft,status=stato, $
	      gti=turbo_gti,/nocrossgti,quality=q,exposure=e,minfracexp=1.0d0
   endif else begin
	radps,filenames,np*tres_fft,dps,steptime=np*tres_fft/proliferation, $
	      avgtime=np*tres_fft/proliferation,tbinsize=tres_fft,climits=canali, $
	      quiet=1,p0=licu,time=times,freqavg=nu,nspecsum=numerofft,status=stato, $
	      /nocrossgti,quality=q,exposure=e,minfracexp=1.0d0
   endelse
;
;  Cleaning of bad (not completely exposed intervals)

buoni     = where(q eq 0)   ; identify good spectra
dps       = dps(*,buoni)
times     = times(buoni)
licu      = licu(buoni)
numerofft = n_elements(times)
;
; Now needs to accumulate the housekeeping information, not included in
; radps

; WORK IN PROGRESS
vle_times       = dblarr(MAX_VLE)
vle_counts      = fltarr(MAX_VLE)
npcus           = intarr(MAX_VLE)
rate1 	    	= fltarr(MAX_VLE)
rate2 	     	= fltarr(MAX_VLE)
rate3 		    = fltarr(MAX_VLE)
nvle_total      = 0
;
; Major constants
;
MAX_TYPES    =   5
C_TYPES      =   6
MAX_COL      = 256
MAX_FILES    =  40
;
; Prepare keyword names for VLE rates (per PCU)
;
datatypes    = strarr(MAX_TYPES)
datatypes(0) = 'VLECntPcu0'
datatypes(1) = 'VLECntPcu1'
datatypes(2) = 'VLECntPcu2'
datatypes(3) = 'VLECntPcu3'
datatypes(4) = 'VLECntPcu4'
;
; Prepare keyword names for PCU2  columns
;
datatypesc   = strarr(C_TYPES)
datatypesc(0) = 'X1LSpecPcu2'
datatypesc(1) = 'X1RSpecPcu2'
datatypesc(2) = 'X2LSpecPcu2'
datatypesc(3) = 'X2RSpecPcu2'
datatypesc(4) = 'X3LSpecPcu2'
datatypesc(5) = 'X3RSpecPcu2'
;
; Array to contain column numbers corresponding to VLE rates
;
vlecol       = intarr(MAX_TYPES)
;
; Array to contain column numbers corresponding to PCU2 rates
;
ccol         = intarr(C_TYPES)

; Added 7-FEB-2012 
dirname=file_dirname(filenames(0,0))
std2_filelist   = file_search(dirname+'/FS4a*')

; check whether any STD2 files have been found. If not, skip the section and have
; all data seto to 0 to avoid a crash
if std2_filelist[0] ne '' then begin
   po              = strpos(std2_filelist,'_bkg')
   std2_filelist    = std2_filelist(where(po eq -1)) ; remove bkg files

   for istd2=0,n_elements(std2_filelist)-1 do begin ; loop over the std2 files
   		;
		; Open the file
		;
		;print,std2_filelist(istd2)  ; test
		fxbopen,unit,std2_filelist(istd2),1,header,errmsg=errmsg
		;
		; Read in basic header information
		;
		t1        = fxpar(header,'TSTART')
		t2        = fxpar(header,'TSTOP')
		timedel   = fxpar(header,'TIMEDEL')
		datamode  = strtrim(fxpar(header,'DATAMODE'))
		tdim      = fxpar(header,'TDIM*')
		nfields   = strtrim(fxpar(header,'TFIELDS'))
		ttype     = fxpar(header,'TTYPE*')
		nrows     = fxpar(header,'NAXIS*')
		idltype   = intarr(nfields)
		;
		; Identifies which column numbers correspond to the VLE and Time columns
		;
		timecol = -100  ; marker for column number corresponding to Time
		fxbtform,header,tbcol,idltype,formato,numval,maxval
		ttype       = fxpar(header,'TTYPE*')
		for i=0,nfields-1 do begin
	   		for j=0,4 do begin
	      		if(ttype(i) eq datatypes(j)) then vlecol(j) = i+1
	   		endfor
	   		if(ttype(i) eq 'Time    ') then timecol = i+1
		endfor
		;
		; Identifies which column numbers correspond to the PCU2 columns
		;
		for i=0,nfields-1 do begin
	   		for j=0,5 do begin
	      		if(ttype(i) eq datatypesc(j)) then ccol(j) = i+1
	   		endfor
		endfor

    	nvle       = nrows(1)
		if(nvle_total+nvle gt MAX_VLE) then begin
	   		massage,'Insufficient array size for vle counts'
	   		retall
		endif
		;
		; Now loop over the STD2 rows
		;
		for irow=0,nvle-1 do begin
	   		npcus(nvle_total+irow)      = 0
	   		vle_counts(nvle_total+irow) = 0.0
			;
			; For each row, read in time and vle rates
			; The assumption is that if the VLE counts for a PCU are > 1, the PCU is on
	   		fxbread,unit,timel,timecol,irow+1
	   		vle_times(nvle_total+irow) = timel
	   		for j=0,4 do begin
	      		fxbread,unit,icounts,vlecol(j),irow+1
	      	if (icounts gt 1) then begin
	         	npcus(nvle_total+irow) = npcus(nvle_total+irow)+1
	      	endif
	      	vle_counts(nvle_total+irow) = vle_counts(nvle_total+irow) + icounts
	   	endfor
		;
		;  Here accumulate the three PCU2 light curves
		;
	   	spettro = intarr(129)
	   	for j=0,5 do begin
			;       Add up all layers, left and right
	      	fxbread,unit,spe,ccol(j),irow+1
	      	spettro =  spettro + spe
	   	endfor
	   	rate1(nvle_total+irow) = float(total(spettro(BANDE(0):BANDE(1))))
	   	rate2(nvle_total+irow) = float(total(spettro(BANDE(2):BANDE(3))))
	   	rate3(nvle_total+irow) = float(total(spettro(BANDE(4):BANDE(5))))
		;
		;  Convert VLE counts into counts per second (VLE rate)
		;
	   	vle_counts(nvle_total+irow) = vle_counts(nvle_total+irow)/timedel
		endfor
	
		nvle_total = nvle_total + nvle
		;
		; Close standard2 file
		;
		fxbclose,unit
	endfor
endif else begin
	       print,'No STD2 data found!'
endelse

vle_times  = vle_times(0:nvle_total-1)
vle_counts = vle_counts(0:nvle_total-1)
npcus      = npcus(0:nvle_total-1)
rate1      = rate1(0:nvle_total-1)
rate2      = rate2(0:nvle_total-1)
rate3      = rate3(0:nvle_total-1)

;====================================================================
	
; Assemble mjd
	sanno   = strmid(startdate,0,4)
	smese   = strmid(startdate,5,2)
	sgiorno = strmid(startdate,8,2)
	sora    = strmid(startdate,11,2)
	sminuto = strmid(startdate,14,2)
    ssecondo = strmid(startdate,17,2)
    reads,sanno,anno
    reads,smese,mese
    reads,sgiorno,giorno
    reads,sora,ora
    reads,sminuto,minuto
    reads,ssecondo,secondo
    quandoera = [anno,mese,giorno,ora,minuto,secondo]
    juldate,quandoera,rmjd
    rmjd=double((rmjd)-0.5)
; Array handling
    n_frequenze = 2*long(n_elements(nu))
    n_tempi     = long(n_elements(licu))
    telescopio  = string(telescope,format='(a16)')
    sorgente    = string(sourcename,format='(a16)')
    strumento   = string(instrument,format='(a16)')
	tempo0  = times(0)
; Main loop to write th PDS file
	for ncraig=0,n_elements(licu)-1 do begin
		;
		; Here accumulate HK information and rates for each line
		;
		tempo_centrale = times(ncraig)+T/2.0
		mu_get_vle_rate,tempo_centrale,current_vle_rate, $
		                ndet,vle_times,vle_counts,npcus,nvle_total
		mu_read_std2_rates,tempo_centrale,std21,std22,std23, $
		                   vle_times,rate1,rate2,rate3,nvle_total
		;
		cct     = float(licu(ncraig)*np*tres_fft)
		potenza = float(dps(*,ncraig))
		ttt     = rmjd+(times(ncraig)-tempo0)/86400.0d
		
		mu_zhang,cct,T,ndet,n_frequenze,current_vle_rate,i_vle,poi
		
	mu_write_single_power,2.0*cct,potenza,n_frequenze,poi,cct, $
	       T,ttt,oufilename,ndet,current_vle_rate,new_output_flag, $ 
	       output_unit,n_tempi,canali, $
	       telescopio, $
	       sorgente,$
	       strumento,std21,std22,std23,  $
		   BANDE,turbo,gti_flag,proliferation,i_vle     ; MU6
	endfor
	close,/all
    if(stato eq 1) then begin
	    goto,fine_dps
	endif else begin 
	   print,'RADPS error!'
	   return
	endelse
endif
;========================NON TURBO===============================
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
;	if(keyword_set(usergti)) then begin
;		mu_intersect_gtis,valtimes1,valtimes2,u_gti1,u_gti2,nvaltimes,  $
;		                  MAX_N_VALTIMES,nfiles,nmetafiles
;	endif

; Identify GTIs shorter than the user-specified limit
   if(keyword_set(mingap) eq 0) then begin
      mingap=0
   endif
   gti_duration = u_gti2 - u_gti1
   to_cut  = where(gti_duration lt mingap)
   to_cut2 = where(gti_duration lt mingap)+1
; Reconstruct GTI array from user GTIs, throwing away original GTIs
   valtimes1 = dblarr(MAX_N_VALTIMES,nfiles,nmetafiles)
   valtimes2 = valtimes1

   IF(to_cut(0) ge 0) THEN BEGIN
   	toomuch = where(to_cut2 gt (n_elements(u_gti1)-1))
   	to_cut2(toomuch) = where(toomuch-1)

   	gap1 =  u_gti1(to_cut)
   	gap2 =  u_gti2(to_cut)

   	u_gti2(to_cut)   = -u_gti2(to_cut)
   	u_gti1(to_cut+1) = -u_gti1(to_cut+1)
   
; Cut short GTIs
   	u_gti1 = u_gti1(where(u_gti1 ge 0))
   	u_gti2 = u_gti2(where(u_gti2 ge 0))

   ENDIF ELSE BEGIN

    for jjj = 0,nmetafiles-1 do begin
       for kkk = 0,nfiles-1 do begin
           valtimes1(0:n_elements(u_gti1)-1,kkk,jjj) = u_gti1
           valtimes2(0:n_elements(u_gti1)-1,kkk,jjj) = u_gti2
           nvaltimes(kkk,jjj)   = n_elements(u_gti1)
       endfor
    endfor
; Trim
    mass = max(nvaltimes)
    valtimes1 = valtimes1(0:mass-1,*,*)
    valtimes2 = valtimes2(0:mass-1,*,*)

    gap1 = -1
    gap2 = -1

   ENDELSE

;---------------------------------
;
for i=0,nfiles-1 do begin    ;         ******** MAIN LOOP ********
;--------------------------------------------------------------------------
;    Obtain information for the header in output
;
   if (i eq 0) then begin
      for j=0,nmetafiles-1 do begin
         sw_first_file = 1
         mu_init_xmm,filenames(i,j),time_offset,sw_first_file, $
                    fields_selected,j,source,observatory, $
                    instrument,mjdrefi,mjdreff,timezero
      endfor
   endif ELSE begin
         sw_first_file = 0
         mu_init_xmm,filenames(i,0),time_offset,sw_first_file, $
                    fields_selected,0,source,observatory, $
                    instrument,mjdrefi,mjdreff,timezero
   endelse
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
      mu_read_swift,unita,tag,rdata,np,cannali, $
                  fields_selected,j,tress(i,j),tres_fft,   $
                  tstart_fft,tend_fft,$
                  riga,channels_std21,channels_std22,channels_std23,nffts,energy
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
   crate = total(rdata)/n_elements(rdata)/tres_fft
   mu_compute_fft_xmm,rdata,pwr,np,T,start_time,source,observatory,           $
          instrument,oufilename,chout,1,0,                   $
          0.0,tstart_fft,new_output_flag,                    $
      output_unit,ntotal_ffts,canali,                        $
	  std21,std22,std23,                                     $
	  BANDE,turbo,gti_flag,proliferation,                    $
	  gap1,gap2,mingap
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
   openu,uu,oufilename,/get_lun
   point_lun,uu,84
   writeu,uu,nffts
   free_lun,uu
endif
;--------------------------------------------------------------------------
; End of program
;
close,/all
fine_dps:
print,'gh_xmm: normal termination'
;beep
;
; Standard plot on a PS file
;
if(keyword_Set(gh)) then begin
   ghats_all,oufilename,/PS,/poisson
endif

end
