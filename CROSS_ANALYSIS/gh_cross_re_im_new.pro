;+
; NAME: 
;      GH_LOWCOH_SSQ
; PURPOSE: 
;      Calculate the necesarry input to obtain Coherence error estimates
;      in the High Powers, Low Measured Coherence case of VN97.
; NOTES:
;      This functions are used by GH_CROSS_RE_IM_NEW.
; MODIFICATION HISTORY: 
;      F. Garcia  22 Mar 2024  from scratch


FUNCTION PROBDIST_CROSS, ssq, asq=asq, nsq=nsq
    ; Probability distribution of SSQ (VN97 eq 6), typo corrected (s*n)^-1.
    a = sqrt(asq)
    n = sqrt(nsq)
    s = sqrt(ssq)
    axn = 2.0*a*s/nsq < 700
    eexp = EXP(-0.5*(asq+2*ssq)/nsq)
    b1 = BESELI(axn, 0)
    maann = 0.5*asq/nsq < 700
    b2 = BESELI(maann, 0)
    ret = eexp*b1/b2/sqrt(!pi)/(n*s)
    return, ret
END

FUNCTION INVERSE_CDF, ssq, asq=asq, nsq=nsq, LVL=LVL
    ; Calculate CDF confidence interval by integrating a PDS.    
    S = 0
    FOR i = 0, 33 DO BEGIN
        TRAPZD,'PROBDIST_CROSS', 1e-6, ssq, S, 1000, asq=asq, nsq=nsq
    ENDFOR
    return, S - LVL
END

FUNCTION gh_lowcoh_ssq_CI,CL,asq,nsq
    ; Calculate the Confidence Limits of SSQ (VN97 eqn 9)
    if (asq/nsq GT 300) then begin
    ; if S/N is large, use the Gaussian approximation (VN97 eqn 7)
       lower_bound = asq-sqrt(2*nsq*asq)
       upper_bound = asq+sqrt(2*nsq*asq)
    endif else begin
    ; Confidence levels as roots of the CDF (VN97 eqn 6)
        lower_bound = ZBRENT(1e-6,asq+100*nsq,FUNC='INVERSE_CDF', asq=asq, nsq=nsq, LVL=0.5*(1-CL))
        upper_bound = ZBRENT(1e-6,asq+100*nsq,FUNC='INVERSE_CDF', asq=asq, nsq=nsq, LVL=0.5*(1+CL))
    endelse
    return, [lower_bound, upper_bound]
END

FUNCTION gh_lowcoh_ssq,asq,nsq
     ; Analytic formula for average SSQ (eqn 3.20 Chakrabarty 1995)
     xsq = asq/nsq    
     mxsq = 0.5*xsq < 700
     return,0.5*nsq*(1+xsq+xsq*BESELI(mxsq,1)/BESELI(mxsq,0))  
END

pro gh_cross_re_im_new,filename1,filename2,             $
             frequency_reb,lag,lag_err,coherence,coherence_err, $
             irf,tlag=tl, $
             repart,repart_err,impart,impart_err,modcs,modcs_err,mw,qf, $
             index=index,time=time,sel=sel,rate=rate,poisson=poisson,poivalue=poivalue, $
             manual_poi=manual_poi, $
             rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh,rawcoh=rawcoh, $
             status_quiet=status_quiet,help=help
;+
; NAME: 
;      GH_CROSS_RE_IM_NEW
; PURPOSE: 
;      Reads is two FFTs from FFT files, selecting on a range of FFTs
;			in time, range, or selection of indices
;      This procedure extracts two FFTs from two compatible files. The
;      selection is possible only on a continuous range of FFTs.
;      Produces phase lags (in radians) or time lags, and coherence.
;      It returns also the Re/Im parts of the CS
;
; CALLING SEQUENCE: 
; (This is obsolete; see examples below)
;       MUCROSS,FILENAME1,FILENAME2,FREQUENCY,LAG,LAG_ERR,
;                 COHERENCE,COHERENCE_ERR,IRF,
;			      [,INDEX=INDEX][,TIME=TIME][,SEL=SEL]
; INPUTS:
;       FILENAME1 = name of the first  input FFT file
;       FILENAME2 = name of the second input FFT file
;       STARTING = index of first FFT to be selected
;       ENDING   = index of last  FFT to be selected
;       IRF      = rebin factor (negative for log rebinning)
;
; OUTPUTS:
;       FREQUENCY    = Frequency array
;       LAG          = Phase/timelag array
;       LAG_ERR      = Array of errors on phase/timelags
;       COHERENCE    = Coherence function
;       COHERENCE_ERR= Error on coherence function
;       REPART       = Real part of the cross spectrum
;       REPART_ERR   = Error of the eeal part of the cross spectrum
;       IMPART       = Imaginary part of the cross spectrum
;       IMPART_ERR   = Error of the imaginaryi part of the cross spectrum
;
; KEYWORDS:
;       TIME         = Switch for time lags (default phase lags)
;       POIVALUE     = Real part of noise average cross-spectrum (for subtraction)
;       POISSON      = MM: I think this is a 2D array that has start/end freq. to 
;                      compute the Poisson level; for instance [800,1200]
;       MANUAL_POI   = Fixed PDS1/PDS2 noise levels in input-file units.
;       RAWCOH       = Return measured/raw coherence, using the observed PDS
;                      in the denominator and no Poisson/noise correction.
;       RMS1         = background value for rms2 conversion for band1. If specified, 
;                      rms2 conversion will be performed
;       RMS2         = background value for rms2 conversion for band2. If specified, 
;                      rms2 conversion will be performed
;       ROTATE       = Switch to rotate the cross vect5or by X (currently pi/4) rad
;                      The user should remember to subtract X rad from the phase lags.
;                      Default is no rotate
;
; EXAMPLE:
;       Read in two FFTs and compute phaselag spectrum
;
;       MU> MUCROSS_N,'data/4u1630a.fft','data/4u1630b.fft',nu,pha,pha_e, $
;                      coh,coh_e,1,10000,-100
;
; EXAMPLE 2:
;       First read in two FFTs (R=Reference,band  S=Subject band)
;
;       Ghats> filer='fileR.fft' & files='fileS.fft''
;       Ghats> gh_cross_re_im_new,filer,files,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,rms1=0.001,rms2=0.001,/rotate
;
;       Then find the average of the Real part of the cross spectrum over a frequency range dominated by Poisson noise (in this case 400-800 Hz)
;
;       Ghats> poivalue=mean(real(where(freq ge 400 and freq le 800)))
;
;       And finally subtract that average real value from the Real part of the cross spectrum
;       (This corrects both for cross talk and for the partial correlation if the subject band is part of te reference band)
;
;       Ghats> gh_cross_re_im_new,filer,files,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,rms1=0.001,rms2=0.001,poivalue=poivalue,/rotate
;
;       After that, for instance, write the Real and Imaginary parts of the cross spectrum, the lags and the coherence function in Xspec
;
;       Ghats> gh_xspec,freq,real,real_err,'real'
;       Ghats> gh_xspec,freq,imag,imag_err,'imag'
;       Ghats> gh_xspec,freq,lag,lag_err,'lag'
;       Ghats> gh_xspec,freq,coh,coh_err,'coh'
;
;       Notice that For non-PCA data, no noise level is subtracted unless POISSON=[f1,f2] or MANUAL_POI=[p1,p2] is given.
;       To use the Poisson level correctly, add the frequency range to average the PDS, e.g.:
;       
;       Ghats> gh_cross_re_im_new,filer,files,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,rms1=0.001,rms2=0.001,poivalue=poivalue,poisson=[400,800],/rotate
;       
;       (MM: I am not 100% sure this is how it works; I deduced this from the code, but I would have to test it)

;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       MUXANA_OPENFFT:    Opens a FFT file
;       MUXANA_GETHEADER:  Reads in the header
;       READ_FFT_LINE:     Reads in next line from FFT file
;       REBINCROSS:        Array rebinning
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  30 Sep 2002  from MUXANA_N
;       T. Belloni  01 Oct 2002  coherence and time lag switch added
;       T. Belloni  20 Jan 2003  changed ATAN for lag determination in IDL 5.5
;       T. Belloni  11 Dec 2003  mu 3.0 with colors
;       T. Belloni  29 Mar 2009  I_VLE now read in
;		T. Belloni  01 Apr 2010  new version from Uttley's code
;	    T. Belloni  07 Apr 2010  mucross from mucross_n and mux
;		T. Belloni  03 Dec 2010  from mu6
;		T. Belloni   9 Feb 2012  time range bug fixed
;		T. Belloni   9 Feb 2012  /TIMELAG keyword changed to /TLAG
;		T. Belloni  20 Feb 2012  no floor rounding for times array
;		T. Belloni   7 Jun 2012  coherence calculation had Poisson estimate only valid for RXTE/PCA
;       T. Belloni  12 Jun 2012  increased tolrance for time check 
;		T. Belloni   5 Dec 2013  free_lun replaces close
;		T. Belloni   7 Feb 2017  poisson keyword added (and poissonian bug fixed)
;		T. Belloni  24 May 2017  corrected poilevel=2
;		T. Belloni  30 Jun 2020  added rate selection
;		T. Belloni  16 Jul 2020  absolute value of difference for time compatibility check
;  		T. Belloni  26 Feb 2013  added keyword for input real part to subtract (a la gh_cross_range)
;       M. Mendez   21 Apr 2021 Output Re/Im parts and coresponding errors (in Leahy units)
;       M. Mendez    6 Jun 2023 Option to output Re/Im parts and errors in rms units
;       M. Mendez    9 Jun 2023 Option to rotate the cross vector to make fit of CDS more stable
;       F. Garcia   20 Mar 2024 Slight corrections to coherence error calculation from VN97 (eq. 8)
;       F. Garcia   20 Mar 2024 Added option to calculate low coherence case from VN97 (eq. 9)
;       M. Mendez   15 May 2026 Added /HELP and returned cross-spectrum modulus
;                               with covariance-propagated error; developed with
;                               assistance from ChatGPT.
;       M. Mendez 16 May 2026 renamed manual noise keyword to MANUAL_POI
;                              and clarified that values are in input-file
;                              units before RMS1/RMS2 conversion;
;                              developed with assistance from ChatGPT;
;       M. Mendez 15 Jul 2026 added /RAWCOH to force measured coherence
;                              with observed PDS denominator and no noise
;                              correction; developed with assistance
;                              from ChatGPT
;       M. Mendez 16 Jul 2026 added private STATUS_QUIET keyword for
;                              suppressing noise-prescription status during
;                              GH_CS preliminary calls; developed with
;                              assistance from ChatGPT
;--------------------------------------------------------------------------
; MM
if(keyword_set(help)) then begin
   print,''
   print,'GH_CROSS_RE_IM_NEW'
   print,''
   print,'Compute lags, coherence, Re/Im parts, and modulus'
   print,'of the cross spectrum from two FFT files.'
   print,''
   print,'Usage:'
   print,"  gh_cross_re_im_new,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf"
   print,"  gh_cross_re_im_new,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf,poisson=[400,800]"
   print,"  gh_cross_re_im_new,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf,manual_poi=[p1,p2]"
   print,''
   print,'Keywords:'
   print,'  /TLAG                 Return time lags'
   print,'  POISSON=[f1,f2]       Frequency range for PDS noise'
   print,'                         subtraction and Real correction'
   print,'  MANUAL_POI=[p1,p2]    Manual PDS1/PDS2 levels in input-file units'
   print,'                         before RMS1/RMS2 conversion'
   print,'                         (coherence only; no Real correction)'
   print,'  POIVALUE=             Subtract constant Real component'
   print,'  RMS1=, RMS2=          Convert Re/Im/Mod to rms^2'
   print,'  /ROTATE               Rotate cross vector by pi/4'
   print,'  /LOWCOH               Use VN97 low-coherence formalism'
   print,'  /RAWCOH               Return measured/raw coherence: observed PDS'
   print,'                         denominator, no Poisson/noise correction'
   print,'  INDEX=, TIME=, RATE=, SEL='
   print,'                         Data selection'
   print,'  /HELP                 Print this message'
   print,''
   return
endif
;
; MM 16 May 2026
; Do not allow two user-requested noise prescriptions.
;
if(keyword_set(poisson) and keyword_set(manual_poi)) then begin
   print,'ERROR: Use only one of POISSON or MANUAL_POI.'
   print,'POISSON=[f1,f2]: frequency range for PDS + Real correction.'
   print,'MANUAL_POI=[p1,p2]: fixed PDS1/PDS2 levels only.'
   retall
endif
; MM
;
; RAWCOH is a separate scientific definition: measured coherence with
; observed PDS in the denominator. Keep it mutually exclusive with
; noise-corrected coherence modes and the low-coherence estimator.
;
if(keyword_set(rawcoh) and keyword_set(poisson)) then begin
   print,'ERROR: /RAWCOH cannot be combined with POISSON.'
   print,'RAWCOH uses the observed PDS denominator with no noise correction.'
   retall
endif

if(keyword_set(rawcoh) and keyword_set(manual_poi)) then begin
   print,'ERROR: /RAWCOH cannot be combined with MANUAL_POI.'
   print,'RAWCOH uses the observed PDS denominator with no noise correction.'
   retall
endif

if(keyword_set(rawcoh) and keyword_set(lowcoh)) then begin
   print,'ERROR: /RAWCOH cannot be combined with /LOWCOH.'
   print,'LOWCOH requires non-zero noise terms for the low-coherence estimator.'
   retall
endif
; MM
;
; Open first FFT file
;
ghats_openfft,filename1,unit1,/dialog
ghats_openfft,filename2,unit2,/dialog
;
; Read in fft file header
;
ntrafos1   = 0l
ntrafos2   = 0l
gh_version_string    = '                '
observatory1         = '                '
instrument1          = '                '
target1              = '                '

observatory2         = '                '
instrument2          = '                '
target2              = '                '
dummy     = bytarr(100)
rmjd01              = 0.0D0
rmjd02              = 0.0D0
;ghats_getheader,unit1,nft1,T1,observatory1,target1,rmjd01,instrument1, $
;                 ntrafos1,e1,dummy
;mufla1 = dummy(0)
;ghats_getheader,unit2,nft2,T2,observatory2,target2,rmjd02,instrument2, $
;                 ntrafos2,e2,dummy
;mufla2 = dummy(0)

ghats_getheader,unit1,gh_version_string,observatory1,instrument1,target1,rmjd01, $
                     nft1,T1,ntrafos1,e1,proliferation1,baryflag1,n_spectral_bins1, $
                     background_flag1,dummy
mufla1 = dummy(0)
ghats_getheader,unit2,gh_version_string,observatory2,instrument2,target2,rmjd02, $
                     nft2,T2,ntrafos2,e2,proliferation2,baryflag2,n_spectral_bins2, $
                     background_flag2,dummy
mufla2 = dummy(0)

i_vle = fix([dummy(17:18)],0)+1
;
; Checks, checks, checks
;
if(strcmp(observatory1,observatory2) eq 0) then begin
   massage,'Observatory keyword not compatible'
   retall
endif
if(strcmp(target1,target2) eq 0) then begin
   massage,'Target keyword not compatible'
   retall
endif
if(strcmp(instrument1,instrument2) eq 0) then begin
   massage,'Instrument keyword not compatible'
   massage,'continue at your own risk'
;   retall
endif

tol_t = 1.0d-6

if(abs(T1-T2) gt tol_t) then begin
   print,T1,T2,abs(T1-T2),tol_t
   massage,'Time resolution not compatible'
   retall
endif


;if(abs(T1-T2) ne 0.0) then begin
;   massage,'Time resolution not compatible'
;   print,T1,T2
;   retall
;endif

if(nft1 ne nft2) then begin
   massage,'Frequency resolution not compatible'
   retall
endif
if(abs(rmjd01-rmjd02) ne 0.0) then begin
   massage,'Start time not compatible'
   retall
endif
if(ntrafos1 ne ntrafos2) then begin
   massage,'Number of FFT not compatible'
   retall
endif
; MM
;
; MM 16 May 2026
; Zhang Poisson correction is valid only when both files are PCA.
;
is_pca = (strcmp(instrument1,'PCA             ') eq 1) and $
         (strcmp(instrument2,'PCA             ') eq 1)
; MM

nft     = nft1/2
ntrafos = ntrafos1
T       = T1
;
rdata1    = complexarr(nft)
rdata2    = complexarr(nft)
ccsumjk   = complexarr(nft)

; MM 15 May 2026
; Accumulate second moments of the unrebinned cross vectors.
; These are used later to estimate the covariance between Re(C) and Im(C)
; after frequency rebinning. This is needed for a proper propagated error
; on |C| = sqrt(Re(C)^2 + Im(C)^2).
rrsumjk = dblarr(nft)*0.0d0
iisumjk = dblarr(nft)*0.0d0
risumjk = dblarr(nft)*0.0d0
; MM

coherence = ccsumjk
coherence_err = coherence

power1     = fltarr(nft)*0.0
power1_err = fltarr(nft)*0.0
power2     = fltarr(nft)*0.0
power2_err = fltarr(nft)*0.0
;
; MM
; Define the vectors to hole the Re/Im parts of the CDS and corresponding errors.
; The Re' and Im' are the c omponents of the cross vector when it is rotated by X degrees.
; This may be needed to make the new Real/Imag parts more or less of the same magnitude
; cause otherwise the fits have problems.
; Notice that a rotation of the cross vector does not affect the modulus of the vector, it only
; adds a constant to the lags, which is harmless because the reference band of the lags is
; arbitrary
;
; Real part of CDS to pass back
repart     = fltarr(nft)*0.0
; (Real part)' after rotating the cross vector by X degrees
repart_prime = fltarr(nft)*0.0
repart_err = fltarr(nft)*0.0
; Imag part of CDS to pass back
impart     = fltarr(nft)*0.0
; (Imag part)' after rotating the cross vector by X degrees
impart_prime = fltarr(nft)*0.0
impart_err = fltarr(nft)*0.0
; MM

;
; Loop over the FFTs
;
n_selected       = 0l
itrafos          = 0l
flux             = 0.0d0
; MM
flux1            = 0.0d0
flux2            = 0.0d0
; MM
tdead            = 1.0e-5
poilevel1        = fltarr(nft)*0.0
poilevel2        = fltarr(nft)*0.0

; --------------- Data selection ------------------

if(keyword_set(index)) then begin
    firstpds = index[0]
    lastpds  = index[1]
endif else begin
    firstpds = 0L
    lastpds  = 100000L
endelse
if(keyword_set(time)) then begin
    firsttime = time[0]
    lasttime  = time[1]
endif else begin
    firsttime = 0.0D0
    lasttime  = 1.0e10
endelse
if(keyword_set(rate)) then begin
    firstrate = rate[0]
    lastrate  = rate[1]
endif else begin
    firstrate = 0.0
    lastrate  = 2000000.0
endelse
if(keyword_set(sel)) then begin
    goodarray = sel  ; used to be sel-1
endif else begin
    goodarray = indgen(ntrafos)
endelse

goodarray = goodarray(where((goodarray ge firstpds) and (goodarray le lastpds)))
lastindex = max(goodarray)
a0 = 0.0
rmjd1 = 0d0
rmjd2 = 0d0

;
;  Now fill in frequency array
;
frequency = (findgen(nft)+1.0) * T

; MM
;
; POISSON range bookkeeping.
; The last Fourier bin is the Nyquist bin. It has different statistics
; for a real time series and should not be used to estimate the noise
; level or Real-part correction.
;
n_frequency = n_elements(frequency)
last_non_nyq = n_frequency - 2L

if(last_non_nyq LT 0L) then begin
   massage,'Frequency array too short to exclude Nyquist bin'
   retall
endif

;
; Check MANUAL_POI.
;
if(keyword_set(manual_poi)) then begin
   if(n_elements(manual_poi) ne 2) then begin
      massage,'MANUAL_POI must have two elements: [pds1,pds2]'
      retall
   endif
endif

;
; If POISSON was requested, build a safe index array once.
; We clip the requested range to the available non-Nyquist frequencies.
;
if(keyword_set(poisson)) then begin

   if(n_elements(poisson) ne 2) then begin
      massage,'POISSON must have two elements: [f1,f2]'
      retall
   endif

   fpoi1_user = double(poisson[0])
   fpoi2_user = double(poisson[1])

   if(fpoi2_user lt fpoi1_user) then begin
      massage,'POISSON range invalid: f2 < f1'
      retall
   endif

   fmin_valid = double(frequency[0])
   fmax_valid = double(frequency[last_non_nyq])

   fpoi1 = fpoi1_user > fmin_valid
   fpoi2 = fpoi2_user < fmax_valid

   if(fpoi2 lt fpoi1) then begin
      massage,'POISSON range has no overlap with available non-Nyquist frequencies'
      retall
   endif

   if((fpoi1 ne fpoi1_user) or (fpoi2 ne fpoi2_user)) then begin
      print,'WARNING: POISSON range clipped to available non-Nyquist frequencies.'
      print,'         Requested: ',fpoi1_user,' - ',fpoi2_user,' Hz'
      print,'         Used     : ',fpoi1,' - ',fpoi2,' Hz'
   endif

   wpoi = where((frequency ge fpoi1) and $
                (frequency le fpoi2) and $
                (indgen(n_frequency) le last_non_nyq), npoi)

   if(npoi le 0) then begin
      massage,'POISSON range has no valid non-Nyquist frequency bins'
      retall
   endif

endif
; MM

printed_zhang_message = 0

for itrafos=0l,ntrafos-1l do begin
;
read_fft_line,unit1,mufla1,rmjd1,cnts1,poisson1,current_vle_rate1,fndet1,rdata1
   read_fft_line,unit2,mufla2,rmjd2,cnts2,poisson2,current_vle_rate2,fndet2,rdata2
   cnts = cnts1+cnts2
   timediff = abs(rmjd1-rmjd2)
   if(timediff gt 1.0e-4) then begin
       massage,'Times of FFTs not compatible!'
       retall
   endif   

   t_1 = ((rmjd1-rmjd01)*86400.0)
   t_2 = t_1 + 1.0/T1
   gotcha = where(goodarray eq itrafos)
;
;  Accumulate sum and sum of squares
;
  if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (gotcha[0] ge 0)) then begin  
;  Here JK -----------------------------
      cc       = conj(rdata1) * rdata2 * 2.0 /sqrt(cnts1*cnts2)
      ccsumjk  = ccsumjk + cc       ; it remains a complex number
; MM 15 May 2026
; Second moments of individual cross vectors before averaging.
; These are accumulated before frequency rebinning. Later, the same
; frequency rebinning is applied to the moments.
rrsumjk = rrsumjk + double(float(cc))^2
iisumjk = iisumjk + double(imaginary(cc))^2
risumjk = risumjk + double(float(cc)) * double(imaginary(cc))
; MM
;--Here powers--------------------------
      pwr1      = abs(rdata1)^2 * 2.0/cnts1
      power1   = power1 + pwr1
      pwr2      = abs(rdata2)^2 * 2.0/cnts2
      power2   = power2 + pwr2
	 
      n_selected = n_selected + 1
      flux      = flux + cnts * T
; MM
      flux1     = flux1 + cnts1 * T
      flux2     = flux2 + cnts2 * T
; MM
     
; MM
; Here Poissonian/noise estimates for the two bands.
;
; Priority:
;   1. MANUAL_POI=[p1,p2]
;   2. POISSON=[f1,f2]
;   3. RAWCOH: no subtraction, even for RXTE/PCA
;   4. RXTE/PCA Zhang correction
;   5. no subtraction, raw coherence
;
      if(keyword_set(manual_poi)) then begin

         ; User gives fixed PDS noise levels by hand.
         ; These affect the intrinsic-coherence correction only.
         poilevel1 = float(manual_poi[0]) + frequency*0.0
         poilevel2 = float(manual_poi[1]) + frequency*0.0

      endif else begin

         if(keyword_set(poisson)) then begin

            ; Estimate PDS noise levels from the requested frequency range,
            ; excluding the Nyquist bin.
            poilevel1 = mean(pwr1[wpoi]) + frequency*0.0
            poilevel2 = mean(pwr2[wpoi]) + frequency*0.0

         endif else begin

            if(keyword_set(rawcoh)) then begin

               ; Explicit measured/raw coherence: use observed PDS in the
               ; denominator and do not apply automatic PCA/Zhang correction.
               poilevel1 = frequency*0.0
               poilevel2 = frequency*0.0

            endif else begin

               if(is_pca) then begin
                  ; BUG FIX:
                  ; Accumulate the Zhang Poisson estimate over selected FFTs.
                  ; The old version overwrote POILEVEL1/2 for each FFT and
                  ; later divided only the last value by N_SELECTED.
                  poisson_estimate,poitmp1,cnts1,1.0/T,current_vle_rate1,nft,fndet1,i_vle,tdead
                  poisson_estimate,poitmp2,cnts2,1.0/T,current_vle_rate2,nft,fndet2,i_vle,tdead

                  poilevel1 = poilevel1 + poitmp1
                  poilevel2 = poilevel2 + poitmp2

                  if(printed_zhang_message eq 0) then begin
                     if(~keyword_set(status_quiet)) then begin
                        massage,'This is RXTE/PCA: Poisson level computed from Zhang (1995)'
                     endif
                     printed_zhang_message = 1
                  endif

               endif else begin

                  ; No noise subtraction: raw coherence.
                  poilevel1 = frequency*0.0
                  poilevel2 = frequency*0.0

               endelse

            endelse

         endelse

      endelse
; MM

   endif
   
   if(itrafos gt (lastindex-1)) then goto,finito
;
endfor

finito:

;
;  Compute average and standard deviation
;
   if(n_selected gt 0) then begin
      flux      = flux      / n_selected
; MM
      flux1     = flux1     / n_selected
      flux2     = flux2     / n_selected
; MM
      poilevel1 = poilevel1 / n_selected
      poilevel2 = poilevel2 / n_selected
      power1    = power1    / n_selected
      power2    = power2    / n_selected
      power1_err= power1    / sqrt(n_selected)
      power2_err= power2    / sqrt(n_selected)
	  
;     compute again poissonian level, but from the average
; MM      IF(keyword_set(poisson)) THEN BEGIN
; MM        massage,'Poisson level recomputed from average'
; MM        p1 = where(frequency ge poisson[0])
; MM        p1 = p1[0]
; MM        p2 = where(frequency le poisson[1])
; MM        p2 = p2(n_elements(p2)-1)
; MM        poilevel1 = mean(power1[p1:p2])+frequency*0
; MM        poilevel2 = mean(power2[p1:p2])+frequency*0
; MM	  ENDIF

;MM
      if(keyword_set(manual_poi)) then begin

         ; Keep the manually supplied PDS noise levels.
         poilevel1 = float(manual_poi[0]) + frequency*0.0
         poilevel2 = float(manual_poi[1]) + frequency*0.0

      endif else begin

         if(keyword_set(poisson)) then begin
            massage,'Poisson level recomputed from average'
            poilevel1 = mean(power1[wpoi]) + frequency*0.0
            poilevel2 = mean(power2[wpoi]) + frequency*0.0
         endif

      endelse
;MM
     
        power1    = power1 - poilevel1 
        power2    = power2 - poilevel2
      
        ccsumjk   = ccsumjk   / n_selected
; MM 15 May 2026
; Convert accumulated second moments to per-frequency averages over FFTs.
; Do not form the covariance yet, because frequency rebinning still has to
; combine neighbouring Fourier bins. The covariance must be computed after
; the final binning.
        rrsumjk = rrsumjk / double(n_selected)
        iisumjk = iisumjk / double(n_selected)
        risumjk = risumjk / double(n_selected)
; MM
      
     endif else begin
      massage,'Warning! No FFTs retrieved for input selection!'
      retall
   endelse
free_lun,unit1
free_lun,unit2

;
;  Now fill in frequency array
;
   frequency = (findgen(nft)+1.0) * T
;
   ps = ' FFTs '
   if(n_selected eq 1) then ps = ' FFT '
   if(~keyword_set(status_quiet)) then begin
      print,'  ',strtrim(string(n_selected),1),' '+ps+'selected'
   endif
;
;  Now rebinning in frequency
;
;    Here JK ------------------------------------------------
   cr      = float(ccsumjk)
   ci      = imaginary(ccsumjk)
   rebincross,frequency,cr,irf,frequency_reb,cr_reb,nn
   rebincross,frequency,ci,irf,frequency_reb,ci_reb,nn
   ccsumjk_reb  = complex(cr_reb,ci_reb)
; MM 15 May 2026
; Rebin the second moments in frequency with the same binning as Re and Im.
; r2_reb, i2_reb, and ri_reb are the means of R^2, I^2, and R*I over all
; samples contributing to each rebinned frequency bin.
   rebincross,frequency,float(rrsumjk),irf,zzz,r2_reb,nn2
   rebincross,frequency,float(iisumjk),irf,zzz,i2_reb,nn2
   rebincross,frequency,float(risumjk),irf,zzz,ri_reb,nn2
; MM
;------------------------------------------------------------
;    Here the two PDS ---------------------------------------

   gh_reb,frequency,power1,power1_err,irf,x1,y1,y1e
   gh_reb,frequency,power2,power2_err,irf,x2,y2,y2e
;------------------------------------------------------------
;    Here the poissonian levels, which will be used later
   gh_reb,frequency,double(poilevel1),poilevel1*0,irf,zzz,nslev1,www
   gh_reb,frequency,double(poilevel2),poilevel2*0,irf,zzz,nslev2,www
   
; Factor MW
   mw = float(nn * n_selected)
   
; Computation of lags and coherence (from Phil's formulas)

; Calculate noise correction to coherence (Vaughan & Nowak 1997, see paragraph
; following Eqn. 4)

nsq     = ((y1*nslev2)+(y2*nslev1)+(nslev1*nslev2))/mw

; Calculate the coherence corrected for noise (VN97 eqn 8)

csq     = abs(ccsumjk_reb)^2
coh     = (csq - nsq)/(y1*y2) 

; Calculate the error on coherence (VN97 eqn 8)
coherr1 = (2.0*(nsq^2)*mw)/(csq-nsq)^2
coherr2 = (nslev1/y1)^2 + (nslev2/y2)^2
coherr3 = 2.0*(1-coh)^2/coh ; Corrected by FG reinterpreting VN97 definition
; coherr3 = 2.0*((1-coh)/(coh^1.5))^2 ; TMB original
coherr  = coh/sqrt(mw) * sqrt(coherr1+coherr2+coherr3)   


; FG 22 Mar 2024 - Added "lowcoh" option for Case II in VN97 (eqn 9).
IF (keyword_set(lowcoh)) THEN BEGIN
    ; When power < nslev, y<0, thus we take "absolute values" here.
    nsq     = (abs(y1*nslev2)+abs(y2*nslev1)+abs(nslev1*nslev2))/mw
    csq     = abs(ccsumjk_reb)^2
    coh     = abs(csq - nsq)/abs(y1*y2)

    nfreqs = N_ELEMENTS(csq)
    qf = FLTARR(nfreqs)*0
    ssq = FLTARR(nfreqs)
    coherr1 = FLTARR(nfreqs)        

    FOR i = 0, nfreqs-1 DO BEGIN
        ; Pre-define channels marked as bad for gh_xspec with NaN.
        qf[i] = 1
        coherr1[i] = sqrt(-1)          
        ; Check wether we are in the "useful" range to estimate instrinsic coherence
        IF ((csq[i] GT 1*nsq[i]) && $
            (y1[i]/nslev1[i] GT 1/sqrt(mw[i])) && $
            (y2[i]/nslev2[i] GT 1/sqrt(mw[i]))) THEN BEGIN
            ; Check wether we are in the Gaussian limit (VN97 eqn 8)
            IF ((coh[i] GT 10*nsq[i]/((y1[i]+nslev1[i])*(y2[i]+nslev2[i]))) && $
                (y1[i] GT 10*nslev1[i]/sqrt(mw[i])) && $
                (y2[i] GT 10*nslev2[i]/sqrt(mw[i]))) THEN BEGIN
;                  coherr1[i] = 2.0*nsq[i]^2*mw[i]/(csq[i]-nsq[i])^2 ; VN97 eqn 8
;                  coherr1[i] = 2.0*nsq[i]*mw[i]/(csq[i]-nsq[i])
                  coherr1[i] = 2.0*nsq[i]/(csq[i]-nsq[i]) ; FG correction to VN97 eqn 8
                  qf[i] = 0
            ; Or otherwise use the low measured coherence case (VN97 eqn 9)
            ENDIF ELSE BEGIN
                  ssq[i] = gh_lowcoh_ssq(csq[i], nsq[i])
                  coh[i] = ssq[i]/abs(y1[i]*y2[i])
                  CI = gh_lowcoh_ssq_CI(0.68, csq[i], nsq[i])
;                  coherr1[i] = mw[i]*(0.5*(CI[1]-CI[0])/ssq[i])^2
                  coherr1[i] = (0.5*(CI[1]-CI[0])/ssq[i])^2 ; FG error propagation
                  qf[i] = 2
            ENDELSE
        ENDIF
    ENDFOR
 ; Calculate the full error on the coherence (NV97 eqns 8-9, typo check)
    coherr2 = (nslev1/y1)^2 + (nslev2/y2)^2
    coherr3 = 2.0*(1-coh)^2/coh ; Corrected by FG reinterpreting notation in VN97
    coherr  = coh/sqrt(mw)*sqrt(coherr1+coherr2+coherr3)        
ENDIF
; FG 22 Mar 2024


; If an average real part of the noise cross-spectrum (evaluated in a previous run of gh_cross_range)
; is set, subtract it before computing the lag

 if(keyword_set(poivalue)) then begin
      ccsumjk_reb = ccsumjk_reb - complex(poivalue,0)
 endif

; MM 15 May 2026
; Recompute Re, Im, and modulus after optional subtraction of the average
; noise/crosstalk real component. The covariance matrix is unchanged by
; subtracting a constant from Re, but the propagation direction changes
; because the mean vector has changed.
repart = float(ccsumjk_reb)
impart = imaginary(ccsumjk_reb)
modcs  = sqrt(repart^2 + impart^2)

; Initialise before assigning into modcs_err[wmod]
modcs_err = fltarr(n_elements(modcs))*0.0

; Number of independent cross vectors contributing to each rebinned bin.
; This is the same MW used elsewhere in the routine.
mw = float(nn * n_selected)

; Empirical covariance matrix of the mean cross vector.
; r2_reb, i2_reb, ri_reb are averages of R^2, I^2, R*I over all samples.
; Therefore:
;     Var(mean R)   = ( <R^2>  - <R>^2 ) / (N-1)
;     Var(mean I)   = ( <I^2>  - <I>^2 ) / (N-1)
;     Cov(mean R,I) = ( <RI>   - <R><I> ) / (N-1)
; with N = MW. This is equivalent to sample covariance divided by N.
var_re_mean  = fltarr(n_elements(modcs))*0.0
var_im_mean  = fltarr(n_elements(modcs))*0.0
cov_reim_mean = fltarr(n_elements(modcs))*0.0

wgoodmw = where(mw gt 1.0, ngoodmw)
if(ngoodmw gt 0) then begin
   var_re_mean[wgoodmw]   = (r2_reb[wgoodmw] - repart[wgoodmw]^2) / (mw[wgoodmw]-1.0)
   var_im_mean[wgoodmw]   = (i2_reb[wgoodmw] - impart[wgoodmw]^2) / (mw[wgoodmw]-1.0)
   cov_reim_mean[wgoodmw] = (ri_reb[wgoodmw] - repart[wgoodmw]*impart[wgoodmw]) / (mw[wgoodmw]-1.0)

   ; Numerical protection. Small negative variances can appear from roundoff.
   badv = where(var_re_mean lt 0.0, nbadv)
   if(nbadv gt 0) then var_re_mean[badv] = 0.0

   badv = where(var_im_mean lt 0.0, nbadv)
   if(nbadv gt 0) then var_im_mean[badv] = 0.0
endif

; Propagate the full 2D covariance matrix to |C|.
; This includes the covariance term that is missing from the simpler
; quadrature expression.
wmod = where(modcs gt 0.0, nmod)
if(nmod gt 0) then begin
   modvar = $
      (repart[wmod]/modcs[wmod])^2 * var_re_mean[wmod] + $
      (impart[wmod]/modcs[wmod])^2 * var_im_mean[wmod] + $
      2.0*(repart[wmod]/modcs[wmod])*(impart[wmod]/modcs[wmod]) * cov_reim_mean[wmod]

   ; Numerical protection against tiny negative values from roundoff.
   badm = where(modvar lt 0.0, nbadm)
   if(nbadm gt 0) then modvar[badm] = 0.0

   modcs_err[wmod] = sqrt(modvar)
endif

; MM

; Phase and phase error (from Nowak et al. 1999). Note that the phase error
; use coherence uncorrected for noise levels, so we add the noise back on to 
; the powers and we use csq in the numerator since we do not want to subtract nsq

lag       = atan(imaginary(ccsumjk_reb),float(ccsumjk_reb))
phaserr1  = csq/((y1+nslev1)*(y2+nslev2))
lag_err   = sqrt((1.0-phaserr1)/(2.0*phaserr1*mw))

; MM

; Re/Im here are in Leahy units
; Re/Im here are in Leahy units. MODCS and MODCS_ERR were already
; computed above, after the optional POIVALUE subtraction.

repart    = float(ccsumjk_reb)
impart    = imaginary(ccsumjk_reb)

; Rotate the cross vector by X=45 degrees.
; The modulus and its covariance-propagated error are not changed by this
; rotation, so MODCS and MODCS_ERR are left untouched.

   if(keyword_set(rotate)) then begin
     theta=!DPI/4.0
     repart_prime = repart*cos(theta) - impart*sin(theta)
     impart_prime = repart*sin(theta) + impart*cos(theta)
     
     repart = repart_prime
     impart = impart_prime

   endif


; add noise back to the powers to get the errors of the Re and Im parts

; Errors are taken from eq. 13 in Ingram 2019. Here errors are in Leahy units
; (notice that for the moment I assumed that the errors are independent of the rotation the cross vector.
; I need to check this)
repart_err= sqrt(((y1 + nslev1)*(y2 + nslev2) + repart*repart - impart*impart)/(2.0*mw))
impart_err= sqrt(((y1 + nslev1)*(y2 + nslev2) - repart*repart + impart*impart)/(2.0*mw))



; check that is a background is given, both are given.
   if((keyword_set(back1)) and ~(keyword_set(back2))) then begin
         massage,'Error: if you give the background you must give it in both bands! (use rms1=back1,rms2=back2)'
         retall
   endif

   if((keyword_set(back2)) and ~(keyword_set(back1))) then begin
         massage,'Error: if you give the background you must give it in both bands! (use rms1=back1,rms2=back2)'
         retall
   endif

;
;  Rms^2 conversion
;

   if(keyword_set(back1)) then begin
;      print,'Rms normalisation'

      if(back1 ge cnts1) then begin
         massage,'Error: background flux higher than source+bkg flux in band1'
;       retall
      endif
      if(back2 ge cnts2) then begin
         massage,'Error: background flux higher than source+bkg flux in band2!'
;         retall
      endif

; This part is experimental.
; Convert the Re/Im to rms^2.
; For this multiply the Re/Im parts and errors by sqrt(C1*C2)/(S1*S2)
; where: 
; S1 = C1 - B1 
; S2 = C2 - B2
; with:
; C1 = total counts (source+background) in band 1
; B1 = background in band 1
; C2 = total counts (source+background) in band 2
; B2 = background in band 2



;      norm1         = sqrt(cnts1 * cnts2)/((cnts1 - back1) * (cnts2 - back2))
      norm_flux     = sqrt(flux1 * flux2)/((flux1 - back1) * (flux2 - back2))
; This norm is not working yet
;      massage,'The rms norm is not working yet. Use Leahy'
;      norm1         = sqrt(sqrt(cnts1*cnts2)/((cnts1 - back1) * (cnts2 - back2)))
;      norm1         = 0.00065
;      print,'cnts: ',cnts,' cnts1: ',cnts1,' cnts2: ',cnts2
;      print,'T: ',T
;      print,'norm_flux: ',norm_flux
;      print,'flux1=cnts1*T: ',flux1,' flux2=cnts2*T: ',flux2,' norm_flux: ',norm_flux
;      print,'norm1: ',norm1
      repart        = repart    * norm_flux
      impart        = impart    * norm_flux
      repart_err    = repart_err * norm_flux
      impart_err    = impart_err * norm_flux
; MM 15 May 2026
; The modulus has the same units as Re and Im, so apply the same conversion.
      modcs     = modcs     * norm_flux
      modcs_err = modcs_err * norm_flux
; MM

   endif

; MM

; Here integrate the lags in the Poisson range, if specified

; MM
IF(keyword_set(poisson)) THEN BEGIN

   ;
   ; Apply the Real-part correction using the same requested POISSON range,
   ; but now on the rebinned frequency grid. Exclude the rebinned Nyquist bin
   ; if it is present as the last element.
   ;
   nfreb = n_elements(frequency_reb)
   last_reb_non_nyq = nfreb - 2L

   if(last_reb_non_nyq LT 0L) then begin
      massage,'Rebinned frequency array too short for POISSON Real correction'
      retall
   endif

   fmin_reb = double(frequency_reb[0])
   fmax_reb = double(frequency_reb[last_reb_non_nyq])

   fpoi1_reb = double(poisson[0]) > fmin_reb
   fpoi2_reb = double(poisson[1]) < fmax_reb

   if(fpoi2_reb lt fpoi1_reb) then begin
      massage,'POISSON range has no overlap with rebinned non-Nyquist frequencies'
      retall
   endif

   wpoir = where((frequency_reb ge fpoi1_reb) and $
                 (frequency_reb le fpoi2_reb) and $
                 (indgen(nfreb) le last_reb_non_nyq), npoir)

   if(npoir le 0) then begin
      massage,'POISSON range has no valid rebinned bins for Real correction'
      retall
   endif

   noisereal = mean(float(ccsumjk_reb[wpoir]))
   noiseimag = mean(imaginary(ccsumjk_reb[wpoir]))

   lag = atan(imaginary(ccsumjk_reb)-noiseimag, $
              float(ccsumjk_reb)-noisereal)

ENDIF

; MM


; MM IF(keyword_set(poisson)) THEN BEGIN
; MM      p1 = where(frequency_reb ge poisson[0])
; MM      p1 = p1[0]
; MM      p2 = where(frequency_reb le poisson[1])
; MM      p2 = p2(n_elements(p2)-1)
; MM 	noisereal = mean(float(ccsumjk_reb[p1:p2]))
; MM 	noiseimag = mean(imaginary(ccsumjk_reb[p1:p2]))
; MM 	lag       = atan(imaginary(ccsumjk_reb)-noiseimag,float(ccsumjk_reb)-noisereal)
; MM ENDIF

coherence     = coh
coherence_err = coherr
;
;
;  Switch for time lags
;
   if(keyword_set(tl)) then begin
      lag     = lag/2.0/!PI/frequency_reb
      lag_err = lag_err/2.0/!PI/frequency_reb
   endif

; MM
;
; Tell the user what noise prescription was used for coherence.
;

if(~keyword_set(status_quiet)) then begin

if(keyword_set(manual_poi)) then begin

   print,'Using MANUAL_POI for intrinsic coherence.'
   print,'PDS1/PDS2 levels: ',manual_poi[0],manual_poi[1]
   print,'No Real-part correction is applied by MANUAL_POI.'
   print,'Assuming values are in input-file PDS units, before RMS1/RMS2 conversion.'

   if(keyword_set(back1) or keyword_set(back2)) then begin
      print,'If the constants were fitted from RMS-normalised PDS,'
      print,'rerun GHX without RMS to measure them in input-file units:'
      print,"   ghx,'file1.pds',nu,pow,powe"
      print,"   ghx,'file2.pds',nu,pow,powe"
   endif

endif else begin

   if(keyword_set(poisson)) then begin

      print,'Using POISSON range for intrinsic coherence.'
      print,'Frequency range used for PDS + Real correction: ', $
            fpoi1_user,' - ',fpoi2_user,' Hz'

   endif else begin

      if(keyword_set(rawcoh)) then begin
         print,'RAWCOH: raw/measured coherence returned.'
         print,'Observed PDS used in the denominator; no Poisson/noise correction.'
         print,'Coherence errors use nslev1=nslev2=0.'
      endif else begin

         if(is_pca) then begin
            print,'PCA data: using Zhang (1995) Poisson correction.'
            print,'Use POISSON=[f1,f2] or MANUAL_POI=[p1,p2] to override.'
            print,'POISSON: f1-f2 used for PDS + Real correction.'
            print,'MANUAL_POI: fixed PDS1/PDS2 levels only.'
         endif else begin
            print,'WARNING: Raw coherence returned.'
            print,'Use POISSON=[f1,f2] or MANUAL_POI=[p1,p2]'
            print,'for intrinsic coherence.'
            print,'POISSON: f1-f2 used for PDS + Real correction.'
            print,'MANUAL_POI: fixed PDS1/PDS2 levels only.'
         endelse

      endelse

   endelse

endelse

endif

; MM

end
