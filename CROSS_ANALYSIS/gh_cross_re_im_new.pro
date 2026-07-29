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
                 repart,repart_err,impart,impart_err,mw,qf,$
	         index=index,time=time,sel=sel,rate=rate,poisson=poisson,poivalue=poivalue, $
                 rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh
;+
; NAME: 
;      GH_CROSS
; PURPOSE: 
;      Reads is two FFTs from FFT files, selecting on a range of FFTs
;			in time, range, or selection of indices
;      This procedure extracts two FFTs from two compatible files. The
;      selection is possible only on a continuous range of FFTs.
;      Produces phase lags (in radians) or time lags, and coherence.
;
; CALLING SEQUENCE: 
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
;       Notice that, except for RXTE, in this case the intrinsic coherence is computed assuming that Poisson=2
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
;			
;
;--------------------------------------------------------------------------
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
   massage,'Time resolution not compatible'
   print,T1,T2
   retall
endif
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
nft     = nft1/2
ntrafos = ntrafos1
T       = T1
;
rdata1    = complexarr(nft)
rdata2    = complexarr(nft)
ccsumjk   = complexarr(nft)
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

IF(keyword_set(poisson)) THEN BEGIN
   n_frequency = n_elements(frequency)
   last_non_nyq = n_frequency - 2L
   if(last_non_nyq LT 0L) then begin
      massage,'Poisson frequency range outside frequency array'
      retall
   endif
   fpoi1 = double(poisson[0]) > double(frequency[0])
   fpoi2 = double(poisson[1]) < double(frequency[last_non_nyq])
   if(fpoi2 lt fpoi1) then begin
      massage,'Poisson frequency range outside frequency array'
      retall
   endif
   wpoi = where((frequency ge fpoi1) and $
                (frequency le fpoi2) and $
                (lindgen(n_frequency) le last_non_nyq), npoi)
   if(npoi le 0) then begin
      massage,'Poisson frequency range outside frequency array'
      retall
   endif
ENDIF

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
     
; Here Poissonian estimates for the two bands
      if(strcmp(instrument1,'PCA             ') eq 1) then begin
          poisson_estimate,poilevel1,cnts1,1.0/T,current_vle_rate1,nft,fndet1,i_vle,tdead
          poisson_estimate,poilevel2,cnts2,1.0/T,current_vle_rate2,nft,fndet2,i_vle,tdead
          massage,'This is RXTE/PCA: Poisson level computed from formulae in Zhang (1995)'
      endif else begin
		  IF(keyword_set(poisson)) THEN BEGIN
			  p1 = where(frequency ge poisson[0])
			  p1 = p1[0]
			  p2 = where(frequency le poisson[1])
			  p2 = p2(n_elements(p2)-1)
			  poilevel1 = mean(pwr1[wpoi])+frequency*0
			  poilevel2 = mean(pwr2[wpoi])+frequency*0
                          ;massage,'Poisson range given'
		  ENDIF ELSE BEGIN
	         poilevel1 = 2.0*frequency
	         poilevel2 = 2.0*frequency
			 ;massage,'WARNING: assuming poissonian level is 2.0!'
	      ENDELSE
      endelse
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
      IF(keyword_set(poisson)) THEN BEGIN
        massage,'Poisson level recomputed from average'
        poilevel1 = mean(power1[wpoi])+frequency*0
        poilevel2 = mean(power2[wpoi])+frequency*0
	  ENDIF
     
        power1    = power1 - poilevel1 
        power2    = power2 - poilevel2
      
        ccsumjk   = ccsumjk   / n_selected
      
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
   print,'  ',strtrim(string(n_selected),1),' '+ps+'selected'
;
;  Now rebinning in frequency
;
;    Here JK ------------------------------------------------
   cr      = float(ccsumjk)
   ci      = imaginary(ccsumjk)
   rebincross,frequency,cr,irf,frequency_reb,cr_reb,nn
   rebincross,frequency,ci,irf,frequency_reb,ci_reb,nn
   ccsumjk_reb  = complex(cr_reb,ci_reb)
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

; Phase and phase error (from Nowak et al. 1999). Note that the phase error
; use coherence uncorrected for noise levels, so we add the noise back on to 
; the powers and we use csq in the numerator since we do not want to subtract nsq

lag       = atan(imaginary(ccsumjk_reb),float(ccsumjk_reb))
phaserr1  = csq/((y1+nslev1)*(y2+nslev2))
lag_err   = sqrt((1.0-phaserr1)/(2.0*phaserr1*mw))

; MM

; Re/Im here are in Leahy units
repart    = float(ccsumjk_reb)
impart    = imaginary(ccsumjk_reb)

; Rotate the cross vector by X=45 degrees (in the future maybe allow the user to choose the angle)

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

   endif

; MM

; Here integrate the lags in the Poisson range, if specified

IF(keyword_set(poisson)) THEN BEGIN
     nfreb = n_elements(frequency_reb)
     last_reb_non_nyq = nfreb - 2L
     if(last_reb_non_nyq LT 0L) then begin
        massage,'Poisson frequency range outside frequency array'
        retall
     endif
     fpoi1_reb = double(poisson[0]) > double(frequency_reb[0])
     fpoi2_reb = double(poisson[1]) < double(frequency_reb[last_reb_non_nyq])
     if(fpoi2_reb lt fpoi1_reb) then begin
        massage,'Poisson frequency range outside frequency array'
        retall
     endif
     wpoir = where((frequency_reb ge fpoi1_reb) and $
                   (frequency_reb le fpoi2_reb) and $
                   (lindgen(nfreb) le last_reb_non_nyq), npoir)
     if(npoir le 0) then begin
        massage,'Poisson frequency range outside frequency array'
        retall
     endif
	noisereal = mean(float(ccsumjk_reb[wpoir]))
	noiseimag = mean(imaginary(ccsumjk_reb[wpoir]))
	lag       = atan(imaginary(ccsumjk_reb)-noiseimag,float(ccsumjk_reb)-noisereal)
ENDIF

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

end
