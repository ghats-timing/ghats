pro gh_cs,filename1,filename2,             $
             frequency_reb,lag,lag_err,coherence,coherence_err, $
             irf,tlag=tl, $
             repart,repart_err,impart,impart_err,modcs,modcs_err,mw,qf, $
             index=index,time=time,sel=sel,rate=rate,poisson=poisson, $
             manual_poi=manual_poi, $
             rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh,rawcoh=rawcoh,help=help
;+
; NAME: 
;      GH_CS
; PURPOSE: 
;      This is a wrapper that takes the user's input and runs gh_cross_re_im
;      with proper parameters, the necessary number of times to compute and
;      subtract the Real part of the Poisson-dominated part of the cross spectrum
;      to correct for cross talk and for partial correlation if one uses the full 
;      band as reference.
;
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
;       MODCS        = Modulus of the cross spectrum, sqrt(REPART^2+IMPART^2)
;       MODCS_ERR    = Error on MODCS, propagated including empirical Re-Im covariance
;
; KEYWORDS:
;       TIME         = Switch for time lags (default phase lags)
;       POISSON      = This is an array that has start/end freq. to 
;                      compute the Poisson level; for instance 'poisson=[800,1200]'
;       RMS1         = background value for rms2 conversion for band1. If specified, 
;                      rms2 conversion will be performed
;       RMS2         = background value for rms2 conversion for band2. If specified, 
;                      rms2 conversion will be performed
;       ROTATE       = Switch to rotate the cross vector by X (currently pi/4) rad
;                      The user should remember to subtract X rad from the phase lags.
;                      Default is no rotate
;       LOWCOH       = Switch to use Eq 9 from VN97 to calculate the Coherence in the
;                      High-Power Low-Coherence limit.
;       RAWCOH       = Switch to return measured/raw coherence, using the observed
;                      PDS in the denominator and no Poisson/noise correction.
;
; EXAMPLE:
;       First read in two FFTs (R=Reference,band  S=Subject band)
;
;       Ghats> filer='fileR.fft' & files='fileS.fft'
;       Ghats> gh_cs,filer,files,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,poisson=[400,800],rms1=0.001,rms2=0.001,/rotate
;       Ghats> gh_cs,filer,files,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,poisson=[400,800],rms1=0.001,rms2=0.001,/rotate,/lowcoh
;
;       This will find the average of the Real part of the cross spectrum over a frequency range dominated by Poisson noise (in this case 400-800 Hz)
;       and subtract that average real value from the Real part of the cross spectrum
;       (This corrects both for cross talk and for the partial correlation if the subject band is part of te reference band)
;
;       After that, for instance, write the Real and Imaginary parts of the cross spectrum, the lags and the coherence function in Xspec
;
;       Ghats> gh_xspec,freq,real,real_err,'real'
;       Ghats> gh_xspec,freq,imag,imag_err,'imag'
;       Ghats> gh_xspec,freq,lag,lag_err,'lag'
;       Ghats> gh_xspec,freq,coh,coh_err,'coh'
;       Ghats> gh_xspec,freq,mod,mod_err,'mod'
;
;       NOTES:
;       The keyword poisson=[freq1,freq2] indicates 2 things:
;       1. Correct the cross spectrum by cross talk/partial correlation (for lags)
;       2. Subtract the Poisson level from the PDS in the 2 individual bands (for coherence function)
;       If the keyword is given, the same frequency range (freq1-freq2) is used for both corrections.

;       For RXTE the intrinsic coherence is computed assuming Zhang 1995 formulae for the Poisson level.
;       If no noise keyword is given, for any other mission the program returns raw/measured
;       coherence with the observed PDS in the denominator. /RAWCOH forces this same raw
;       coherence definition also for RXTE/PCA data.


;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       GH_CROSS_RE_IM     Does all the dirty work
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       M. Mendez 13 Aug 2023 from GH_CROSS_RE_IM
;       F. Garcia 22 Mar 2024 added LOWCOH variant for VN97 Eq 9
;       M. Mendez 15 May 2026 added /HELP and returned cross-spectrum modulus
;                              with covariance-propagated error; developed with
;                              assistance from ChatGPT
;       M. Mendez 16 May 2026 renamed manual noise keyword to MANUAL_POI
;                              and clarified that values are in input-file
;                              units before RMS1/RMS2 conversion;
;                              developed with assistance from ChatGPT
;       M. Mendez 15 Jul 2026 added /RAWCOH to force measured coherence
;                              with observed PDS denominator and no noise
;                              correction; developed with assistance
;                              from ChatGPT
;       M. Mendez 16 Jul 2026 suppressed noise-prescription status output
;                              from the preliminary internal pass used for
;                              POISSON Real-part correction; developed
;                              with assistance from ChatGPT
;
;--------------------------------------------------------------------------
;
; MM
if(keyword_set(help)) then begin
   print,''
   print,'GH_CS'
   print,''
   print,'Wrapper around GH_CROSS_RE_IM_NEW.'
   print,'Computes lags, coherence, Re/Im parts, and modulus'
   print,'of the cross spectrum.'
   print,''
   print,'Usage:'
   print,"  gh_cs,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf"
   print,"  gh_cs,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf,poisson=[400,800]"
   print,"  gh_cs,file1,file2,freq,lag,lag_err,coh,coh_err,-100,real,real_err,imag,imag_err,modcs,modcse,mw,qf,manual_poi=[p1,p2]"
   print,''
   print,'Keywords:'
   print,'  /TLAG                 Return time lags'
   print,'  POISSON=[f1,f2]       Frequency range for PDS noise'
   print,'                         subtraction and Real correction'
   print,'  MANUAL_POI=[p1,p2]    Manual PDS1/PDS2 levels in input-file units'
   print,'                         before RMS1/RMS2 conversion'
   print,'                         (coherence only; no Real correction)'
   print,'  RMS1=, RMS2=          Convert Re/Im/Mod to rms^2'
   print,'  /ROTATE               Rotate cross vector by pi/4'
   print,'  /LOWCOH               Use VN97 low-coherence formalism'
   print,'  /RAWCOH               Return measured/raw coherence: observed PDS'
   print,'                         denominator, no Poisson/noise correction'
   print,'  INDEX=, TIME=, RATE=, SEL='
   print,'                         Data selection'
   print,'  /HELP                 Print this message'
   print,''
   print,'Note: with POISSON, GH_CS makes a preliminary internal pass'
   print,'to estimate the Real-part correction, then a final pass that'
   print,'returns the requested products. Only the final pass reports'
   print,'the noise prescription used for the output arrays.'
   print,''
   return
endif
; MM
;
; Do not allow two different user-requested noise prescriptions.
; POISSON uses a frequency range to estimate both PDS noise levels and
; the Real-part correction. MANUAL_POI uses fixed values for PDS1
; and PDS2 only, and does not correct the Real part.
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
; the observed PDS in the denominator and no noise subtraction. Do not
; silently choose between incompatible user requests.
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
; check that if a background is given, both are given.
if((keyword_set(back1)) and ~(keyword_set(back2))) then begin
      massage,'Error: if you give the background you must give it in both bands! (use rms1=back1,rms2=back2)'
      retall
endif

if((keyword_set(back2)) and ~(keyword_set(back1))) then begin
      massage,'Error: if you give the background you must give it in both bands! (use rms1=back1,rms2=back2)'
      retall
endif

;
; First call without "poisson", "rotate" or "back1"/"back2" to get the Real part in Leahy units.
; Suppress informational status from this preliminary pass.

; It now returns the modulus of the CS and error 

gh_cross_re_im_new,                                                                        $
                   filename1,filename2,frequency_reb,lag,lag_err,coherence,coherence_err,  $
                   irf,                                                                    $
                   repart,repart_err,impart,impart_err,modcs,modcs_err,mw,qf,              $
                   index=index,time=time,sel=sel,rate=rate,                                $
                   manual_poi=manual_poi,rawcoh=rawcoh,/status_quiet
poivalue = 0.0

IF(keyword_set(poisson)) THEN BEGIN

   w1 = where(frequency_reb ge poisson[0], n1)
   w2 = where(frequency_reb le poisson[1], n2)

   IF(n1 GT 0 AND n2 GT 0) THEN BEGIN

      p1 = w1[0]
      p2 = w2[n2-1]

      IF(p2 GE p1) THEN BEGIN

         ; --- New (correct) selection ---
         poivalue = mean(repart[p1:p2])

;         print,'NEW uses ',p2-p1+1,' bins: ', $
;               frequency_reb[p1],' - ',frequency_reb[p2], $
;               ' Hz; indices: ',p1,'-',p2, $
;               ' poivalue= ',poivalue
;
;         ; --- What the old code actually did ---
;         wold = where(frequency_reb ge p1 and frequency_reb le p2, nold)
;
;         IF(nold GT 0) THEN BEGIN
;            poivalue_old = mean(repart[wold])
;
;            print,'OLD actually used ',nold,' bins: ', $
;                  frequency_reb[wold[0]],' - ', $
;                  frequency_reb[wold[nold-1]], $
;                  ' Hz; poivalue_old= ',poivalue_old
;         ENDIF ELSE BEGIN
;            print,'OLD selection found no bins'
;         ENDELSE

      ENDIF ELSE BEGIN
         massage,'Poisson frequency range is empty'
         retall
      ENDELSE

   ENDIF ELSE BEGIN
      massage,'Poisson frequency range outside frequency array'
      retall
   ENDELSE

ENDIF

; MM: chatGPT tells me that there was a bug here, that instead of using frequencies I was using indexes.
; Given my littel experience with IDL, I may have been making amkistake here. I'll trust chatGPT
; poivalue = 0.0
; IF(keyword_set(poisson)) THEN BEGIN
;          p1 = where(frequency_reb ge poisson[0])
;          p1 = p1[0]
;          p2 = where(frequency_reb le poisson[1])
;          p2 = p2(n_elements(p2)-1)
;          poivalue=mean(repart(where(frequency_reb ge p1 and frequency_reb le p2)))
; ; debugging
;          print,'Poisson range: ',frequency_reb[p1],'-',frequency_reb[p2], $
;          ' Hz; indices: ',p1,'-',p2,' poivalue= ',poivalue
; 
; ENDIF


; It now returns the modulus of the CS and error 

gh_cross_re_im_new,                                                                           $
                   filename1,filename2,                                                       $
                   frequency_reb,lag,lag_err,coherence,coherence_err,                         $
                   irf,tlag=tl,                                                               $
                   repart,repart_err,impart,impart_err,modcs,modcs_err,mw,qf,                 $
                   index=index,time=time,sel=sel,rate=rate,poisson=poisson,poivalue=poivalue, $
                   manual_poi=manual_poi,                                                     $
                   rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh,rawcoh=rawcoh

end
