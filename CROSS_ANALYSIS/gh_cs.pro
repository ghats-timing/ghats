pro gh_cs,filename1,filename2,             $
             frequency_reb,lag,lag_err,coherence,coherence_err, $
	         irf,tlag=tl, $
                 repart,repart_err,impart,impart_err,mw,qf, $
	         index=index,time=time,sel=sel,rate=rate,poisson=poisson, $
                 rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh
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
;
;       NOTES:
;       The keyword poisson=[freq1,freq2] indicates 2 things:
;       1. Correct the cross spectrum by cross talk/partial correlation (for lags)
;       2. Subtract the Poisson level from the PDS in the 2 individual bands (for coherence function)
;       If the keyword is given, the same frequency range (freq1-freq2) is used for both corrections.

;       For RXTE the intrinsic coherence is computed assuming Zhang 1995 formulae for the Poisson level.
;       If no keyword is given, for any other mission the program assumes that Poisson=2.


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
;			
;
;--------------------------------------------------------------------------
;
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
; First call without "poisson", "rotate" or "back1"/"back2" to get the Real part in Leahy units

gh_cross_re_im_new,                                                                        $
                   filename1,filename2,frequency_reb,lag,lag_err,coherence,coherence_err,  $
                   irf,                                                                    $
                   repart,repart_err,impart,impart_err,                                    $
                   index=index,time=time,sel=sel,rate=rate                                 


poivalue = 0.0
IF(keyword_set(poisson)) THEN BEGIN
         p1 = where(frequency_reb ge poisson[0])
         p1 = p1[0]
         p2 = where(frequency_reb le poisson[1])
         p2 = p2(n_elements(p2)-1)
         poivalue=mean(repart(where(frequency_reb ge p1 and frequency_reb le p2)))
ENDIF


gh_cross_re_im_new,                                                                           $
                   filename1,filename2,                                                       $
                   frequency_reb,lag,lag_err,coherence,coherence_err,                         $
                   irf,tlag=tl,                                                               $
                   repart,repart_err,impart,impart_err,mw, qf,                                   $
                   index=index,time=time,sel=sel,rate=rate,poisson=poisson,poivalue=poivalue, $
                   rms1=back1,rms2=back2,rotate=rotate,lowcoh=lowcoh


end
