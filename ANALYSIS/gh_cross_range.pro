pro gh_cross_range,filename1,filename2,nu0,width,             $
             frequency_reb,lag,lag_err,irfplot,                       $
	         tlag=tl, $
	         index=index,time=time,sel=sel,plot=pl,poivalue=poivalue
;+
; NAME: 
;      GH_CROSS_RANGE
; PURPOSE: 
;      Reads is two FFTs from FFT files, selecting on a range of FFTs
;			in time, range, or selection of indices, and produces the
;			averaged time/phase lag around a given centroid frequency.
; CALLING SEQUENCE: 
;       MUCROSS,FILENAME1,FILENAME2,CENTROID,WIDTH,FREQUENCY,LAG,LAG_ERR,
;                 IRF,[,/TIMELAG][,INDEX=INDEX][,TIME=TIM][,SEL=SEL],[/PLOT][,POIVALUE=POIVALUE]
; INPUTS:
;       FILENAME1 = name of the first  input FFT file
;       FILENAME2 = name of the second input FFT file
;       STARTING = index of first FFT to be selected
;       ENDING   = index of last  FFT to be selected
;       IRF      = rebin factor (negative for log rebinning) for plot only
;
; OUTPUTS:
;       FREQUENCY    = Frequency array
;       LAG          = Phase/timelag array
;       LAG_ERR      = Array of errors on phase/timelags
;
; KEYWORDS:
;       TIMELAG       = Switch for time lags (default phase lags)
;       TIME          = Array for time selection
;		    INDEX         = Array for index range selection
;		    SEL           = Array for index selection
;       PLOT          = Flag for plot of lag spectrum with indicated range
;       POIVALUE      = Real part of noise average cross-spectrum (for subtraction)
;
; EXAMPLE:
;       Read in two FFTs and compute phaselag spectrum
;
;       MU> GH_CROSS_RANGE'data/4u1630a.fft','data/4u1630b.fft',5.5,0.5,nu,pha,pha_e,-60,/plot
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
;		T. Belloni  27 Jan 2012  adapted from gh_cross
;		T. Belloni   9 Feb 2012  check for zero or one bins in selection
;		T. Belloni   9 Feb 2012  time range bug fixed
;		T. Belloni   9 Feb 2012  /TIMELAG keyword changed to /TLAG
;		T. Belloni  20 Feb 2012  no floor rounding for times array
;		T. Belloni   7 Jun 2012  coherence calculation had Poisson estimate only valid for RXTE/PCA
;  		T. Belloni  12 Jun 2012  increased tolrance for time check 
;   	T. Belloni  26 Feb 2013  added output of average cross spectrum value
;  		T. Belloni  26 Feb 2013  added keyword for input real part to subtract
;		T. Belloni  19 Jun 2013  bug for time lag extraction corrected
;		T. Belloni   5 Dec 2013  free_lun replaces close
;		T. Belloni  20 Feb 2017  fixed /tlag behavior
;		T. Belloni  27 Aug 2020  output cross-spectral amplitude with error
;-
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
;if(strcmp(instrument1,instrument2) eq 0) then begin
;   massage,'Instrument keyword not compatible'
;   retall
;endif
if((T1-T2) gt 1.e-9) then begin
   massage,'Time resolution not compatible'
   retall
endif
if(nft1 ne nft2) then begin
   massage,'Frequency resolution not compatible'
   retall
endif
if((rmjd01-rmjd02) ge 1.0e-8) then begin
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
; Loop over the FFTs
;
n_selected       = 0l
itrafos          = 0l
flux             = 0.0d0
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
if(keyword_set(sel)) then begin
    goodarray = sel-1
endif else begin
    goodarray = indgen(ntrafos)
endelse

goodarray = goodarray(where((goodarray ge firstpds) and (goodarray le lastpds)))
lastindex = max(goodarray)

a0 = 0.0
rmjd1 = 0d0
rmjd2 = 0d0

for itrafos=0l,ntrafos-1l do begin
;
   read_fft_line,unit1,mufla1,rmjd1,cnts1,poisson1,current_vle_rate1,fndet1,rdata1
   read_fft_line,unit2,mufla2,rmjd2,cnts2,poisson2,current_vle_rate2,fndet2,rdata2
   cnts = cnts1+cnts2
   timediff = abs(rmjd1-rmjd2)

   if(timediff gt 1.e-4) then begin
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
      pwr      = abs(rdata1)^2 * 2.0/cnts1
      power1   = power1 + pwr
      pwr      = abs(rdata2)^2 * 2.0/cnts2
      power2   = power2 + pwr

      n_selected = n_selected + 1
      flux      = flux + cnts * T
      
; Here Poissonian estimates for the two bands
      if(strcmp(instrument1,'PCA             ') eq 0) then begin
         poisson_estimate,poilevel1,cnts1,1.0/T,current_vle_rate1,nft,fndet1,i_vle,tdead
         poisson_estimate,poilevel2,cnts2,1.0/T,current_vle_rate2,nft,fndet2,i_vle,tdead
        endif 	else begin
		 fndet1 = 2.0
		 fndet2 = 2.0
	  endelse
   endif
   if(itrafos ge (lastindex-1)) then goto,finito
;
endfor

finito:

print,n_selected
;
;  Compute average and standard deviation
;
   if(n_selected gt 0) then begin
      flux      = flux      / n_selected
      poilevel1 = poilevel1 / n_selected
      poilevel2 = poilevel2 / n_selected
      power1    = power1    / n_selected
      power2    = power2    / n_selected
      power1_err= power1    / sqrt(n_selected)
      power2_err= power2    / sqrt(n_selected)
      
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
; Here the cross spectrum is available. We need to cut only the part that we want
;
   nu1 = nu0 - width/2.0
   nu2 = nu0 + width/2.0

   buoni = where((frequency ge nu1) and (frequency le nu2))
   nbuoni = n_elements(buoni)

   if(nbuoni eq 0) then begin
	    print,'No bins within frequency selection. Increase the frequency resolution'
	    print,'by using longer time intervals or less rebinning'
    	return
   endif
   if(nbuoni eq 1) then begin
	    print,'Only one bin within frequency selection. Increase the frequency resolution'
	    print,'by using longer time intervals or less rebinning'
    	return
   endif


   frequency = frequency(buoni)
   flux      = flux(buoni)
   poilevel1 = poilevel1(buoni)
   poilevel2 = poilevel1(buoni)
   power1    = power1(buoni)
   power2    = power2(buoni)
   power1_err    = power1_err(buoni)
   power2_err    = power2_err(buoni)
   ccsumjk       = ccsumjk(buoni)

   irf = nbuoni
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
;coh     = (csq - nsq)/(y1*y2)

; Calculate the error on coherence (VN97 eqn 8)

;coherr1 = (2.0*(nsq^2)*mw)/(csq-nsq)^2
;coherr2 = (nslev1/y1)^2 + (nslev2/y2)^2
;coherr3 = 2.0*((1-coh)/(coh^1.5))^2
;coherr  = coh/sqrt(mw) * sqrt(coherr1+coherr2+coherr3)

; Phase and phase error (from Nowak et al. 1999). Note that the phase error
; use coherence uncorrected for noise levels, so we add the noise back on to 
; the powers and we use csq in the numerator since we do not want to subtract nsq

; If an average real part of the noise cross-spectrum (evaluated in a previous run of the same
; function) is set, subtract it before computing the lag

 if(keyword_set(poivalue)) then begin
      ccsumjk_reb = ccsumjk_reb - complex(poivalue,0)
   endif

lag       = atan(imaginary(ccsumjk_reb),float(ccsumjk_reb))
phaserr1  = csq/((y1+nslev1)*(y2+nslev2))
lag_err   = sqrt((1.0-phaserr1)/(2.0*phaserr1*mw))

;coherence     = coh
;coherence_err = coherr
;
;
;  Switch for time lags
;
   uscita     = 'Phase'
   if(keyword_set(tl)) then begin
      lag     = lag/2.0/!PI/frequency_reb
      lag_err = lag_err/2.0/!PI/frequency_reb
      uscita  = 'Time'
   endif

; Terminal output
print,uscita,' lag between ',nu1, ' and ',nu2,' Hz : ',lag, ' +/- ',lag_err
print,uscita,' Complex cross-spectrum value : ',ccsumjk_reb    ; NEW 2013 Feb 13
print,uscita,' Cross-spectrum amplitude : ',abs(ccsumjk_reb), ' +/- ',sqrt(y1*y2/mw) ; NEW 2020 Aug 27

   if(keyword_set(pl)) then begin
	  if(keyword_set(tl)) then begin
		if(keyword_set(time)) then begin
			gh_cross,filename1,filename2,f,l,l_e,c,ce,irfplot,/tlag,time=time
 	       endif else begin
        	gh_cross,filename1,filename2,f,l,l_e,c,ce,irfplot,/tlag
         endelse
        gh_plot_lag,f,l,l_e,/time,/lin
        oplot,[nu1,nu1],[-!PI,!PI],linestyle=1
        oplot,[nu2,nu2],[-!PI,!PI],linestyle=1
        oplot,[nu0,nu0],[-!PI,!PI],linestyle=2
      endif else begin
	    if(keyword_set(time)) then begin
		    gh_cross,filename1,filename2,f,l,l_e,c,ce,irfplot,time=time
		 endif else begin
    	    gh_cross,filename1,filename2,f,l,l_e,c,ce,irfplot
        endelse
    	gh_plot_lag,f,l,l_e,/lin
        oplot,[nu1,nu1],[-!PI,!PI],linestyle=1
        oplot,[nu2,nu2],[-!PI,!PI],linestyle=1
        oplot,[nu0,nu0],[-!PI,!PI],linestyle=2
	  endelse
   endif

end
