pro ghx_nustar,filename1,filename2,             $
             frequency,power1,power1_err, $
	         irf, index=index,time=time,sel=sel,rate=rate,rms=back
;+
; NAME: 
;      GHX_NUSTAR
; PURPOSE: 
;      Reads is two FFTs from FFT files of both NuSTAR telescopes and 
;			produces a cross-spectrum between them. Then, the real part
;			of the cross-spectrum is taken, to isolate the correlated
;			part and leave out the noise. See Bachetti et al. (2015, ApJ, 800, 109).
;			The rms calculation however is incorrect and it will have to be
;			estimated differently (see Bachetti et al. 2015).
;
; CALLING SEQUENCE: 
;        GHX_NUSTAR,FILENAME1,FILENAME2,FREQUENCY,POW,POW_ERR,IRF,
;			      [,INDEX=INDEX][,TIME=TIME][,SEL=SEL]
; INPUTS:
;       FILENAME1 = name of the first  input FFT file
;       FILENAME2 = name of the second input FFT file
;       IRF      = rebin factor (negative for log rebinning)
;
; OUTPUTS:
;       FREQUENCY    = Frequency array
;       POW          = Phase/timelag array
;       POW_ERR      = Array of errors on phase/timelags
;
; KEYWORDS:
;       INDEX    = range of selected indices
;       TIME     = Switch for time lags (default phase lags)
;		SEL      = array with index selection
;		RMS      = background rate for rms conversion
;
; EXAMPLE:
;       Read in two FFTs and compute phaselag spectrum
;
;       MU> GHX_MUSTAR,'tel1.fft','tel2.fft',nu,pow,pow_e,-100
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       MUXANA_OPENFFT:    Opens a FFT file
;       MUXANA_GETHEADER:  Reads in the header
;       READ_FFT_LINE:     Reads in next line from FFT file
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni   9 Dec 2013  from GH_CROSS
;		T. Belloni  22 May 2017  final implementation for testing
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
if(strcmp(strmid(instrument1,0,3),strmid(instrument2,0,3)) eq 0) then begin
   massage,'Instrument keyword not compatible'
   retall
endif
if((T1-T2) ne 0.0) then begin
   massage,'Time resolution not compatible'
   retall
endif
if(nft1 ne nft2) then begin
   massage,'Frequency resolution not compatible'
   retall
endif
if((rmjd01-rmjd02) ne 0.0) then begin
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
if(keyword_set(rate)) then begin
    firstrate = rate[0]
    lastrate  = rate[1]
endif else begin
    firstrate = 0.0
    lastrate  = 100000.0
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
   rate = cnts*T1
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
      (rate ge firstrate) and (rate le lastrate) and $
      (gotcha[0] ge 0)) then begin  
;  -----------------------------
      cc         = conj(rdata1) * rdata2 * 2.0 /sqrt(cnts1*cnts2) ; cross spectrum
;      pwr1       = abs(cc)^2
      pwr1       = float(cc)     
      cc         = conj(rdata2) * rdata1 * 2.0 /sqrt(cnts1*cnts2) ; cross spectrum
;      pwr2       = abs(cc)^2
      pwr2    = float(cc)
	  
      power1     = power1 + pwr1
      power1_err = power1_err + pwr1*pwr1
      power2     = power2 + pwr2
      power2_err = power2_err + pwr2*pwr2


      n_selected = n_selected + 1
      flux      = flux + cnts * T
   endif

   if(itrafos ge (lastindex-1)) then goto,finito
;
endfor

finito:

;
;  Compute average and standard deviation
;
   if(n_selected gt 0) then begin
      flux      = flux      / n_selected

      power1    = power1    / n_selected
      power2    = power2    / n_selected
      power1_err= power1    / sqrt(n_selected)
      power2_err= power2    / sqrt(n_selected)
      
      power1    = power1 - poilevel1 
      power2    = power2 - poilevel2
      
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

   if(keyword_set(back)) then begin
;
;  Rms^2 conversion
;
      if(back ge flux) then begin
         massage,'Error: background flux higher than source+bkg flux!'
	 retall
      endif
      normalization = flux /(flux - back)^2
      power1         = power1     * normalization
      power1_err     = power1_err * normalization
      power2         = power2     * normalization
      power2_err     = power2_err * normalization
   endif

end
