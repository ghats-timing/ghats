pro ghx_fft,filename,              $
             frequency,power,power_err, $
             index=index,time=time,rate=rate,sel=sel	     
;+
; NAME: 
;      ghx_fft
; PURPOSE: 
;      Reads is an FFT from an FFT file, selecting on a range of FFTs if desired
; EXPLANATION:
;      This procedure extracts an average FFT from an FFT file. 
;
; CALLING SEQUENCE: 
;       GHX_FFT,FILENAME,FREQUENCY,FFTX,FFTX_ERR,
;                [,INDEX=INDEX][,TIME=TIME][,RATE=RATE][,SEL=SEL]
; INPUTS:
;       FILENAME = name of the input FFT file
;
; OUTPUTS:
;       FREQUENCY = Frequency array
;       POWER     = FFT array
;       POWER_ERR = Array of errors on FFT
;
; KEYWORDS:
;       INDEX    = range of selected indices
;       TIME     = range of selected times
;       RATE     = range of selected rates
;       SEL      = array with index selection
;			
; EXAMPLE:
;       Read in all FFT in an FFT file
;
;       MU> GHX_FFT,'data/4u1630-47.fft',nu,ff,ff_e
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       GHATS_OPENFFT:    Opens an FFT file
;       GHATS_GETHEADER:  Reads in the header
;       READ_FFT_LINE:     Reads in next line from FFT file
;       POISSON_ESTIMATE   Computes Poissonian level
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  04 Jun 2021  from GHX
;-
;--------------------------------------------------------------------------
;
; Open FFT file
;
ghats_openfft,filename,unit,/dialog
;inizio_run = systime(0)
;
; Read in PDS file header
;
ntrafos             = 0l
dummy               = bytarr(100)
gh_version_string   = '                '
observatory         = '                '
instrument          = '                '
target              = '                '
rmjd0               = 0.0D0
goodarray           = 0l

ghats_getheader,unit,gh_version_string,observatory,instrument,target,rmjd0, $
                     nft,T,ntrafos,e,proliferation,baryflag,n_spectral_bins, $
                     background_flag,dummy
muflag = dummy(0)
i_vle = fix([dummy(17:18)],0)+1
;
;   nft is number of TIME points. Must be divided by 2 to get frequencies
;
nft       = nft/2

;
pwr       = complexarr(nft)*0.0
power     = complexarr(nft)*0.0
power_err = complexarr(nft)*0.0
;
; Loop over the trafos
;
n_selected       = 0l
itrafos          = 0l
tdead            = 1.0e-5
flux             = 0.0d0
fondo            = 0.0d0

; --------------- Data selection ------------------

if(keyword_set(index)) then begin
    firstpds = index[0]
    lastpds  = index[1]
endif else begin
    firstpds = 0L
    lastpds  = 10000000L
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
    goodarray = sel   ; used to be sel-1
endif else begin
    goodarray = lindgen(ntrafos)
endelse

goodarray = goodarray(where((goodarray ge firstpds) and (goodarray le lastpds)))
lastindex = max(goodarray)

for itrafos=0l,ntrafos-1l do begin
;
   read_fft_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,pwr
;  Accumulate average power
;
   t_1 = ((rmjd-rmjd0)*86400.0)
   t_2 = t_1 + 1.0/t
   rate = cnts*T
   gotcha = where(goodarray eq itrafos)

   if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (rate ge firstrate) and (rate le lastrate) $
       and (gotcha[0] ge 0)) then begin
      n_selected = n_selected + 1
      power     = power     + pwr
      power_err = power_err + pwr*pwr
      flux      = flux + cnts * T
	  fondo     = fondo + current_vle_rate * T
   endif
   if(itrafos gt (lastindex-1)) then goto,finito
;
endfor

finito:
;
;  Compute average and standard deviation
;
   if(n_selected gt 0) then begin
      power     = power     / n_selected
      power_err = power / sqrt(n_selected)
      flux      = flux     / n_selected           ; and so must the flux
     endif else begin
      print,'Warning! No FFT retrieved for input selection!'
   endelse
close,unit
free_lun,unit
;
;  Now fill in frequency array
;
   frequency = (findgen(nft)+1.0) * T
;
   print,'  ',strtrim(string(n_selected),1),' FFT '+'selected'
end
