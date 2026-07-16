pro ghx_poisrange,filename,              $
             frequency,power,power_err, $
             index=index,time=time,rate=rate,sel=sel, $
             poisson=poi,rms=back,poisrange=poisrange,help=help
;+
; NAME:
;      GHX_POISRANGE
; PURPOSE:
;      Reads a PDS from a PDS file, selecting on a range of PDS if desired.
;      Optionally subtracts an empirical constant Poisson/noise level
;      estimated from a user-given frequency range.
;
; EXPLANATION:
;      This procedure extracts an average PDS from a PDS file.
;
;      It is based on GHX, but adds the POISRANGE keyword. This allows
;      the user to reproduce the usual interactive sequence:
;
;          ghx,'file.pds',nu,pow,powe,rms=0.001
;          z  = mean(pow(where(nu gt f1 and nu lt f2)))
;          po = pow - z
;
;      directly inside the routine, with a safer check that the selected
;      frequency range is not empty.
;
;      The POISSON keyword keeps its original meaning: it applies the
;      RXTE/PCA Zhang Poisson formulae. Since that correction is only valid
;      for RXTE/PCA data, this routine stops if POISSON is requested for a
;      non-RXTE file.
;
; CALLING SEQUENCE:
;      GHX_POISRANGE,FILENAME,FREQUENCY,POWER,POWER_ERR, $
;             [,INDEX=INDEX][,TIME=TIME][,RATE=RATE][,SEL=SEL] $
;             [,POISSON][,RMS=RMS][,POISRANGE=POISRANGE][,/HELP]
;
; INPUTS:
;      FILENAME = name of the input PDS file
;
; OUTPUTS:
;      FREQUENCY = frequency array
;      POWER     = power array
;      POWER_ERR = array of errors on power
;
; KEYWORDS:
;      INDEX     = range of selected PDS indices
;      TIME      = range of selected times
;      RATE      = range of selected rates
;      SEL       = array with index selection
;
;      RMS       = background value for rms^2 conversion. If specified,
;                  rms^2 conversion will be performed.
;
;      POISSON   = if set, the RXTE/PCA Zhang Poisson level is computed
;                  and subtracted. This is valid only for RXTE data.
;                  If POISSON is set for a non-RXTE file, the routine stops.
;
;      POISRANGE = two-element frequency range [f1,f2]. If set, the mean
;                  power in this range is subtracted from the whole PDS
;                  after optional rms conversion. This is useful to subtract
;                  a residual constant Poisson/noise level empirically.
;
;      HELP      = if set, print a short usage message and return.
;
; EXAMPLE:
;      Read in all PDS in a PDS file:
;
;      IDL> ghx_poisrange,'data/4u1630-47.pds',nu,pow,pow_e
;
;      Read in all PDS, convert to rms units, and subtract an empirical
;      Poisson level estimated between 400 and 800 Hz:
;
;      IDL> ghx_poisrange,'event_30-1200.pds',nu,pow,pow_e, $
;             rms=0.001,poisrange=[400,800]
;
;      Then rebin and write in XSPEC format:
;
;      IDL> gh_reb,nu,pow,pow_e,-100,x,y,ye
;      IDL> gh_xspec,x,y,ye,'event_30-1200'
;
;      RXTE/PCA case with Zhang Poisson subtraction plus an empirical
;      residual correction:
;
;      IDL> ghx_poisrange,'event_30-1200.pds',nu,pow,pow_e, $
;             poisson=1,rms=0.001,poisrange=[400,800]
;
; COMMON BLOCKS:
;      None
;
; ROUTINES USED:
;      GHATS_OPENPDS
;      GHATS_GETHEADER
;      READ_PDS_LINE
;      POISSON_ESTIMATE
;
; NOTES:
;      POISRANGE is applied after optional rms conversion, so that it
;      reproduces the usual interactive workflow in which the high-frequency
;      constant level is estimated from the rms-normalised PDS.
;
; MODIFICATION HISTORY:
;      T. Belloni  20 Aug 2001  implementation of original GHX
;      T. Belloni  05 May 2010  GHX from MUX
;      M. Mendez   14 May 2026  GHX_POISRANGE, based on GHX; added
;                               POISRANGE keyword, safer empty-range check,
;                               RXTE-only check for POISSON, and /HELP,
;                               with help from ChatGPT
;-
;--------------------------------------------------------------------------

;
; Help screen
;
if(keyword_set(help)) then begin
   print,' '
   print,'GHX_POISRANGE'
   print,'Read and average a GHATS PDS file.'
   print,' '
   print,'Usage:'
   print,"  ghx_poisrange,'file.pds',nu,pow,pow_e"
   print,"  ghx_poisrange,'file.pds',nu,pow,pow_e,rms=0.001"
   print,"  ghx_poisrange,'file.pds',nu,pow,pow_e,rms=0.001,poisrange=[400,800]"
   print,' '
   print,'Keywords:'
   print,'  INDEX     selected PDS index range'
   print,'  TIME      selected time range'
   print,'  RATE      selected count-rate range'
   print,'  SEL       selected PDS indices'
   print,'  RMS       background level for rms^2 normalisation'
   print,'  POISSON   RXTE/PCA Zhang Poisson subtraction only'
   print,'  POISRANGE empirical constant level from [f1,f2] Hz'
   print,' '
   print,'Example:'
   print,"  ghx_poisrange,'event_30-1200.pds',nu,pow,powe,rms=0.001,poisrange=[400,800]"
   print,'  gh_reb,nu,pow,powe,-100,x,y,ye'
   print,"  gh_xspec,x,y,ye,'event_30-1200'"
   print,' '
   return
endif

;--------------------------------------------------------------------------
; Open PDS file
;--------------------------------------------------------------------------
ghats_openpds,filename,unit,/dialog

;--------------------------------------------------------------------------
; Read in PDS file header
;--------------------------------------------------------------------------
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
i_vle  = fix([dummy(17:18)],0)+1

;--------------------------------------------------------------------------
; The POISSON keyword uses RXTE/PCA-specific Zhang formulae.
; Do not allow silent use for other missions.
;--------------------------------------------------------------------------
use_poisson = keyword_set(poi)

if(use_poisson) then begin
   if(strtrim(observatory,2) ne 'XTE') then begin
      massage,'POISSON keyword uses RXTE/PCA Zhang formulae and is only valid for RXTE data. Use POISRANGE for empirical subtraction.'
      retall
   endif
endif

if(use_poisson and keyword_set(poisrange)) then begin
   print,'Applying RXTE Poisson subtraction + residual correction'
endif

;--------------------------------------------------------------------------
; nft is number of TIME points. Divide by 2 to get number of positive
; Fourier frequencies stored in the PDS.
;--------------------------------------------------------------------------
nft = nft/2

pwr       = fltarr(nft)*0.0
power     = fltarr(nft)*0.0
power_err = fltarr(nft)*0.0
poilevel  = fltarr(nft)*0.0

;--------------------------------------------------------------------------
; Loop over the trafos
;--------------------------------------------------------------------------
n_selected = 0l
itrafos    = 0l
tdead      = 1.0e-5
flux       = 0.0d0
fondo      = 0.0d0

;--------------------------------------------------------------------------
; Data selection
;--------------------------------------------------------------------------
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
   goodarray = sel
endif else begin
   goodarray = lindgen(ntrafos)
endelse

goodarray = goodarray(where((goodarray ge firstpds) and (goodarray le lastpds)))
lastindex = max(goodarray)

for itrafos=0l,ntrafos-1l do begin

   read_pds_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,a0,pwr

   ; Accumulate average power only for selected spectra.
   t_1    = ((rmjd-rmjd0)*86400.0)
   t_2    = t_1 + 1.0/T
   cratem = cnts*T
   gotcha = where(goodarray eq itrafos)

   if((t_1 ge firsttime) and (t_2 le lasttime) and $
      (cratem ge firstrate) and (cratem le lastrate) and $
      (gotcha[0] ge 0)) then begin

      n_selected = n_selected + 1
      power      = power     + pwr
      power_err  = power_err + pwr*pwr
      flux       = flux + cnts * T
      fondo      = fondo + current_vle_rate * T

      ; RXTE/PCA Poisson estimate from Zhang formulae.
      ; This block is reached only if the header check above confirmed XTE.
      if(use_poisson) then begin
         poisson_estimate,poilevel,cnts,1.0/T,current_vle_rate,nft,fndet,i_vle,tdead
      endif

   endif

   if(itrafos gt (lastindex-1)) then goto,finito

endfor

finito:

;--------------------------------------------------------------------------
; Compute average and errors
;--------------------------------------------------------------------------
if(n_selected gt 0) then begin
   power     = power / n_selected
   power_err = power / sqrt(n_selected)
   poilevel  = poilevel / n_selected
   flux      = flux / n_selected
endif else begin
   print,'Warning! No power spectra retrieved for input selection!'
endelse

close,unit
free_lun,unit

;--------------------------------------------------------------------------
; Frequency array
;--------------------------------------------------------------------------
frequency = (findgen(nft)+1.0) * T

ps = ' spectra '
if(n_selected eq 1) then ps = ' spectrum '
print,'  ',strtrim(string(n_selected),1),' power'+ps+'selected'

;--------------------------------------------------------------------------
; If requested, subtract RXTE/PCA Poissonian component.
;--------------------------------------------------------------------------
if(use_poisson) then begin
   power = power - poilevel
endif

;--------------------------------------------------------------------------
; Rms^2 conversion.
;--------------------------------------------------------------------------
if(keyword_set(back)) then begin

   if(back ge flux) then begin
      massage,'Error: background flux higher than source+bkg flux!'
      retall
   endif

   normalization = flux /(flux - back)^2
   power         = power     * normalization
   power_err     = power_err * normalization

endif

;--------------------------------------------------------------------------
; If requested, estimate a residual constant Poisson/noise level from a
; user-given frequency range and subtract it from the whole PDS.
;
; This is done after optional rms conversion, to reproduce the usual
; interactive workflow:
;
;      ghx,...,rms=0.001
;      z  = mean(pow(where(nu gt f1 and nu lt f2)))
;      po = pow - z
;
; The WHERE call uses the frequency values directly. The npoi check avoids
; the dangerous IDL behaviour in which WHERE returns -1 when the selection
; is empty.
;--------------------------------------------------------------------------
if(keyword_set(poisrange)) then begin

   if(n_elements(poisrange) ne 2) then begin
      massage,'POISRANGE must have two elements: [f1,f2]'
      retall
   endif

   wpoi = where(frequency gt poisrange[0] and frequency lt poisrange[1], npoi)

   if(npoi gt 0) then begin
      poivalue = mean(power[wpoi])
      power    = power - poivalue

      print,'Poisson level estimated as the average between ', $
            frequency[wpoi[0]],' and ',frequency[wpoi[npoi-1]], $
            ' Hz: ',poivalue

   endif else begin
      massage,'POISRANGE outside frequency array'
      retall
   endelse

endif

end
