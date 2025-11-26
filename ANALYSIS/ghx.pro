pro ghx,filename,              $
             frequency,power,power_err, $
             index=index,time=time,rate=rate,sel=sel, $
	         poisson=poi,rms=back;,czti=czti	     
;+
; NAME: 
;      ghx
; PURPOSE: 
;      Reads is a PDS from a PDS file, selecting on a range of PDS if desired
; EXPLANATION:
;      This procedure extracts an average PDS from a PDS file. 
;
; CALLING SEQUENCE: 
;       GHX,FILENAME,FREQUENCY,POWER,POWER_ERR,
;                [,INDEX=INDEX][,TIME=TIME][,RATE=RATE][,SEL=SEL]
;                [,POISSON][,RMS=RMS]
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       FREQUENCY= Frequency array
;       POWER    = Power array
;       POWER_ERR= Array of errors on power
;
; KEYWORDS:
;       INDEX    = range of selected indices
;       TIME     = range of selected times
;       RATE     = range of selected rates
;       SEL      = array with index selection
;       RMS      = background value for rms2 conversion. If specified, 
;                  rms2 conversion will be performed
;       POISSON  = i_vle for Poissonian calculation. If set, the 
;                  poissonian contribution will be subtracted
;
; EXAMPLE:
;       Read in all PDS in a PDS file
;
;       MU> GHX,'data/4u1630-47.pds',nu,pow,pow_e
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       GHATS_OPENPDS:    Opens a PDS file
;       GHATS_GETHEADER:  Reads in the header
;       READ_PDS_LINE:     Reads in next line from PDS file
;       POISSON_ESTIMATE   Computes Poissonian level
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  11 Nov 2001  modular version + rearranged order of pars +
;                                keyword+fast exit
;       T. Belloni  12 May 2002  adapted to MUFFT format
;       T. Belloni  17 May 2002  Poissonian subtraction added
;       T. Belloni  10 Jun 2002  rms conversion added
;       T. Belloni  08 Aug 2002  added skip to end when last good spectrum read
;       T. Belloni  17 Dec 2002  added variable i_vle
;       T. Belloni  11 Dec 2003  Mu3 version (colors)
;       T. Belloni  29 Mar 2009  I_VLE now read in + conditional poisson
;       T. Belloni  07 Apr 2010  mux from muxana*
;		T. Belloni  05 May 2010  ghx from mux
;		T. Belloni  03 Dec 2010  default rate limit increased
;		T. Belloni  20 Feb 2012  no floor rounding for times array
;		T. Belloni  19 Dec 2013  free_lun added
;       T. Belloni  12 Jun 2015  increased count rate range to allow for Sco X-1
;		T. Belloni  05 May 2017  SEL indices changed: now from 0 and not from 1
;		T. Belloni  09 Nov 2017  CZTI background added - currently 
;		T. Belloni  29 May 2019  goodarray now long, not limited to 32768 spectra
;-
;--------------------------------------------------------------------------
;
; Open PDF file
;
ghats_openpds,filename,unit,/dialog
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

;muxana_getheader,unit,nft,T,observatory,target,rmjd0,instrument,ntrafos,e,dummy *MU*

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
pwr       = fltarr(nft)*0.0
power     = fltarr(nft)*0.0
power_err = fltarr(nft)*0.0

poilevel  = fltarr(nft)*0.0
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
   read_pds_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,a0,pwr
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

      if(keyword_set(poi)) then begin
       poisson_estimate,poilevel,cnts,1.0/T,current_vle_rate,nft,fndet,i_vle,tdead
      endif
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
      poilevel  = poilevel / n_selected           ; poilevel must be average
      flux      = flux     / n_selected           ; and so must the flux
     endif else begin
      print,'Warning! No power spectra retrieved for input selection!'
   endelse
close,unit
free_lun,unit
;
;  Now fill in frequency array
;
   frequency = (findgen(nft)+1.0) * T
;
   ps = ' spectra '
   if(n_selected eq 1) then ps = ' spectrum '
   print,'  ',strtrim(string(n_selected),1),' power'+ps+'selected'
;
;  If requested, subtract Poissonian component
;
   if(keyword_set(poi)) then begin
      power = power - poilevel
   endif
   if(keyword_set(back)) then begin
;
;  Rms^2 conversion
;
      if(back ge flux) then begin
         massage,'Error: background flux higher than source+bkg flux!'
	 retall
      endif
      normalization = flux /(flux - back)^2
      power         = power     * normalization
      power_err     = power_err * normalization
   endif
   
;   if(keyword_set(czti)) then begin
;
;  Rms^2 conversion for CZTI, background comes from the file itself, a correction factor is given as keyword czti
; for the moment the correction is not used
;
;      czti = 1  ; TEMPORARY!
;      normalization = flux /((flux - fondo)*czti)^2
;      power         = power     * normalization
;      power_err     = power_err * normalization
;   endif
   
;print,'Inizio: ',inizio_run  ; for test purposes
;print,'Fine  : ',systime(0)  ; for test purposes

end
