pro ghx,filename,              $
             frequency,power,power_err, $
             index=index,time=time,rate=rate,sel=sel, $
             poisson=poi,manual_poi=manual_poi, $
             rms=back,help=help;,czti=czti
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
;                [,/POISSON] [,POISSON=[f1,f2]]
;                [,MANUAL_POI=value] [,RMS=RMS]
; INPUTS:
;       FILENAME = name of the input PDS file
;                  If FILENAME starts with '@', it is interpreted as a
;                  metafile containing one PDS filename per line.
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
;       POISSON  = /POISSON uses Zhang-style Poisson subtraction
;                  (legacy behaviour)
;                  POISSON=[f1,f2] estimates a constant noise level
;                  from a frequency range and subtracts it
;       MANUAL_POI = Manual constant noise level in input-file units
;                    before RMS conversion
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
;       M. Mendez 16 May 2026 added /HELP, MANUAL_POI keyword,
;                              Poisson range handling, and checks against
;                              incompatible noise-subtraction options;
;                              fixed Poisson-level accumulation bug;
;                              developed with assistance from ChatGPT
;       M. Mendez  08 Jun 2026 added @metafile support for compatible PDS files
;                              developed with assistance from ChatGPT
;-
;--------------------------------------------------------------------------
;MM
if(keyword_set(help)) then begin
   print,''
   print,'GHX'
   print,''
   print,'Extract average PDS from a GHATS PDS file.'
   print,''
   print,'Usage:'
   print,"  ghx,'file.pds',freq,pow,powe"
   print,"  ghx,'@pds_files.lis',freq,pow,powe"
   print,"  ghx,'file.pds',freq,pow,powe,/poisson"
   print,"  ghx,'file.pds',freq,pow,powe,poisson=[400,800]"
   print,"  ghx,'file.pds',freq,pow,powe,manual_poi=p"
   print,''
   print,'Keywords:'
   print,'  /POISSON              Zhang-style Poisson subtraction'
   print,'  POISSON=[f1,f2]       Estimate constant level from range'
   print,'  MANUAL_POI=p          Manual constant level in input-file units'
   print,'                         before RMS conversion'
   print,'  RMS=                  Convert to rms^2 units'
   print,'  INDEX=, TIME=, RATE=, SEL='
   print,'                         Data selection'
   print,'  /HELP                 Print this message'
   print,''
   return
endif
;MM;
; Build list of input PDS files. A scalar FILENAME keeps the historical
; behaviour; '@file' is a metafile with one PDS filename per line.
;
if(strpos(filename,'@') eq 0) then begin
   ghx_read_pds_metafile,filename,pdsfiles
   metafile_flag = 1
endif else begin
   pdsfiles = [filename]
   metafile_flag = 0
endelse

n_pds_files = n_elements(pdsfiles)
ntrafos_files = lonarr(n_pds_files)
rmjd0_files   = dblarr(n_pds_files)

ntrafos             = 0l
dummy               = bytarr(100)
gh_version_string   = '                '
observatory         = '                '
instrument          = '                '
target              = '                '
rmjd0               = 0.0D0
goodarray           = 0l

for ipds=0,n_pds_files-1 do begin

   if(metafile_flag eq 1) then begin
      ghats_openpds,pdsfiles[ipds],unit
   endif else begin
      ghats_openpds,pdsfiles[ipds],unit,/dialog
   endelse

   this_dummy               = bytarr(100)
   this_gh_version_string   = '                '
   this_observatory         = '                '
   this_instrument          = '                '
   this_target              = '                '
   this_rmjd0               = 0.0D0
   this_ntrafos             = 0l
   this_e                   = intarr(2)

   ghats_getheader,unit,this_gh_version_string,this_observatory, $
                        this_instrument,this_target,this_rmjd0, $
                        this_nft,this_T,this_ntrafos,this_e, $
                        this_proliferation,this_baryflag, $
                        this_n_spectral_bins,this_background_flag, $
                        this_dummy

   close,unit
   free_lun,unit

   if(ipds eq 0) then begin
      gh_version_string = this_gh_version_string
      observatory       = this_observatory
      instrument        = this_instrument
      target            = this_target
      rmjd0             = this_rmjd0
      nft               = this_nft
      T                 = this_T
      e                 = this_e
      proliferation     = this_proliferation
      baryflag          = this_baryflag
      n_spectral_bins   = this_n_spectral_bins
      background_flag   = this_background_flag
      dummy             = this_dummy
   endif else begin
      if(this_nft ne nft) then begin
         massage,'Incompatible PDS files in metafile: different NFT'
         retall
      endif
      if(this_T ne T) then begin
         massage,'Incompatible PDS files in metafile: different frequency resolution'
         retall
      endif
      if(total(abs(this_e-e)) ne 0) then begin
         massage,'Incompatible PDS files in metafile: different channel range'
         retall
      endif
      if(this_proliferation ne proliferation) then begin
         massage,'Incompatible PDS files in metafile: different proliferation'
         retall
      endif
      if(this_baryflag ne baryflag) then begin
         massage,'Incompatible PDS files in metafile: different barycenter flag'
         retall
      endif
      if(this_n_spectral_bins ne n_spectral_bins) then begin
         massage,'Incompatible PDS files in metafile: different spectral-bin layout'
         retall
      endif
      if(this_background_flag ne background_flag) then begin
         massage,'Incompatible PDS files in metafile: different background flag'
         retall
      endif
      if(this_dummy[0] ne dummy[0]) then begin
         massage,'Incompatible PDS files in metafile: different GHATS mode flags'
         retall
      endif
      if(this_dummy[1] ne dummy[1]) then begin
         massage,'Incompatible PDS files in metafile: different window flag'
         retall
      endif
      if(total(abs(this_dummy[17:18]-dummy[17:18])) ne 0) then begin
         massage,'Incompatible PDS files in metafile: different VLE flag'
         retall
      endif
   endelse

   rmjd0_files[ipds]   = this_rmjd0
   ntrafos_files[ipds] = this_ntrafos
   ntrafos             = ntrafos + this_ntrafos

endfor

if(metafile_flag eq 1) then begin
   sort_order = sort(rmjd0_files)
   if(total(abs(sort_order-lindgen(n_pds_files))) ne 0) then begin
      print,'WARNING: PDS files in metafile are not in time order.'
      print,'         GHX will process them sorted by header start time.'
      pdsfiles       = pdsfiles[sort_order]
      rmjd0_files    = rmjd0_files[sort_order]
      ntrafos_files  = ntrafos_files[sort_order]
   endif

   trafo_duration = 1.0d0 / double(T)
   for ipds=1,n_pds_files-1 do begin
      prev_end = rmjd0_files[ipds-1] + $
                 double(ntrafos_files[ipds-1])*trafo_duration/86400.0d0
      if(rmjd0_files[ipds] le rmjd0_files[ipds-1]) then begin
         massage,'PDS files in metafile have repeated or reversed start times'
         retall
      endif
      if(rmjd0_files[ipds] lt prev_end) then begin
         massage,'PDS files in metafile overlap in time'
         retall
      endif
   endfor

   rmjd0 = rmjd0_files[0]
endif

muflag = dummy(0)
i_vle = fix([dummy(17:18)],0)+1
;
;   nft is number of TIME points. Must be divided by 2 to get frequencies
;
nft       = nft/2

; MM
;
; MM 16 May 2026
; Decide how, if at all, the noise level is subtracted.
;
; poimode:
;   0 = no subtraction
;   1 = Zhang-style correction, old /POISSON behaviour
;   2 = estimate constant level from POISSON=[f1,f2]
;   3 = manually supplied constant level in input-file units
;
poimode = 0

if(keyword_set(manual_poi)) then begin
   poimode = poimode + 3
endif

if(keyword_set(poi)) then begin
   if(n_elements(poi) eq 2) then begin
      poimode = poimode + 2
   endif else begin
      poimode = poimode + 1
   endelse
endif

if(poimode gt 3) then begin
   print,'ERROR: Choose only one Poisson/noise option.'
   print,'Use /POISSON, POISSON=[f1,f2], or MANUAL_POI=p.'
   retall
endif

if(poimode eq 3) then begin
   manual_poi = float(manual_poi)
endif
; MM

;
pwr       = fltarr(nft)*0.0
power     = fltarr(nft)*0.0
power_err = fltarr(nft)*0.0

poilevel  = fltarr(nft)*0.0

; MM
;
; MM 16 May 2026
; Frequency grid. The last bin is the Nyquist bin. It is excluded
; when estimating a constant noise level from a frequency range,
; because the Nyquist bin has different statistics for a real time series.
;
frequency = (findgen(nft)+1.0) * T

n_frequency  = n_elements(frequency)
last_non_nyq = n_frequency - 2L

if(last_non_nyq lt 0L) then begin
   massage,'Frequency array too short to exclude Nyquist bin'
   retall
endif

if(poimode eq 2) then begin

   fpoi1_user = double(poi[0])
   fpoi2_user = double(poi[1])

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
;
; Loop over the trafos
;
n_selected       = 0l
itrafos          = 0l
tdead            = 1.0e-5
flux             = 0.0d0
fondo            = 0.0d0

; MM print,n_selected


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

itrafos = 0L
for ipds=0,n_pds_files-1 do begin

   if(metafile_flag eq 1) then begin
      ghats_openpds,pdsfiles[ipds],unit
   endif else begin
      ghats_openpds,pdsfiles[ipds],unit,/dialog
   endelse

   ; Skip header; it was already validated above.
   this_dummy               = bytarr(100)
   this_gh_version_string   = '                '
   this_observatory         = '                '
   this_instrument          = '                '
   this_target              = '                '
   this_rmjd0               = 0.0D0
   this_e                   = intarr(2)
   ghats_getheader,unit,this_gh_version_string,this_observatory, $
                        this_instrument,this_target,this_rmjd0, $
                        this_nft,this_T,this_ntrafos,this_e, $
                        this_proliferation,this_baryflag, $
                        this_n_spectral_bins,this_background_flag, $
                        this_dummy

   for local_trafo=0L,ntrafos_files[ipds]-1L do begin
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

      ;
      ; MM 16 May 2026
      ; BUG FIX:
      ; The old code overwrote POILEVEL for each selected PDS and then
      ; divided by N_SELECTED. Here POILEVEL is accumulated properly.
      ;
      case poimode of

         1: begin
            ; Old /POISSON behaviour: Zhang-style estimate for this PDS.
            poisson_estimate,poitmp,cnts,1.0/T,current_vle_rate,nft,fndet,i_vle,tdead
            poilevel = poilevel + poitmp
         end

         2: begin
            ; Estimate constant level from requested frequency range.
            ; The range excludes the Nyquist bin.
            poilevel = poilevel + (mean(pwr[wpoi]) + frequency*0.0)
         end

         3: begin
            ; User gives one constant level by hand.
            poilevel = poilevel + (manual_poi + frequency*0.0)
         end

         else: begin
            ; No subtraction requested.
         end

      endcase
; MM
   endif
   if(itrafos gt (lastindex-1)) then begin
      close,unit
      free_lun,unit
      goto,finito
   endif
   itrafos = itrafos + 1L
;
   endfor

   close,unit
   free_lun,unit

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
;
;  Now fill in frequency array
;
; MM 16 May 2026
; Frequency was already defined before the loop for POISSON range handling.
; frequency = (findgen(nft)+1.0) * T
;
   ps = ' spectra '
   if(n_selected eq 1) then ps = ' spectrum '
   print,'  ',strtrim(string(n_selected),1),' power'+ps+'selected'
;
;  If requested, subtract Poissonian component
;
; MM
   if(poimode gt 0) then begin
      power = power - poilevel
   endif
; MM

; MM   if(keyword_set(poi)) then begin
; MM      power = power - poilevel
; MM   endif
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

; MM
;
; MM 16 May 2026
; Tell the user which noise prescription was used.
;
case poimode of
   0: begin
      print,'WARNING: raw PDS returned; no noise level subtracted.'
      print,'Use /POISSON, POISSON=[f1,f2], or MANUAL_POI=p to subtract one.'
   end

   1: begin
      print,'Using /POISSON Zhang-style noise subtraction.'
   end

   2: begin
      print,'Using POISSON range: ',fpoi1,' - ',fpoi2,' Hz'
   end

   3: begin
      print,'Using MANUAL_POI.'
      print,'Assuming value is in input-file PDS units, before RMS conversion.'
      if(keyword_set(back)) then begin
         print,'If the constant was fitted from an RMS-normalised PDS,'
         print,'rerun GHX without RMS to measure it in input-file units:'
         print,"   ghx,'file.pds',nu,pow,powe"
      endif
   end
endcase
; MM
end
