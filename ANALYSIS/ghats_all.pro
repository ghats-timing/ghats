pro ghats_all,filename,ps=psopt,poisson=poisso,rebin=rebin,help=help
; MM pro ghats_all,filename,ps=psopt,poisson=poisso,rebin=rebin
;+
; NAME:
;      GHATS_ALL
; PURPOSE:
;      Produces and plots light curve and average PDS from a pds file
; EXPLANATION:
;      This procedure is intended as a first check of a PDS file.
;      It plots the light curve and an average PDS, but does not
;      produce any usable data.
;
; CALLING SEQUENCE:
;       GHATS_ALL[,FILENAME][,/PS] ; MM: chnaged for the one below
;       GHATS_ALL[,FILENAME][,/PS][,/POISSON][,REBIN=REBIN][,/HELP]
; INPUTS:
;       FILENAME = name of the input PDS file. If filename is not
;                  defined, a dialog window will appear
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       PS       = If set, output goes to a PS file
;       POISSON  = Automatic Poissonian subtraction
;	REBIN    = flag for rebinnin factor
;       HELP     = If set, print usage information and return
;
; EXAMPLE:
;       NONE ; MM: Changed for the ones below
; EXAMPLES:
;       GHATS_ALL,'file.pds'
;       GHATS_ALL,'file.pds',/POISSON [MM: this probably only works for RXTE]
;       GHATS_ALL,'file.pds',/PS,REBIN=-100
;       GHATS_ALL,/HELP
;
; COMMON BLOCKS:
;       None
; ROUTINES USED:
;       MUXANA_INFO: Get basic information from PDS file
;       MUXANA_LICU: Extracts light curve
;       GHX        : Extracts total average PDS
;       GHREBIN    : Rebins the PDS
;       GH_PLOT_LICU : Plots the light curve
;       GH_PLOT_POWER : Plots the PDS
; NOTES:
;       The procedure produces no output besides the plots.
; MODIFICATION HISTORY:
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni   9 Nov 2001  clean commented version + mu_plot_licu
;       T. Belloni   9 Nov 2001  dialog window added
;       T. Belloni  20 Nov 2001  complex one-window output
;       T. Belloni  21 Nov 2001  PS option
;       T. Belloni  17 May 2002  Adapted to MUFFT
;       T. Belloni  16 Jul 2002  Poissonian subtraction added
;       T. Belloni  15 Nov 2003  removed Unix-like syntaxes
;       T. Belloni  03 Apr 2009  plot_oo removed
;	    T. Belloni  05 Apr 2009  gmu version adapted
;		T. Belloni  28 Apr 2010  weird bug solved by switching to mux 
;		T. Belloni  06 May 2010  from muxana_all 
;		T. Belloni  07 May 2010  added rebin flag
;		T. Belloni  02 Jan 2012  removed font=0 for GDL compatibility
;		T. Belloni  25 Jan 2012  fixed time per interval output bug
;		T. Belloni  17 Oct 2017  added abs() to channels to intercept negative values
;       M. Mendez  14 May 2026  added /HELP, prints Nyquist frequency in the plot, updated the examples; with assistance from ChatGPT
;-
;--------------------------------------------------------------------------
common sis,sistema
common vers,versione, data_versione
; MM
;
; Print concise usage information and return.
;
if(keyword_set(help)) then begin
   print,''
   print,'GHATS_ALL'
   print,''
   print,'Usage:'
   print,"  GHATS_ALL,'file.pds'"
   print,"  GHATS_ALL,'file.fft'  ; PDS panel only; ignore light curve"
   print,"  GHATS_ALL,'file.pds',/POISSON [MM: this probably only works for RXTE]"
   print,"  GHATS_ALL,'file.pds',/PS,REBIN=-100"
   print,''
   print,'Notes:'
   print,'  For FFT input, the PDS panel is computed from the FFT records.'
   print,'  The light-curve panel is not reliable until GH_LICU supports FFT records.'
   print,''
   print,'Keywords:'
   print,'  /PS       Write plot to ghats.ps'
   print,'  /POISSON  Subtract estimated Poisson level'
   print,'  REBIN=    Rebin factor passed to GH_REB'
   print,'  /HELP     Print this message'
   print,''
   return
endif
; MM
;
; Identify product type from the extension. GHATS_ALL remains a PDS quick-look
; tool, but for FFT files the PDS panel must be computed from FFT records, not
; through GHX/read_pds_line.
;
filetype = 'PDS'
if(n_elements(filename) ne 0) then begin
   fname_low = strlowcase(strtrim(filename,2))
   if(strlen(fname_low) ge 4) then begin
      ext = strmid(fname_low,strlen(fname_low)-4,4)
      if(ext eq '.fft') then filetype = 'FFT'
   endif
endif
; MM
;
; Open tra file
;
if(filetype eq 'FFT') then begin
   ghats_openfft,filename,unit,/dialog
endif else begin
   ghats_openpds,filename,unit,/dialog
endelse
ntrafos   = 0l
dummy     = bytarr(100)
gh_version_string = '                '
observatory       = '                '
instrument        = '                '
target            = '                '
rmjd0             = 0.0D0

ghats_getheader,unit,gh_version_string,observatory,instrument,target,rmjd0, $
				                     nft,T,ntrafos,e,proliferation,baryflag,n_spectral_bins, $
				                     background_flag,dummy
; MM
;
; Command-line gh_xte versions before the np=long(np) fix could write NFT as
; a 16-bit integer. That shifts the remaining header fields by two bytes. Use
; the same diagnostic recovery as GH_INFO so the printed/graphic header summary
; does not silently display nonsense for those files.
;
header_ok = 1
short_nft_header = 0
if(nft le 0L) then header_ok = 0
if(T le 1.0d-100) then header_ok = 0

if(header_ok eq 0) then begin
   point_lun,unit,72L
   nft_short = 0
   T_short = 0.0D0
   ntrafos_short = 0L
   e_short = intarr(2)
   proliferation_short = 0
   baryflag_short = 0
   n_spectral_bins_short = 0
   background_flag_short = 0
   dummy_short = bytarr(100)
   readu,unit,nft_short,T_short,ntrafos_short,e_short, $
        proliferation_short,baryflag_short,n_spectral_bins_short, $
        background_flag_short,dummy_short
   if((nft_short gt 0) and (T_short gt 1.0d-100)) then begin
      nft = long(nft_short)
      T = T_short
      ntrafos = ntrafos_short
      e = e_short
      proliferation = proliferation_short
      baryflag = baryflag_short
      n_spectral_bins = n_spectral_bins_short
      background_flag = background_flag_short
      dummy = dummy_short
      header_ok = 1
      short_nft_header = 1
   endif
endif
close,unit
; MM
;
; The header stores the frequency resolution T = 1 / segment_length,
; not the Nyquist frequency explicitly. Since the stored powers run from
; frequency bin 1 to nft/2, the Nyquist frequency is:
;
;        nu_nyq = (nft/2) * T
;
if(header_ok eq 1) then begin
   nyquist = 0.5d0 * double(nft) * double(T)
endif else begin
   nyquist = 0.0d0
endelse
; MM
;
;  PS output setup
;
csi = 1.3
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='ghats.ps'
   csi=0.8
endif
;
; Prints out basic info
;
gh_info,filename
;
observatory = strtrim(observatory,2)
instrument  = strtrim(instrument,2)
target      = strtrim(target,2)
;
; Extract the light curve
;
if(filetype eq 'FFT') then begin
   print,'GHATS_ALL WARNING: FFT input: light-curve panel is not reliable; ignore it.'
endif
gh_licu,filename,times,licu,deltat=t,mjd=rmjd0
tmin        = times(0)
tmax        = times(n_elements(times)-1)+1.0/t
;
; Plot light curve
;
;window,0,TITLE=observatory+'/'+instrument+'  '+target
plot,times,licu,psym=10,xstyle=1,ystyle=1,             $
     xtitle='Time (s)',ytitle='Rate (cts/s)',position=[0.15,0.7,0.95,0.95 ]
if(filetype eq 'FFT') then begin
   xyouts,0.16,0.64,'FFT input: light-curve panel is not reliable', $
          charsize=0.8*csi,/normal
endif
;
; Extract power spectrum from full dataset and rebin it at -100
;
if(filetype eq 'FFT') then begin
   ;
   ; Compute the quick-look PDS directly from FFT records. This follows the
   ; existing cross-spectrum routines: PDS = abs(FFT)^2 * 2 / counts, averaged
   ; over selected transforms. Do not use GHX_FFT here, because it averages
   ; complex Fourier amplitudes rather than powers.
   ;
   ghats_openfft,filename,unit,/dialog

   header_dummy               = bytarr(100)
   header_gh_version_string   = '                '
   header_observatory         = '                '
   header_instrument          = '                '
   header_target              = '                '
   header_rmjd0               = 0.0D0
   header_ntrafos             = 0L
   ghats_getheader,unit,header_gh_version_string,header_observatory, $
                    header_instrument,header_target,header_rmjd0, $
                    header_nft,header_T,header_ntrafos,header_e, $
                    header_proliferation,header_baryflag, $
                    header_n_spectral_bins,header_background_flag, $
                    header_dummy

   if(short_nft_header eq 1) then begin
      ;
      ; Broken pre-fix command-line FFT headers are two bytes shorter because
      ; NFT was written as an INT. Start reading records at the actual data
      ; offset for that specific layout.
      ;
      point_lun,unit,198L
      header_nft = nft
      header_T = T
      header_ntrafos = ntrafos
      header_dummy = dummy
   endif

   muflag = header_dummy(0)
   i_vle = fix([header_dummy(17:18)],0)+1
   nfreq = long(header_nft)/2L
   rdata = complexarr(nfreq)*0.0
   frequency = (findgen(nfreq)+1.0) * header_T
   power = fltarr(nfreq)*0.0
   power_err = fltarr(nfreq)*0.0
   poilevel = fltarr(nfreq)*0.0
   n_selected = 0L
   tdead = 1.0e-5

   for itrafos=0L,header_ntrafos-1L do begin
      read_fft_line,unit,muflag,rmjd_fft,cnts_fft,poisson_fft, $
                    current_vle_rate_fft,fndet_fft,rdata
      pwr_fft = abs(rdata)^2 * 2.0/cnts_fft
      power = power + pwr_fft
      n_selected = n_selected + 1L

      if(keyword_set(poisso)) then begin
         poisson_estimate,poitmp,cnts_fft,1.0/header_T, $
                          current_vle_rate_fft,nfreq,fndet_fft,i_vle,tdead
         poilevel = poilevel + poitmp
      endif
   endfor

   if(n_selected gt 0L) then begin
      power = power / n_selected
      power_err = power / sqrt(n_selected)
      if(keyword_set(poisso)) then begin
         poilevel = poilevel / n_selected
         power = power - poilevel
      endif
   endif else begin
      print,'Warning! No FFTs retrieved for quick-look PDS!'
   endelse

   close,unit
   free_lun,unit
   print,'  ',strtrim(string(n_selected),1),' FFTs selected for PDS'
endif else begin
   if(keyword_set(poisso)) then begin
      ;muxana_n,filename,frequency,power,power_err,1L,100000L,/poisson
      ghx,filename,frequency,power,power_err,/poisson
     endif else begin
      ;muxana_n,filename,frequency,power,power_err,1L,100000L
      ghx,filename,frequency,power,power_err
   endelse
endelse

if(keyword_set(rebin)) then begin
	rebin_factor = rebin
endif else begin
	rebin_factor = -100
endelse
if(rebin_factor ne 1) then begin
	gh_reb,frequency,power,power_err,rebin_factor,nu,pow,powerr
endif else begin
	nu     = frequency
	pow    = power
endelse
;murebin,frequency,frequency,power,power_err,-100,nu,nuerr,pow,powerr,nrd
;
; Plot power spectrum
;
buoni = where(pow gt 0.0, nbuoni)

if(nbuoni gt 0) then begin
   plot,nu(buoni),pow(buoni),psym=10,xstyle=1,ystyle=1,     $
        xtitle='Frequency (Hz)',ytitle='Power',/noerase,  $
        position=[0.50,0.14,0.95,0.60],/xlog,/ylog

   if(rebin_factor ne 1) then begin
   	if(sistema eq 'IDL') then begin
      		muerrplot,nu,pow-powerr,pow+powerr,width=0.0
   	endif else begin
   		muploterr,nu,pow,powerr,/xlog,/ylog,psym=3
   	endelse
   endif
endif else begin
   print,'GHATS_ALL WARNING: no positive powers available for log PDS plot.'
   plot,[1.0,10.0],[1.0,10.0],xstyle=1,ystyle=1,/nodata,/noerase, $
        xtitle='Frequency (Hz)',ytitle='Power', $
        position=[0.50,0.14,0.95,0.60],/xlog,/ylog
   xyouts,0.60,0.36,'No positive powers to plot',charsize=csi,/normal
endelse
;
; Label information
;
xyouts,0.05,0.57,versione,charsize=csi+0.2,/normal;,font=0
xyouts,0.05,0.50,'Observatory   ',charsize=csi,/normal;,font=0
xyouts,0.05,0.47,'Instrument    ',charsize=csi,/normal;,font=0
xyouts,0.05,0.44,'Target        ',charsize=csi,/normal;,font=0
xyouts,0.05,0.41,'N. of trafos  ',charsize=csi,/normal;,font=0
xyouts,0.05,0.38,'Channels      ',charsize=csi,/normal;,font=0
xyouts,0.05,0.35,'Trafo duration',charsize=csi,/normal
xyouts,0.05,0.32,'N. of powers  ',charsize=csi,/normal
;MM
xyouts,0.05,0.29,'Nyquist freq. ',charsize=csi,/normal
xyouts,0.05,0.26,'Starting MJD  ',charsize=csi,/normal
xyouts,0.05,0.23,'Obs. Date     ',charsize=csi,/normal
xyouts,0.05,0.17,'Filename      ',charsize=csi,/normal
xyouts,0.05,0.14,'User          ',charsize=csi,/normal
xyouts,0.05,0.11,'Date          ',charsize=csi,/normal
xyouts,0.05,0.08,'Time          ',charsize=csi,/normal

; MM xyouts,0.05,0.32,'N. of powers  ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.29,'Starting MJD  ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.26,'Obs. Date     ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.20,'Filename      ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.17,'User          ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.14,'Date          ',charsize=csi,/normal;,font=0
; MM xyouts,0.05,0.11,'Time          ',charsize=csi,/normal;,font=0
;
xyouts,0.20,0.50,observatory     ,charsize=csi,/normal;,font=0
xyouts,0.20,0.47,instrument      ,charsize=csi,/normal;,font=0
xyouts,0.20,0.44,target          ,charsize=csi,/normal;,font=0
xyouts,0.20,0.41,strtrim(string(ntrafos),2) ,charsize=csi,/normal;,font=0
xyouts,0.20,0.38,strtrim(string(abs(e[0])),1)+'-'+strtrim(string(abs(e[1])),1), $
                                             charsize=csi,/normal;,font=0
; MM
xyouts,0.20,0.35,strtrim(string(t),2)+' s' ,charsize=csi,/normal
xyouts,0.20,0.32,strtrim(string(nft),2) ,charsize=csi,/normal
xyouts,0.20,0.29,strtrim(string(nyquist),2)+' Hz' ,charsize=csi,/normal
xyouts,0.20,0.26,strtrim(string(rmjd0),2) ,charsize=csi,/normal
jd=rmjd0+2400000.5d
daycnv,jd,anno,mese,giorno,ora
xyouts,0.20,0.23,strtrim(string(giorno),2)+'/'+   $
                 strtrim(string(mese),2)+'/'+     $
                 strtrim(string(anno),2),charsize=csi,/normal


; MM xyouts,0.20,0.35,strtrim(string(t),2)+' s' ,charsize=csi,/normal;,font=0
; MM xyouts,0.20,0.32,strtrim(string(nft),2) ,charsize=csi,/normal;,font=0
; MM xyouts,0.20,0.29,strtrim(string(rmjd0),2) ,charsize=csi,/normal;,font=0
; MM jd=rmjd0+2400000.5d
; MM daycnv,jd,anno,mese,giorno,ora
; MM xyouts,0.20,0.26,strtrim(string(giorno),2)+'/'+   $
; MM                  strtrim(string(mese),2)+'/'+     $
; MM                  strtrim(string(anno),2),charsize=csi,/normal;,font=0
; MM
;
opsys = !version.os_family
if (opsys eq 'Windows') then begin
    sla = '\'
   endif else begin
    sla = '/'
endelse
ii = strpos(filename,sla,/reverse_search)
filen = strmid(filename,ii+1)
; MM xyouts,0.20,0.20,strtrim(filen,2) ,charsize=csi,/normal;,font=0
opsys = !version.os_family
if (opsys eq 'Windows') then begin
    io = 'Win user'
   endif else begin
    spawn,'whoami',io
endelse
; MM
xyouts,0.20,0.17,strtrim(filen,2) ,charsize=csi,/normal
xyouts,0.20,0.14,strtrim(io[0],2) ,charsize=csi,/normal
s=strsplit(strtrim(systime(0),2),' ',/extract)
xyouts,0.20,0.11,s[0]+' '+s[1]+' '+s[2]+' '+s[4],charsize=csi,/normal
xyouts,0.20,0.08,s[3],charsize=csi,/normal
if(short_nft_header eq 1) then begin
   xyouts,0.05,0.02,'WARNING: 16-bit NFT header; regenerate file', $
          charsize=0.8*csi,/normal
endif
; MM
; MM xyouts,0.20,0.17,strtrim(io[0],2) ,charsize=csi,/normal;,font=0
; MM s=strsplit(strtrim(systime(0),2),' ',/extract)
; MM 
; MM xyouts,0.20,0.14,s[0]+' '+s[1]+' '+s[2]+' '+s[4],charsize=csi,/normal;,font=0
; MM xyouts,0.20,0.11,s[3],charsize=csi,/normal;,font=0
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif

close,/all
end
