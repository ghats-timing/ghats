pro gh_info,filename,help=help
; MM pro gh_info,filename
;+
; NAME: 
;      GH_INFO
; PURPOSE: 
;      Print to terminal basic information on a GHATS PDS or FFT file
; EXPLANATION:
;      This procedure reads the header of a GHATS PDS or FFT file and reports
;      the basic information about it, that is: observatory, instrument,
;      target name, number of trafos, duration of each interval, number of
;      powers or Fourier bins, start-end channels, and starting MJD.
;
; CALLING SEQUENCE: 
;       GH_INFO,FILENAME ; MM: Changed for the ones below
; CALLING SEQUENCE: 
;       GH_INFO,FILENAME
;       GH_INFO,/HELP
; INPUTS:
;       FILENAME = name of the input PDS or FFT file
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       NONE ; MM: Replaced with the one below
;       HELP     = If set, print usage information and return
;
; EXAMPLE:
;       Obtain basic info from an RXTE/PCA file for GRS 1915+105:
;
;       MU> GH_INFO,'prova.pds'
;       MU> GH_INFO,/HELP
;
; muxana_info ------------------------------------------------------
;       PDS file name   :   prova.pds
;
;       Observatory     :   XTE
;       Instrument      :   PCA
;       Target          :   GRS1915+105
;
;       N. of trafos    :   106
;       Trafo duration  :   16.0000
;       N. of powers    :   1024
;       Starting mjd    :   50909.879
;
; ------------------------------------------------------------------
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       MUXANA_GETHEADER: Get information from PDS header
; NOTES:
;       The procedure produces no output besides the printout. To have
;       information usable from the program level, use muxana_header
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  17 May 2002  adapted to MUFFT
;       T. Belloni  16 Jun 2005  close file added
;		T. Belloni  06 May 2010  from muana_info
;		T. Belloni  09 Mar 2012  added TURBO, WINDOW, full BANDS
;		T. Belloni  22 Jan 2017  corrected Time Step computation
;		T. Belloni  17 Oct 2017  added warning for separate channel selections
;		T. Belloni  11 Nov 2018  added more precision for MJD0
;       M. Mendez  14 May 2026  added /HELP and Nyquist-frequency output
;                                with assistance from ChatGPT
;-
;--------------------------------------------------------------------------
; MM
;
; Print concise usage information and return.
;
if(keyword_set(help)) then begin
   print,''
   print,'GH_INFO'
   print,''
   print,'Usage:'
   print,"  GH_INFO,'file.pds'"
   print,'  GH_INFO,/HELP'
   print,''
   print,'Purpose:'
   print,'  Print basic information from a GHATS PDS/FFT file header.'
   print,''
   return
endif
; MM
;
; Identify product type from the extension. The binary header is common to PDS
; and FFT products, but using the matching opener gives the right dialog filter
; and keeps the printed labels honest.
;
filetype = 'GHATS'
if(n_elements(filename) ne 0) then begin
   fname_low = strlowcase(strtrim(filename,2))
   if(strlen(fname_low) ge 4) then begin
      ext = strmid(fname_low,strlen(fname_low)-4,4)
      if(ext eq '.pds') then filetype = 'PDS'
      if(ext eq '.fft') then filetype = 'FFT'
   endif
endif

case filetype of
   'FFT': ghats_openfft,filename,unit,/dialog
   else:  ghats_openpds,filename,unit,/dialog
endcase
;
; Read in GHATS product header
;
dummy     = bytarr(100)
gh_version_string = '                '
observatory       = '                '
instrument        = '                '
target            = '                '
nft               = 0L
ntrafos           = 0L
rmjd0             = 0.0D0

ghats_getheader,unit,gh_version_string,observatory,instrument,target,rmjd0, $
				                     nft,T,ntrafos,e,proliferation,baryflag,n_spectral_bins, $
				                     background_flag,dummy
; MM
;
; Detect obvious header-layout problems before computing derived quantities.
; This avoids misleading infinities if a file was written without the frequency
; resolution field used by the current GHATS header. Do not fall back merely
; because newer auxiliary fields are absent or zero: the main header can still
; contain a valid T, ntrafos, and channel range.
;
header_ok = 1
legacy_header = 0
short_nft_header = 0
aux_header_ok = 1
if(nft le 0L) then header_ok = 0
if(T le 1.0d-100) then header_ok = 0
if(proliferation le 0) then begin
   aux_header_ok = 0
   proliferation = 1
endif

if(header_ok eq 0) then begin
   ;
   ; Command-line gh_xte historically could write NFT as a 16-bit INT when
   ; NPDS was passed as an integer literal. That shifts the rest of the
   ; header by two bytes. Try this specific recovery before the older fallback.
   ;
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
      legacy_header = 0
      short_nft_header = 1
      aux_header_ok = 1
      if(proliferation le 0) then begin
         aux_header_ok = 0
         proliferation = 1
      endif
   endif
endif

if(header_ok eq 0) then begin
   point_lun,unit,0
   dummy = bytarr(100)
   gh_version_string = '                '
   observatory       = '                '
   instrument        = '                '
   target            = '                '
   rmjd0             = 0.0D0
   nft               = 0L
   ntrafos           = 0L
   e                 = intarr(2)
   T                 = 0.0D0
   proliferation     = 1
   baryflag          = -1
   n_spectral_bins   = -1
   background_flag   = -1
   readu,unit,gh_version_string,observatory,instrument,target,rmjd0, $
        nft,ntrafos,e
   legacy_header = 1
endif
; MM
;
; The header stores T as the frequency resolution, i.e. T = 1 / duration.
; The written frequency bins run from 1 to nft/2, so the highest written
; frequency is the Nyquist frequency:
;
;        nu_nyq = (nft/2) * T
;
if(legacy_header eq 0) then begin
   nyquist = 0.5d0 * double(nft) * double(T)
   duration = 1.0d0 / double(T)
   timestep = duration / double(proliferation)
endif
;
; MM
;
print,''
print,'GH_info ------------------------------------------------'
print,''
print,filetype+' file name   :   '+filename
print,''
print,'Observatory     :   '+observatory
print,'Instrument      :   '+instrument
print,'Target          :   '+target
print,''
if(filetype eq 'FFT') then begin
   print,'N. of FFTs      :   ',strtrim(string(ntrafos),1)
   if(legacy_header eq 0) then begin
      print,'FFT duration    :   ',strtrim(string(duration),1)
      print,'N. of FFT bins  :   ',strtrim(string(nft/2),1)
   endif else begin
      print,'N. of time bins :   ',strtrim(string(nft),1)
   endelse
endif else begin
   print,'N. of trafos    :   ',strtrim(string(ntrafos),1)
   if(legacy_header eq 0) then begin
      print,'Trafo duration  :   ',strtrim(string(duration),1)
   endif
   print,'N. of powers    :   ',strtrim(string(nft),1)
endelse
; MM
if(legacy_header eq 0) then begin
   print,'Nyquist freq.   :   ',strtrim(string(nyquist),1),' Hz'
   if(short_nft_header eq 1) then begin
      print,'WARNING: header was written with NFT as a 16-bit integer.'
      print,'         Regenerate this FFT/PDS file before using it for analysis.'
   endif
   if((e[0] ge 0) and (e[1] ge 0) and (e[0] gt e[1])) then begin
      print,'WARNING: channel range in header is inconsistent.'
   endif
   if(aux_header_ok eq 0) then begin
      print,'WARNING: auxiliary header fields after channel range are absent'
      print,'         or inconsistent; assuming proliferation/time step = 1.'
   endif
endif else begin
   print,'WARNING: header lacks the current frequency-resolution field.'
   print,'         Frequency resolution, duration, Nyquist frequency,'
   print,'         barycenter/background flags, and spectral-bin count'
   print,'         are not available from this header.'
endelse
; MM
if(legacy_header eq 0) then begin
   print,'Start-end chans.:   ',strtrim(string(abs(e[0])),1)+'-'+  $
   	                            strtrim(string(abs(e[1])),1)
   IF (e(0) lt 0) THEN BEGIN
   	print,'WARNING: separate channel selection for different instruments!',''  ; GIGUS
   ENDIF
endif else begin
   print,'Start-end chans.:   unavailable'
endelse

print,'Starting mjd    :   ',strtrim(string(rmjd0,FORMAT='(f20.10)'),1)
if(legacy_header eq 0) then begin
   print,'Time step       :   ',strtrim(string(timestep),1)
endif

if(legacy_header eq 0) then begin
   if(baryflag eq 1) then begin
   	print,'Barycentered    :   YES'
   endif else begin
   	print,'Barycentered    :   NO'
   endelse
   print,'Spectral bins   :   ',strtrim(string(n_spectral_bins),1)
   if(background_flag eq 1) then begin
   	print,'Background      :   YES'
   endif else begin
   	print,'Background      :   NO'
   endelse
   pari = (dummy[0] mod 2)
   if (pari eq 1) then begin
   	print,'Turbo           :   YES'
   endif else begin
   	print,'Turbo           :   NO'
   endelse
endif

print,''
print,'------------------------------------------------------------'
print,''

close,unit

end
