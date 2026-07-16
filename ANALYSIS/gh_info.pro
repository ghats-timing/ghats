pro gh_info,filename,help=help
; MM pro gh_info,filename
;+
; NAME: 
;      GH_INFO
; PURPOSE: 
;      Print to terminal basic information on a PDS file
; EXPLANATION:
;      This procedure reads the header of a PDS file and reports
;      the basic information about it, that is: observatory, instrument,
;      target name, number of trafos, duration of each interval, number of
;      powers in each PDS, start-end channels, and starting MJD.
;
; CALLING SEQUENCE: 
;       GH_INFO,FILENAME ; MM: Changed for the ones below
; CALLING SEQUENCE: 
;       GH_INFO,FILENAME
;       GH_INFO,/HELP
; INPUTS:
;       FILENAME = name of the input PDS file
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
; Open pds file
;
ghats_openpds,filename,unit,/dialog
;
; Read in PDS file header
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
; The header stores T as the frequency resolution, i.e. T = 1 / duration.
; The written frequency bins run from 1 to nft/2, so the highest written
; frequency is the Nyquist frequency:
;
;        nu_nyq = (nft/2) * T
;
nyquist = 0.5d0 * double(nft) * double(T)
;
; MM
;
print,''
print,'GH_info ------------------------------------------------'
print,''
print,'PDS file name   :   '+filename
print,''
print,'Observatory     :   '+observatory
print,'Instrument      :   '+instrument
print,'Target          :   '+target
print,''
print,'N. of trafos    :   ',strtrim(string(ntrafos),1)
print,'Trafo duration  :   ',strtrim(string(1.0/T),1)
print,'N. of powers    :   ',strtrim(string(nft),1)
; MM
print,'Nyquist freq.   :   ',strtrim(string(nyquist),1),' Hz'
; MM
print,'Start-end chans.:   ',strtrim(string(abs(e[0])),1)+'-'+  $
	                            strtrim(string(abs(e[1])),1)
IF (e(0) lt 0) THEN BEGIN
	print,'WARNING: separate channel selection for different instruments!',''  ; GIGUS
ENDIF

print,'Starting mjd    :   ',strtrim(string(rmjd0,FORMAT='(f20.10)'),1)
print,'Time step       :   ',strtrim(string(1.0/T/proliferation),1)

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

print,''
print,'------------------------------------------------------------'
print,''

close,unit

end
