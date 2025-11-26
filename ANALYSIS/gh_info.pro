pro gh_info,filename
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
;       GH_INFO,FILENAME
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       NONE
;
; EXAMPLE:
;       Obtain basic info from an RXTE/PCA file for GRS 1915+105:
;
;       MU> GH_INFO,'prova.pds'
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
;-
;--------------------------------------------------------------------------
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
