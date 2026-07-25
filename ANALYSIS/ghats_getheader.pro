pro ghats_getheader,unit,gh_version_string,osserv,strumento,sorgente,rmjd, $
                     nft,T,ntotal_ffts,e,proliferation,baryflag,n_spectral_bins, $
                     background_flag,dummy,help=help

;+
; NAME: 
;      GHATS_GETHEADER
; PURPOSE: 
;      Reads in a header buffer from an opened PDS file
; EXPLANATION:
;      This procedure reads in a header buffer from an opened PDS file.
; CALLING SEQUENCE: 
;       GHATS_GETHEADER,UNIT,BUFFER
; INPUTS:
;       UNIT     = Unit for opened PDS file
;
; OUTPUTS:
;       NFT           = number of points per PDS
;       T             = inverse of light curve length
;       OBSERVATORY   = name of mission
;       TARGET        = name of target
;       RMJD0         = start date in RMJD
;       INSTRUMENT    = instrument used
;       E             = array with [start_channel,end_channel]
;       NTRAFOS       = number of PDS in file
;
; KEYWORDS:
;       HELP          = Print help and return
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       This is a helper routine normally called by higher-level readers.
; MODIFICATION HISTORY: 
;       T. Belloni  10 Nov 2001  implementation
;       T. Belloni  12 Mag 2002  adapted to MUFFT format
;       T. Belloni  11 Mag 2005  Mac version now supported
;       T. Belloni  22 Aug 2008  Mu 4.0
;	    T. Belloni  17 Mar 2009  Mu 6.0
;   	T. Belloni  11 Nov 2009  to Ghats
;       T. Belloni  15 Nov 2009  new ghats version
;		T. Belloni  02 Dec 2009  reflects new header v. 0.0.2
;       2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,'GHATS_GETHEADER'
   print,'Helper routine: normally called by GHATS readers, not run standalone.'
   print,'Calling sequence:'
   print,'  ghats_getheader, unit, gh_version_string, observatory, instrument, $'
   print,'                  target, rmjd, nft, T, ntotal_ffts, channels, $'
   print,'                  proliferation, baryflag, n_spectral_bins, $'
   print,'                  background_flag, dummy'
   print,'Reads the fixed GHATS binary header from an already-open file unit.'
   return
endif
;
;--------------------------------------------------------------------------
;
; Definitions
;
gh_version_string = '                '
nft               = 0L
T                 = 0.0d0
observatory       = '                '
target            = '                '
instrument        = '                '
rmjd0             = 0.0d0
ntrafos           = 0l
e                 = intarr(2)
proliferation     = 0
baryflag          = 0
n_spectral_bins   = 0
background_flag   = 0
;
; Read in header
;
readu,unit,gh_version_string,             $
                   osserv,                $
                   strumento,             $
                   sorgente,              $
                   rmjd,                  $
                   nft,                   $
                   T,                     $
                   ntotal_ffts,           $
                   e,                     $
                   proliferation,         $
                   baryflag,              $
                   n_spectral_bins,       $
                   background_flag,       $
                   dummy

info = strmid(gh_version_string,0,5)
if(info ne 'GHATS') then begin
    print,'Not a valid GHATS PDS file!'
    retall
endif

end
