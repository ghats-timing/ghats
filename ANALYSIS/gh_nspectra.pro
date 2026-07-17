pro gh_nspectra,filename,ntrafos,help=help
;+
; NAME: 
;      GH_NSPECTRA
; PURPOSE: 
;      Computes the number of trafos in a PDS file
; EXPLANATION:
;      This procedure computes the number of trafos contained in a 
;      PDS file. 
;
; CALLING SEQUENCE: 
;       GH_NSPECTRA,filename,ntrafos[,/HELP]
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       NTRAFOS  = Output number of trafos
;
; KEYWORDS:
;       HELP     = If set, print usage information and return
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       The procedure works only for files with only ONE header.
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  17 May 2002  adapted to MUFFT
;		T. Belloni  06 May 2010  from MU
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_NSPECTRA'
   print,''
   print,'Read the number of transforms stored in a GHATS PDS file.'
   print,''
   print,'Usage:'
   print,"  GH_NSPECTRA, 'file.pds', ntrafos"
   print,''
   print,'Output:'
   print,'  ntrafos  Number of spectra/transforms in the file header'
   print,''
   return
endif
;
; Open pds file
;
ghats_openpds,filename,unit,/dialog
;
; Read in pds file header
;
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

end
