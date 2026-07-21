pro gh_colors,filename,times,hr1,hr2,rate1,rate2,rate3,help=help
;+
; NAME: 
;      GH_COLORS
; PURPOSE: 
;      Get colors (hardnesses) from a PDS file
; EXPLANATION:
;      This procedure extracts hardness information from a PDS file 
;
; CALLING SEQUENCE: 
;       GH_COLORS,FILENAME,TIMES,HR1,HR2,RATE1,RATE2,RATE3[,/HELP]
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       TIMES    = Array with times for the HK series
;       HR1      = Array with soft color
;       HR2      = Array with hard color
;       RATE1    = Array with first rate curve
;       RATE2    = Array with second rate curve
;       RATE3    = Array with third rate curve
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
;       GHATS_OPENPDS:    Opens a PDS file
;       GHATS_GETHEADER:  Reads in the header
;       GH_READ_COLORS:   Rein in color information
; NOTES:
;       NONE
; MODIFICATION HISTORY: 
;       T. Belloni  11 Dec 2003  from muxana_hk
;		T. Belloni  25 Jan 2012  from muxana_colors
;		T. Belloni  20 Feb 2012  no floor rounding for times array
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_COLORS'
   print,''
   print,'Extract color/hardness light curves from a GHATS PDS file.'
   print,''
   print,'Usage:'
   print,"  GH_COLORS, 'file.pds', times, hr1, hr2, rate1, rate2, rate3"
   print,''
   print,'Outputs:'
   print,'  times        Time array, in seconds from the file start'
   print,'  rate1-3      Three stored band count-rate curves'
   print,'  hr1,hr2      Hardness ratios rate2/rate1 and rate3/rate1'
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
muflag    = dummy(0)
i_vle = fix([dummy(17:18)],0)+1
;
rate1   = fltarr(ntrafos)
rate2   = fltarr(ntrafos)
rate3   = fltarr(ntrafos)
hr1     = fltarr(ntrafos)
hr2     = fltarr(ntrafos)
times   = fltarr(ntrafos)

nft = nft / 2
;
; Loop over the trafos
;
itrafos          = 0l

for itrafos=0l,ntrafos-1l do begin
;
   gh_read_colors,unit,nft,rmjd,r1,r2,r3
;
   times(itrafos) = ((rmjd-rmjd0)*86400.0)
   rate1(itrafos) = r1
   rate2(itrafos) = r2
   rate3(itrafos) = r3
   hr1(itrafos)   = r2 / r1
   hr2(itrafos)   = r3 / r1
endfor
;
free_lun,unit
;
end
