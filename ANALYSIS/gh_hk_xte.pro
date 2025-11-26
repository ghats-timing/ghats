pro gh_hk_xte,filename,times,vle,ndet,poiss
;+
; NAME: 
;      GH_HK_XTE
; PURPOSE: 
;      Get auxiliary data from a PDS file
; EXPLANATION:
;      This procedure extracts from a PDS file for XTE/PCA all HK information, that is
;      time series of: vle rate, number of detectors on, and estimate of
;      the poissonian level
;
; CALLING SEQUENCE: 
;       GH_HK_XTE,FILENAME,TIMES,VLE,NDET,POISS
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       TIMES    = Array with times for the HK series
;       VLE      = Array with VLE rate
;       NDET     = Array with number of detectors on
;       POISS    = Array with computed Poissonian level
;
; KEYWORDS:
;       NONE
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       MUXANA_OPENPDS:    Opens a PDS file
;       MUXANA_GETHEADER:  Reads in the header
;       MUXANA_SHEADER:    Parses the header
; NOTES:
;       NONE
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  10 Nov 2001  modularization
;       T. Belloni  12 Mag 2002  adapted to MUFFT
;       T. Belloni  11 Dec 2003  Mu3 version (colors)
;		T. Belloni 	01 Dec 2010  from mu6
;		T. Belloni  22 Nov 2011  free_lun
;		T. Belloni  03 Jan 2012  fixed getheader bug
;		T. Belloni  20 Feb 2012  no floor rounding for times array
;		T. Belloni  22 May 2012  changed name to gh_hk_xte
;-
;--------------------------------------------------------------------------
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
;
vle     = fltarr(ntrafos)
ndet    = fltarr(ntrafos)
times   = fltarr(ntrafos)
poiss   = fltarr(ntrafos)
;
; Loop over the trafos
;
;
;   nft is number of TIME points. Must be divided by 2 to get frequencies
;
nft       = nft/2
pwr       = fltarr(nft)*0.0
itrafos          = 0l

for itrafos=0l,ntrafos-1l do begin
;
   read_pds_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,a0,pwr
;
   times(itrafos) = ((rmjd-rmjd0)*86400.0)
   vle(itrafos)   = current_vle_rate / fndet
   ndet(itrafos)  = fndet
   poiss(itrafos) = poisson
endfor
;
free_lun,unit
;
end
