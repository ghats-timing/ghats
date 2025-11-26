pro gh_licu,filename,times,licu,deltat=tt,mjd=rmjd0;,back=back
;+
; NAME: 
;      GH_LICU
; PURPOSE: 
;      Reads in a light curve from a PDS file
; EXPLANATION:
;      This procedure reads the licu info of a PDS file.
;
; CALLING SEQUENCE: 
;       GH_LICU,FILENAME,TIMES,LICU[,DELTAT=T][,MJD=RMJD0]
; INPUTS:
;       FILENAME = name of the input PDS file. If filename is empty
;                  or contains a non-existent file, a window will
;                  prompt for a file, which will be put in this parameter
;
; OUTPUTS:
;       TIMES    = array of time values for licu
;       LICU     = array of rates for licu
;
; KEYWORDS:
;       T        = time bin size for the output light curve
;       MJD      = start modified julian date
;
; EXAMPLE:
;       Produce and plot a light curve from a PDS file:
;
;       MU> GH_LICU,'prova.pds',times,licu,deltat=deltat
;       MU> GH_PLOT_LICU,times,licu
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       GHATS_OPENPDS  : Opens PDS file
;       GHATS_GETHEADER: Reads PDS header and reports number of trafos
;       READ_PDS_LINE   : Reads in a single power spectrum from the file
; NOTES:
;       NONE
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni   9 Nov 2001  keyword version + muxana_sheader/muxana_opentra
;       T. Belloni  10 Nov 2001  muxana_getheader
;       T. Belloni  17 May 2002  adapted to MUFFT
;       T. Belloni  11 Dec 2003  Mu3 version (colors)
;		T. Belloni  06 May 2010  from mu_plot_licu
;		T. Belloni  22 Nov 2011  free_lun
;		T. Belloni  04 Jan 2012  tt introduced
;		T. Belloni  13 Nov 2017  bacground option for Astrosat/CZTI - currently disabled
;		T. Belloni  04 Feb 2019  double forced on time,licu

;			;-
;--------------------------------------------------------------------------
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

muflag    = dummy(0)
;
times=dblarr(ntrafos)
licu =dblarr(ntrafos)
;fondo=fltarr(ntrafos)
;
;   nft is number of TIME points. Must be divided by 2 to get frequencies
;
nft       = nft/2
;
pwr       = fltarr(nft)*0.0
;
; Loop over the trafos
;
itrafos          = 0l

for itrafos=0l,ntrafos-1 do begin
;
   read_pds_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,a0,pwr

   times(itrafos)  = (rmjd-rmjd0)*86400.0
   licu(itrafos)   = cnts*T
;   fondo(itrafos)  = current_vle_rate*T
;
endfor
free_lun,unit

tt = 1.0/t

;if(keyword_set(back)) THEN BEGIN   ; for Astrosat/CZTI, vle_rate contains the background rate
;	licu = fondo
;ENDIF

end
