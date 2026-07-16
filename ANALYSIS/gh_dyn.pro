pro gh_dyn,filename,      $
             tb,licub,nu,dynimage3,index=index,frebin=reb,trebin=treb
;+
; NAME:
;      GH_DYN
; PURPOSE:
;      Reads is a time-frequency image  from a PDS file, 
;      selecting on a range of PDS
; EXPLANATION:
;      This procedure extracts a time-frequency image  from a PDS file.
;      A selection is possible only on a continuous range of PDS.
;
; CALLING SEQUENCE:
;       GH_DYN,FILENAME,TIMES,LICU,FREQUENCY,DYNIMAGE,STARTING,ENDING
; INPUTS:
;       FILENAME = name of the input PDS file
;       STARTING = index of first PDS to be selected (not mandatory)
;       ENDING   = index of last  PDS to be selected (not mandatory)
;       REB      = rebin factor in frequency (can be negative for LOG)
;       TREB     = rebin factor in time (cannot be negative for LOG)
;
; OUTPUTS:
;       TIME     = Array for times (time axis of dynimage)
;       FREQUENCY= Array for frequencies (frequency axis of dynimage)
;       DYNIMAGE = Time-frequency power array (dim1 is time, dim2 is frequency)
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
;       GHATS_OPENPDS:    Opens a PDS file
;       GHATS_GETHEADER:  Reads in the header
;       READ_PDS_LINE:     Reads in next line from PDS file
;       GHREBIN:           Rebinning of power spectra
; NOTES:
;       None
; MODIFICATION HISTORY:
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  12 Nov 2001  modular version + rearranged order of pars +
;                                keyword+fast exit
;       T. Belloni  12 Mag 2002  adapted to MUFFT format + rebin factor
;       T. Belloni  19 Jul 2002  added time rebinning option
;       T. Belloni  11 Dec 2003  Mu3 version (colors)
;       P. Casella   7 Sep 2005  fixed bug in time rebinning
;		T. Belloni  06 May 2010  from mu
;		T. Belloni   9 Feb 2012  changed default frequency rebinning from -100 to 1
;		T. Belloni  19 Feb 2012  no floor rounding for times array
;-
;--------------------------------------------------------------------------
;

if(keyword_set(index)) then begin
    starting = index[0]
    ending  = index[1]
endif else begin
    starting = 0L
    ending   = 100000L
endelse

if(keyword_set(reb)) then begin

endif else begin
    reb=1
endelse

if(keyword_set(treb)) then begin

endif else begin
    treb=1
endelse

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
;   nft is number of TIME points. Must be divided by 2 to get frequencies
;
nft       = nft/2
;
pwr       = fltarr(nft)*0.0
;
times     = fltarr(ntrafos)*0.0
licu      = fltarr(ntrafos)*0.0
dynimage  = fltarr(ntrafos,nft)*0.0
;
; Loop over the trafos
;

n_selected       = 0l
itrafos          = 0l

for itrafos=0l,ntrafos-1 do begin
;
   read_pds_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,a0,pwr
;
   times(itrafos)  = ((rmjd-rmjd0)*86400.0)
   licu(itrafos)   = cnts*T
;
;  Fill in dynamic power spectrum image
;
   dynimage(itrafos,*) = pwr
;
   if(itrafos eq (ending-1)) then goto,finito
endfor
finito:

close,unit
;
;  Now fill in frequency array
;
   frequency = (findgen(nft)+1.0) * T
;
;  Logarithmic rebinning in frequency
;
if(reb ne 1) then begin
   ghrebin,frequency,frequency,dynimage(0,*),dynimage(0,*)*0.0,reb, $
           nu,nue,newdin,newdine,n2
   dynimage2=fltarr(ntrafos,n2)
   for i=0,ntrafos-1 do begin
      ghrebin,frequency,frequency,dynimage(i,*),dynimage(i,*)*0.0,reb, $
              nu,nue,newdyn,newdyne,n2
      dynimage2(i,*) = newdyn
   endfor
  endif else begin
      dynimage2 = dynimage
      nu        = frequency
      n2        = nft
  endelse
if(treb ne 1) then begin
;
;  First time use licu instead of times to rebin rate curve as well
;
   ghrebin,licu,licu,dynimage2(*,0),dynimage2(*,0)*0.0,treb,  $
           licub,licube,newd,newde,n3
   dynimage3=fltarr(n3,n2)
   for i=0,n2-1 do begin
      ghrebin,times,times,dynimage2(*,i),dynimage2(*,i)*0.0,treb,   $
              tb,tbe,newd,newde,n3
      dynimage3(*,i) = newd
   endfor
  endif else begin
      dynimage3 = dynimage2
      tb        = times
      licub     = licu
endelse
;
   n_selected = ntrafos
   ps = ' spectra '
   if(n_selected eq 1) then ps = ' spectrum '
   print,'  ',strtrim(string(n_selected),1),' power'+ps+'read in'
end
