PRO GH_PHA,POWER,POWER_ERR,QQ,PHANAME,RMFNAME,TELESCOPE,INSTRUMENT, $
           extra_keys=extra_keys,extra_values=extra_values,help=help
;+
; NAME: 
;      GH_PHA
; PURPOSE: 
;      Produces a PHA file from a PDS
; EXPLANATION:
;      This procedure is used by MU_XSPEC to produce a FITS PHA file
;      containing a PDS. It should not be used interactively.
;
; CALLING SEQUENCE: 
;       MU_PHA,POWER,POWER_ERR,PHANAME,RMFNAME
; INPUTS:
;       POWER    = Array with powers (watch out, power per bin!)
;       POWER_ERR= Array with power errors
;       PHANAME  = String containing the pha filename
;       RMFNAME  = String containing the rmf filename (for the pha header)
;
; OUTPUTS:
;       The PHANAME.pha file is produced.
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
;       FX* routines from the astron IDL library
; NOTES:
;       NONE
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation: original from T. Yaqoob
;       T. Belloni   9 Nov 2001  removal of remnant arf filename
;       T. Belloni   7 Jun 2005  useless ANCRFILE keyword added
;		T. Belloni  11 Apr 2009  keep END in header and no print,hcr
;		T. Belloni	01 Dec 2010  from mu6
;       M. Mendez   29 Sep 2023  write column with quality of channel

;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_PHA'
   print,''
   print,'Write an XSPEC/OGIP PHA file. Normally called by GH_XSPEC.'
   print,''
   print,'Usage:'
   print,'  GH_PHA, power, power_err, quality, phaname, rmfname, telescope, instrument'
   print,''
   print,'Arguments:'
   print,'  power, power_err  Values and errors to write.'
   print,'  quality           XSPEC quality flags.'
   print,'  phaname           Output PHA filename.'
   print,'  rmfname           Response filename recorded in the PHA header.'
   print,''
   print,'Keywords: EXTRA_KEYS= and EXTRA_VALUES= add FITS header keywords; /HELP prints this message.'
   return
endif

;ncol   = 6l
; MM ncol   = 3l
ncol   = 4l
; MM

nrow   = long(n_elements(POWER))
channel  = fix(findgen(nrow) + 1)
quality  = intarr(nrow)
sys_err  = fltarr(nrow)
grouping = intarr(nrow) + 1

;create primary header

fxhmake,hdr,/extend,/date
fxaddpar,hdr,'TELESCOPE',telescope
fxaddpar,hdr,'INSTRUME',instrument
fxaddpar,hdr,'CONTENT','SPECTRUM'
fxaddpar,hdr,'PHAVERSN','1992a'

; Create the file writing primary FITS header
fxwrite,phaname,hdr

;create the extension header

fxbhmake,hdr,nrow,'SPECTRUM','name of this binary table extension',/initialize
fxaddpar,hdr,'TELESCOPE',telescope
fxaddpar,hdr,'INSTRUME',instrument
fxaddpar,hdr,'FILTER','NONE'
fxaddpar,hdr,'EXPOSURE',double(1.0)
fxaddpar,hdr,'AREASCAL',double(1.0)
fxaddpar,hdr,'BACKSCAL',double(1.0)
fxaddpar,hdr,'CORRSCAL',double(1.0)
fxaddpar,hdr,'BACKFILE','NONE'
fxaddpar,hdr,'CORRFILE','NONE'
fxaddpar,hdr,'ANCRFILE','NONE'
fxaddpar,hdr,'RESPFILE',RMFNAME
fxaddpar,hdr,'POISSERR','F'
fxaddpar,hdr,'CHANTYPE','PHA'
fxaddpar,hdr,'DETCHANS',fix(nrow)
fxaddpar,hdr,'SYS_ERR',0
fxaddpar,hdr,'QUALITY',0
fxaddpar,hdr,'GROUPING',0

fxaddpar,hdr,'HDUCLASS','OGIP'
fxaddpar,hdr,'HDUCLAS1','SPECTRUM'
fxaddpar,hdr,'HDUVERS','1.1.0'

if n_elements(extra_keys) gt 0 then begin
   if n_elements(extra_values) eq n_elements(extra_keys) then begin
      for ih=0L,n_elements(extra_keys)-1L do begin
         fxaddpar,hdr,strtrim(string(extra_keys[ih]),2),strtrim(string(extra_values[ih]),2)
      endfor
   endif
endif

;now create the columns
fxbaddcol,col1,hdr,CHANNEL(0),'CHANNEL',tunit='       '
fxbaddcol,col2,hdr,POWER(0),'COUNTS',tunit='counts   '
fxbaddcol,col3,hdr,POWER_ERR(0),'STAT_ERR ',tunit='counts   '
; MM
fxbaddcol,col4,hdr,QQ(0),'QUALITY',tunit='         '
; MM

nh = WHERE(STRMID(HDR,0,8) EQ 'END     ', nend)
;hdr=hdr(0:nh(0)-1)
;print,'hdr: ',(size(hdr))(1)
;s=strarr(3)
;s(0)='TELESCOP= ''XTE     ''           / mission/satellite name'
;s(1)='INSTRUME= ''PCA     ''           / instrument/detector name'
;s(2)='END                                                         '

;print,'hdr: ',(size(hdr))(1)

;Open a new binary table at the end of a FITS file.
fxbcreate,unit,phaname,hdr

;write  the data

NORM=1.

FOR I=1L,NROW DO BEGIN
  FXBWRITE,UNIT,CHANNEL(I-1),COL1,I
  FXBWRITE,UNIT,POWER(I-1),COL2,I
  FXBWRITE,UNIT,POWER_ERR(I-1),COL3,I
; MM
  FXBWRITE,UNIT,QQ(I-1),COL4,I
; MM
ENDFOR
FXBFINISH,UNIT
END
