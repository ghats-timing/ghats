PRO GH_RMF,F1,F2,RMFNAME,TELESCOPE,INSTRUMENT,help=help
;+
; NAME: 
;      GH_RMF
; PURPOSE: 
;      Produces a RMF file from a PDS
; EXPLANATION:
;      This procedure is used by MU_XSPEC to produce a FITS RMF file
;      in combination with MU_PHA. It should not be used interactively.
;
; CALLING SEQUENCE: 
;       MU_RMF,F1,F2,RMFNAME
; INPUTS:
;       F1       = Array with low boundaries for frequency bins
;       F2       = Array with high boundaries for frequency bins
;       RMFNAME  = String containing the rmf filename
;
; OUTPUTS:
;       The PHANAME.rmf file is produced.
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
;		T. Belloni  11 Apr 2009  keep END in header and no print,hcr
;		T. Belloni  01 Dec 2010  from mu6
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_RMF'
   print,''
   print,'Write a diagonal XSPEC/OGIP RMF/RSP file. Normally called by GH_XSPEC.'
   print,''
   print,'Usage:'
   print,'  GH_RMF, f1, f2, rmfname, telescope, instrument'
   print,''
   print,'Arguments:'
   print,'  f1, f2      Low and high boundaries for frequency bins.'
   print,'  rmfname     Output RMF/RSP filename.'
   print,'  telescope   TELESCOP header value.'
   print,'  instrument  INSTRUME header value.'
   print,''
   print,'Keywords: /HELP prints this message.'
   return
endif

ncol   = 6l
nrow   = long(n_elements(f1))
n_grp  = intarr(nrow) + 1
f_chan = intarr(nrow) + 1
n_chan = intarr(nrow) + fix (nrow)
matrix = fltarr(nrow,nrow)

for ii = 0,nrow-1 do begin
    matrix(ii,ii) = 1.0
endfor

;create primary header

fxhmake,hdr,/extend,/date

fxaddpar,hdr,'TELESCOPE',telescope
fxaddpar,hdr,'INSTRUME',instrument

fxwrite,rmfname,hdr

;create the extension header

fxbhmake,hdr,nrow,'MATRIX','name of this binary table extension'

fxaddpar,hdr,'FILTER','NONE'
fxaddpar,hdr,'CHANTYPE','PHA'
fxaddpar,hdr,'CHANTYPE','PHA'
fxaddpar,hdr,'DETCHANS',nrow
fxaddpar,hdr,'LO_THRES',0.0d0
fxaddpar,hdr,'LO_THRES',0.0d0
fxaddpar,hdr,'EFFAREA',1.0d0
fxaddpar,hdr,'RMFVERSN','1992a'

;now create the columns
fxbaddcol,col1,hdr,F1(0),'ENERG_LO',tunit='keV     '
fxbaddcol,col2,hdr,F2(0),'ENERG_HI',tunit='keV     '
fxbaddcol,col3,hdr,n_grp(0),'N_GRP ',tunit='        '
fxbaddcol,col4,hdr,f_chan(0),'F_CHAN',tunit='        '
fxbaddcol,col5,hdr,N_CHAN(0),'N_CHAN',tunit='        '
fxbaddcol,col6,hdr,MATRIX(*,0),'MATRIX',tunit='        '

nh = WHERE(STRMID(HDR,0,8) EQ 'END     ', nend)

;Open a new binary table at the end of a FITS file.

fxbcreate,unit,rmfname,hdr

;write  the data

NORM=1.

FOR I=1L,NROW DO BEGIN
  FXBWRITE,UNIT,F1(I-1),COL1,I
  FXBWRITE,UNIT,F2(I-1),COL2,I
  FXBWRITE,UNIT,N_GRP(I-1),COL3,I
  FXBWRITE,UNIT,F_CHAN(I-1),COL4,I
  FXBWRITE,UNIT,N_CHAN(I-1),COL5,I
  FXBWRITE,UNIT,MATRIX(*,I-1),COL6,I
ENDFOR

FXBFINISH,UNIT

;create the extension header

fxbhmake,hdr,nrow,'EBOUNDS','name of this binary table extension'

fxaddpar,hdr,'CHANTYPE','PHA'
fxaddpar,hdr,'DETCHANS',nrow
fxaddpar,hdr,'EFFAREA',1.0d0
fxaddpar,hdr,'FILTER','NONE'
fxaddpar,hdr,'RMFVERSN','1992a'

;now create the columns
channel = long(indgen(nrow)+1)

fxbaddcol,col1,hdr,channel(0),'CHANNEL',tunit='     '
fxbaddcol,col2,hdr,F1(0),'E_MIN',tunit='keV     '
fxbaddcol,col3,hdr,F2(0),'E_MAX ',tunit='keV     '

nh = WHERE(STRMID(HDR,0,8) EQ 'END     ', nend)

;Open a new binary table at the end of a FITS file.

fxbcreate,unit,rmfname,hdr

;write  the data

NORM=1.

FOR I=1L,NROW DO BEGIN
  FXBWRITE,UNIT,channel(I-1),COL1,I
  FXBWRITE,UNIT,F1(I-1),COL2,I
  FXBWRITE,UNIT,F2(I-1),COL3,I
ENDFOR

FXBFINISH,UNIT
END
