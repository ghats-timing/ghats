pro mu_read_laxpc_event_new,unit,laxpc_unit,laxpc_layer,tag,rdata,np,energie,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          BANDE,nffts
;
;  Read in LAXPC photon files
;
;-------------------------------------------------------------------------
; Parameters
;
; unit                 I: unit for file reading (already opened)
; tag                  I: file tag identifier
; rdata                O: array with extracted light curve
; np                   I: number of points to accumulate
; energie               I: array of start-end energies
; selected             I: still a mystery to me
; jmeta                I: index of current metafile
; tres                 I: time resolution of current file
; tres_fft             I: time resolution needed for light curve
; tstart_fft           I: fft start time
; tend_fft             I: fft end time
; irowl                I: current row
;-------------------------------------------------------------------------
;common dati,lay,lax,tempi1,energy1,tempi2,energy2,tempi3,energy3,canale_corrente1,canale_corrente2,canale_corrente3,std21,std22,std23 ; common block for keeping the data
common dati,lay,lax,tempi1,energy1,canale_corrente1,std21,std22,std23 ; common block for keeping the data
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)

;must find the columns (time and pi) and read in the whole file
;it should be straightforward
header      = fxbheader(unit)
ttype       = fxpar(header,'TTYPE*')
time_unit   = where(ttype eq 'TIME    ')+1
picol       = where(ttype eq 'Energy  ')+1
laxpc       = where(ttype eq 'LAXPC_No.')+1
layer       = where(ttype eq 'Layer   ')+1

;
; MU6 addition: at first call read in full file to memory using fxbreadm
IF(tag ne canale_corrente1) THEN BEGIN
	IF(baryflag eq 1) THEN BEGIN
		time_unit = fxbcolnum(unit,'barytime')
	ENDIF
	fxbreadm,unit,[time_unit,picol,laxpc,layer],tempi1,energy1,lax,lay  ; energy here is real energy in keV
	nrows = n_elements(tempi1)
	canale_corrente1 = tag
	; Select only unit and layer as desired, throw away other points
	bu1 = (lax eq (1*laxpc_unit[0])) OR (lax eq (2*laxpc_unit[1])) OR (lax eq (3*laxpc_unit[2]))
	bu2 = (lay eq (1*laxpc_layer[0])) OR (lay eq (2*laxpc_layer[1])) OR (lay eq (3*laxpc_layer[2])) OR $$
		  (lay eq (4*laxpc_layer[3])) OR (lay eq (5*laxpc_layer[4]))
	bu  = (bu1 AND bu2)

	tempi1  = tempi1[where(bu eq 1)]
	energy1 = energy1[where(bu eq 1)]
ENDIF
;
; Now we have tstart_fft and tend_fft, we need to select the photons within the time limit and energy;
; BUONI contains the indices of photons to use to build light curve

tt = tempi1 - tstart_fft
buoni = WHERE( $$
	(tt GE 0) AND (tt LT (tend_fft-tstart_fft)) $$
	 AND (energy1 GE energie[0]) AND (energy1 LT energie[1]) $$
	)
buonie = WHERE((tempi1 GE tstart_fft) AND (tempi1 LT tend_fft))
; Binning of photons into RDATA. T is good times with a TSTART_FFT time offset for small numbers (start from 0); E is good energies,
; used for main rates accumulation

;t = tempi1[buoni]-tstart_fft
t  = tt[buoni]
e = energy1[buonie]
rdata = histogram(t,BINSIZE=tres_fft,NBINS=n_elements(rdata))
std21       = total((e ge BANDE[0])*(e le BANDE[1]))
std22       = total((e ge BANDE[2])*(e le BANDE[3]))
std23       = total((e ge BANDE[4])*(e le BANDE[5]))

END