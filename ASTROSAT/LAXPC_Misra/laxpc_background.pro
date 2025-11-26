PRO laxpc_background,en1,en2,back=back,drm=drm

	total_background = 0    ; output value
; read in detector responses
OPENR,lun,drm,/GET_LUN
   rmf = STRARR(3)
   READF,lun,rmf
FREE_LUN,lun

; read in background files
OPENR,lun,back,/GET_LUN
   fondo = STRARR(3)
   READF,lun,fondo
FREE_LUN,lun

; loop on the units
FOR iunit=0,2 DO BEGIN
	; read in channel-energy conversion
	fxbopen,unitr,rmf[iunit],1,hea,errmsg=errmsg
	   ttype       = fxpar(hea,'TTYPE*')
	   cha         = where(ttype eq 'CHANNEL ')+1
	   e1          = where(ttype eq 'E_MIN   ')+1
	   e2          = where(ttype eq 'E_MAX   ')+1
	   fxbreadm,unitr,[cha,e1,e2],can,emin,emax
	fxbclose,unitr
	; read in background spectrum
	fxbopen,unitb,fondo[iunit],1,hea,errmsg=errmsg
	   fxbreadm,unitb,[1,2],canale,bkg
	fxbclose,unitb 
	; find channels corresponding to energy selection for current unit
	c1 = min(where(emin ge en1))
	c2 = max(where(emax le en2))
	
	fondo_unita = total(bkg(c1:c2))
	print,'Unit',(iunit+1)*10,': ',fondo_unita
	total_background = total_background + fondo_unita
ENDFOR
print,'Total background (',en1,' - ',en2,') keV: ',total_background
END
