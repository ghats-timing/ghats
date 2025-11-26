PRO buoni_1348,buoni,flares=flares
	; questi sono i punti buoni per eliminare intervalli senza la type-B
   o0 = indgen(92)+218                                                              ; 107 
   ;o1 = indgen(180)+1    + indgen(68)+244 + indgen(99)+380  + o0[n_elements(o0)-1]  ; 108
   o1 = [indgen(180)+1, indgen(68)+244, indgen(99)+380] + o0[n_elements(o0)-1]      ; 108
   o2 = indgen(650)+1    + o1[n_elements(o1)-1]                                     ; 109   
   o3 = indgen(33)+142   + o2[n_elements(o2)-1]                                     ; 111 
   o4 = indgen(217)+1    + o3[n_elements(o3)-1]                                     ; 112
   o5 = indgen(102)+1    + o4[n_elements(o4)-1]                                     ; 113
   o6 = indgen(201)+104  + o5[n_elements(o5)-1]                                     ; 117
   o7 = indgen(306)+1    + o6[n_elements(o7)-1]                                     ; 124
   o8 = indgen(285)+1    + o7[n_elements(o8)-1]                                     ; 125
   buoni1 = [o0,o1,o2,o3,o4,o5,o6,o7,o8]
   buoni1 = o0
   ;
   ; adesso si selezionano gli intervalli senza flares
   ;
   IF(keyword_set(flares)) THEN BEGIN
	   gh_licu,'B9_onesec.pds',t9,b9   ; read the 12-15 keV light curve
	   gheflare = where(b9 gt 2.0)      ; identify possible flares
	   delta9   = t9[1]-t9[0]           ; binsize at high res  
	   t99      = t9[gheflare]
	   t9e      = t9+delta9              ; array of end times
	   
	   gh_licu,'B2a.pds',time,foo       ; read in one of the light curves to get times
	   deltat = time[1]-time[0]        ; data bin size
       timee  = time+deltat            ; array of end times
	   
	   b1 = []
	   ; loop over 13s bins
	   
	   FOR i=0,n_elements(time)-1 DO BEGIN
		   t  = time[i]
		   te = timee[i]
		   bad1 = where((t99 GE t) AND (t99 LT te))
		   IF(bad1[0] EQ -1) THEN BEGIN
			   b1 = [b1, i+1]
		   ENDIF
	   ENDFOR
	  
	   buoni = intersection(buoni1, b1)
   ENDIF ELSE BEGIN
	   buoni = buoni1
   ENDELSE
END