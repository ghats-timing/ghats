pro mu_intersect_gtis,gti1,gti2,u_gti1,u_gti2,nvaltimes,MAX_N_VALTIMES,nfiles,nmetafiles
	; gti1,gti2 are the input data GTIs
	; u_gti1,u_gti2 are the user's GTIs
	; nvaltimes is an array with the number of valtimes per file and metafile
	; MAX_N_VALTIMES is the maximum number of GTIs accepted (from mufft)
	; The output  intersected GTIs will be written over gti1 and gti2
	n_gti1 = dblarr(MAX_N_VALTIMES,nfiles,nmetafiles)
	n_gti2 = dblarr(MAX_N_VALTIMES,nfiles,nmetafiles)
	ngood  = intarr(nfiles,nmetafiles)*0
	nd     = size(gti1)
	
	for imetafiles=0,nmetafiles-1 do begin            ; loop over metafiles
		for ifiles=0,nfiles-1 do begin                ; loop over files
			
			;for idata=0,n_elements(gti1)-1 do begin   ; loop over the data GTIs OLD incorrect
			for idata=0,nd[1]-1 do begin   ; loop over the data GTIs
               for iuser=0,n_elements(u_gti1)-1 do begin ; loop over the user GTIs
			    	t1 = max([gti1(idata,ifiles,imetafiles),u_gti1(iuser)])
			    	t2 = min([gti2(idata,ifiles,imetafiles),u_gti2(iuser)])
	                   if(t1 lt t2) then begin  ; we found a good one
		                   n_gti1(ngood(ifiles,imetafiles),ifiles,imetafiles) = t1
		                   n_gti2(ngood(ifiles,imetafiles),ifiles,imetafiles) = t2
		                   ngood(ifiles,imetafiles) = ngood(ifiles,imetafiles) + 1
		               endif
		       endfor
		    endfor
		endfor
	endfor

;
;  Now trim arrays
;
    mass       = max(ngood)
    gti1       = n_gti1(0:mass-1,*,*)
    gti2       = n_gti2(0:mass-1,*,*)
	nvaltimes = ngood
		
		
	
;	for i=0,n_elements(gti1)-1 do begin   ; loop over the data GTIs
;		for j=0,n_elements(u_gti1)-1 do begin ; loop over the user GTIs
;			t1 = max([gti1(i),u_gti1(j)])
;			t2 = min([gti2(i),u_gti2(j)])
;			if(t1 lt t2) then begin  ; we found a good one
;				n_gti1(ngood) = t1
;				n_gti2(ngood) = t2
;				ngood = ngood + 1
;			endif
;		endfor
;	endfor
;	ngood = ngood -1
;	; trim output array
;	gti1 = n_gti1(0:ngood)
;	gti2 = n_gti2(0:ngood)
end