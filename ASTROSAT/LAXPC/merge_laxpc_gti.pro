pro merge_laxpc_gti,namelist,output_gti_file UNFINISHED!

	; namelist must be a file containing the paths to orbit folders, the ones that contain
	; lxp*/modeEA
	;
	nlines=file_lines(namelist)
	OPENR, lun, namelist, /GET_LUN
	  osservazioni = STRARR(1, nlines)
	  READF, lun, osservazioni
	CLOSE, lun
	FREE_LUN, lun
	nlines=file_lines(namelist)
	
	FOR i=0,nlines-1 DO BEGIN
		dirname = osservazioni[i]
		lxp1    = FILE_SEARCH(dirname+'/lxp1/modeEA/*.gti')
		lxp2    = FILE_SEARCH(dirname+'/lxp2/modeEA/*.gti')
		lxp3    = FILE_SEARCH(dirname+'/lxp3/modeEA/*.gti')
		;lxp1    = lxp1[0]
		;lxp2    = lxp2[0]
		;lxp3    = lxp3[0]
		
		fxbopen,u1,lxp1,1,h1
		   fxbreadm,u1,[1,2],gstart1,gstop1
		fxbclose,u1
		fxbopen,u2,lxp2,1,h2
		   fxbreadm,u2,[1,2],gstart2,gstop2
		fxbclose,u2
		fxbopen,u3,lxp3,1,h3
		   fxbreadm,u3,[1,2],gstart3,gstop3
		fxbclose,u3
		ooo
	ENDFOR

END