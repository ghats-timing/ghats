PRO concatenate_pds,pdslist,outputfilename

nlines=file_lines(pdslist)
OPENR,lun,pdslist,/GET_LUN
   filenames = STRARR(1, nlines)
   READF,lun,filenames
FREE_LUN,lun
; open output file
openw,output_unit,outputfilename,/GET_LUN

; open first file and get header from that one
openu,unit,filenames[0],/GET_LUN

gh_version_string = '                '
nft               = 0L
T                 = 0.0d0
observatory       = '                '
target            = '                '
instrument        = '                '
rmjd0             = 0.0d0
ntrafos           = 0l
e                 = intarr(2)
proliferation     = 0
baryflag          = 0
n_spectral_bins   = 0
background_flag   = 0
dummy               = bytarr(100)
;
; Read in header
;
readu,unit,gh_version_string,             $
                   observatory,           $
                   instrument,            $
                   target,                $
                   rmjd0,                 $
                   nft,                   $
                   T,                     $
                   ntrafos,               $
                   e,                     $
                   proliferation,         $
                   baryflag,              $
                   n_spectral_bins,       $
                   background_flag,       $
                   dummy
writeu,output_unit,gh_version_string,     $
                   observatory,           $
                   instrument,            $
                   target,                $
                   rmjd0,                 $
                   nft,                   $
                   T,                     $
                   ntrafos,               $
                   e,                     $
                   proliferation,         $
                   baryflag,              $
                   n_spectral_bins,       $
                   background_flag,       $
                   dummy
;
; Loop over the trafos
;
nft       = nft/2

;
pwr       = fltarr(nft)*0.0
rmjd             = 0.0d0
cnts             = 0.0
poisson          = 0.0
current_vle_rate = 0.0
fndet            = 0.0
a0               = 0.0
nbins            = 0
FOR itrafos=0l,ntrafos-1l DO BEGIN

   readu,unit,rmjd,cnts,fndet,a0,poisson,current_vle_rate,   $
	    std21,std22,std23,pwr
   writeu,output_unit,rmjd,cnts,fndet,a0,poisson,current_vle_rate,   $
    	std21,std22,std23,pwr
   
ENDFOR

FREE_LUN,unit
ntotal = ntrafos

;
; Now reading the other files. All parameters are supposed to be the same (and unchecked)
; but what needs to be changed is the rmjd, which must be related to rmjd0 of the
; first file. Since rmjd is already MJD, all I need to do is read and write the files
; (I think)
FOR i=1,nlines-1 DO BEGIN  ; loop over the remaining files
	print,i
	openu,unit,filenames[i],/get_lun
	; read in the header, but ignore it
	readu,unit,gh_version_string,             $
	                   observatory,           $
	                   instrument,            $
	                   target,                $
	                   rmjd0,                 $
	                   nft,                   $
	                   T,                     $
	                   ntrafos,               $
	                   e,                     $
	                   proliferation,         $
	                   baryflag,              $
	                   n_spectral_bins,       $
	                   background_flag,       $
	                   dummy
	
	FOR itrafos=0l,ntrafos-1l DO BEGIN

	   readu,unit,rmjd,cnts,fndet,a0,poisson,current_vle_rate,   $
		    std21,std22,std23,pwr
	   writeu,output_unit,rmjd,cnts,fndet,a0,poisson,current_vle_rate,   $
	    	std21,std22,std23,pwr
   
	ENDFOR
    free_lun,unit
	ntotal = ntotal + ntrafos
ENDFOR
FREE_LUN,output_unit
; reopen output file to reset number of trafos
openu,uu,outputfilename,/get_lun
   point_lun,uu,84
   writeu,uu,ntotal
free_lun,uu
END