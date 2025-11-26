pro mu_write_single_fft,rdata,nft,poi,cnts,T,rmjd,filename,   $
               ndet,               $
	       current_vle_rate,new_output_flag,output_unit,   $
	       ntotal_ffts,canali,osserv,sorgente,strumento,   $
		std21,std22,std23, $
		BANDE,turbo,gti_flag,proliferation,i_vle     ; MU6
;
; Procedure to write a RXTE FFT to disk
;
;---------------------------------------------------------------------
; Parameters written
; HEADER
;   mu_version_string                  S4   string for chek
;   nft                                L   number of frequencies per PDS
;   1.0/T                              D   inverse of licu length
;   osserv                             S16 observatory name
;   sorgente                           S16 source name
;   strumento                          S16 instrument name
;   rmjd                               R   starting time in RMJD 
;   ntotal_ffts                        L   number of light curves
;   canali                             2xI array with start end channels
;   dummy                             B100 spare space
;- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
; DATA
;   rmjd                               D   start time (RMJD)
;   cnts                               R   source rate
;   poi                                R   estimated poissonian level
;   current_vle_rate                   R   estimated VLE rate
;   float(ndet)                        R   # of detectors (real)
;   rdata                              R   array with ntotal_ffts FFTs
;---------------------------------------------------------------------
common sis, sistema   ; common block with system variable
common barycentered,baryflag
common versione,version_id
common finestre,finestra,winn

a0 = 0.0   ; for FFT it does not make sense, but must be added
;
; Open the file if necessary
;
if(new_output_flag eq 1) then begin
   openw,output_unit,filename,error=errore ; ,/compress  ; compressed format
      
   if(errore ne 0) then begin
      massage,'Cannot open output file! ERROR #'+string(errore)
      retall
   endif
   new_output_flag = 0      ; do not try this anymore
;
;  Prepare additional MU6 information to write to output header
    n_spectral_bins = 3    ; number of spectral bins
    background_flag = 0    ; this version does not include background info
    tu     = 0B
    if(turbo eq 1) then tu = 1B
    gti    = byte(gti_flag*2)
    dummy  = bytarr(100)*0B
    dummy(0)= tu+gti             ; this byte codes the turbo flag and i_vle
    dummy(1)= winn               ; window function

      writeu,output_unit,version_id       ,     $
                         osserv,                $
                         strumento,             $
                         sorgente,              $
                         rmjd,                  $
                         nft,                   $
                         1.0/T,                 $
                         ntotal_ffts,           $
                         canali,                $
                         proliferation,         $
                         baryflag,              $
                         n_spectral_bins,       $
                         background_flag,       $
                         dummy
endif
;
; Now write the data to the file
;
writeu,output_unit,rmjd,cnts,float(ndet),a0,poi,float(current_vle_rate), $
            std21,std22,std23,rdata

end
