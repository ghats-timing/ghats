PRO gh_split_fft, file_in, gap_threshold=gap_threshold
  ; file_in: Input .fft filename
  ; gap_threshold: Time in seconds (Default: 1200s)

  IF (N_ELEMENTS(gap_threshold) EQ 0) THEN gap_threshold = 1200.0 
  gap_threshold = gap_threshold/86400.0
  
  ghats_openfft,file_in,unit

  gh_version_string = '                '
  observatory       = '                '
  instrument        = '                '
  target            = '                '
  rmjd0             = 0.0D0
  nft               = 0l
  ntrafos           = 0l
  dummy             = bytarr(100)
  canali            = intarr(2)

  ghats_getheader,unit,gh_version_string,observatory,instrument,target,rmjd0, $
				nft,T,ntrafos,canali,proliferation,baryflag,n_spectral_bins, $
				background_flag,dummy

  nfreqs = nft
  rdata = fltarr(nfreqs)*0.0

  file_base = STRMID(file_in, 0, STRPOS(file_in, '.', /REVERSE_SEARCH))
  
  PRINT, '====================================================================================='
  PRINT, 'Splitting FFT: ', file_in
  PRINT, '  for gaps > ' , gap_threshold*86.400, ' ks'
  PRINT, 'Total segments: ', ntrafos
  PRINT, '====================================================================================='

  orbit_count = 0
  rmjd_prev = -1.0d0
  nffts_in_orbit = 0L
  oufilename = ''
  ounit=98 
  
  FOR i = 0L, ntrafos - 1 DO BEGIN
;      read_fft_line, unit, muflag, rmjd_curr, cnts, poi, vle_rate, ndet, a0, rdata
      
      rmjd_curr        = 0.0d0
      cnts             = 0.0
      poisson          = 0.0
      current_vle_rate = 0.0
      fndet            = 0.0
      a0               = 0.0
      
      readu,unit,rmjd_curr,cnts,fndet,a0,poisson,current_vle_rate,   $
           std21,std22,std23,rdata
           
;      print,rmjd_curr,cnts,fndet,a0,poisson,current_vle_rate,std21,std22,std23
      
      ; Orbit Jump Detection (RMJD is in days)
      IF (rmjd_prev LT 0) OR (rmjd_curr - rmjd_prev GT gap_threshold) THEN BEGIN
          free_lun, ounit
                    
          ; Close and patch previous file header
          IF (nffts_in_orbit GT 0) THEN BEGIN
             PRINT, 'Closing orbit ', orbit_count, ' with ', nffts_in_orbit, $
                    ' segments.    ( Gap = ', (rmjd_curr - rmjd_prev)*86.400 , ' ks ).'
             OPENU, uu, oufilename, /GET_LUN
             POINT_LUN, uu, 84
             WRITEU, uu, LONG(nffts_in_orbit)
             FREE_LUN, uu
          ENDIF
          
          ; Create new orbit file
          orbit_count = orbit_count + 1
          nffts_in_orbit = 0L
          oufilename = file_base + '_orb' + STRING(orbit_count, FORMAT='(I2.2)') + '.fft'
          PRINT, '-> Opening: ', oufilename
          new_output_flag = 1 ; Force new file creation/header
      ENDIF ELSE BEGIN
          new_output_flag = 0 ; Append to existing
      ENDELSE

      if(new_output_flag eq 1) then begin
        openw,ounit,oufilename,error=errore; ,/compress  ; compressed format
        if(errore ne 0) then begin
          massage,'Cannot open output file! ERROR #'+string(errore)
          retall
        endif
        new_output_flag = 0      ; do not try this anymore

        writeu,ounit,gh_version_string,observatory,instrument,target,  $
                     rmjd0,nft,T,ntrafos,canali,proliferation,baryflag, $
                     n_spectral_bins,background_flag,dummy
      endif

      writeu,ounit,rmjd_curr,cnts,fndet,a0,poisson,current_vle_rate, $
            std21,std22,std23,rdata

      rmjd_prev = rmjd_curr
      nffts_in_orbit = nffts_in_orbit + 1
  ENDFOR
  FREE_LUN, ounit
  
  IF (nffts_in_orbit GT 0) THEN BEGIN
     PRINT, 'Closing orbit ', orbit_count, ' with ', nffts_in_orbit, ' segments.    ( END ).'
     OPENU, uu, oufilename, /GET_LUN
     POINT_LUN, uu, 84
     WRITEU, uu, LONG(nffts_in_orbit)
     FREE_LUN, uu
  ENDIF
  PRINT, '====================================================================================='
END
