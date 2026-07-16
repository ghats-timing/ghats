pro read_fft_line,unit,muflag,rmjd,cnts,poisson,current_vle_rate,fndet,rdata

rmjd             = 0.0d0
cnts             = 0.0
poisson          = 0.0
current_vle_rate = 0.0
fndet            = 0.0
a0               = 0.0
;
; actual reading
;
     readu,unit,rmjd,cnts,fndet,a0,poisson,current_vle_rate,  $
          std21,std22,std23,rdata

end
