      pro poisson_estimate,freq_dep_poi,cnts,time_per_spectrum,      $
                               current_vle_rate,nfreq,ndet,i_vle,tdead,  $
			                   differential=diff
;
; INPUTS:
;     CNTS
;     TIME_PER_SPECTRUM  = Length of each time stretch
;     CURRENT_VLE_RATE   = VLE rate corresponding to stretch
;     NFREQ              = Number of frequencies in the power spectrum
;     NDET               = Number of working detectors
;     I_VLE              = Indec for the VLE window (addresses TAU_VLE)
;     TDEAD              = Dead time 
;
; OUTPUTS:
;     FREQ_DEP_POI
;
;
      poisson = fltarr(nfreq)
      tau_vle = [20e-6,55e-6,150e-6,550e-6]
;
      iwarnings = 0
;
; ** Initialize **
;
      if(n_elements(freq_dep_poi) eq 0) then freq_dep_poi = fltarr(nfreq)
      rmin  = 1.0e30
      rmax  = -1.e30
      nft   = nfreq *2
      tbin  = time_per_spectrum / float(nft)      
      ratio = tbin / tdead
      vle   = current_vle_rate / ndet
      r0    = cnts/ time_per_spectrum / ndet
      fvle  = tau_vle(i_vle-1)
      indice= lindgen(nfreq)+1

	  pvle         = 2.0*vle*r0*fvle^2 *                         $
	                 (sin(!pi*fvle*indice/time_per_spectrum)/    $
                     (!pi*fvle*indice/time_per_spectrum))^2

      if(ratio ge 1) then begin
         poisson      = 2*(1.0-2.0*r0*tdead*(1.0-tdead/(2.0*tbin)))  -    $
                        2*(nft-1)/nft*r0*tdead*(tdead/tbin)*cos((2*!pi*indice)/nft)
       endif else begin
	     print,' Formula (27) of Zhang et al. not yet thoroughly tested'
	     m        = int(tdead/tbin)
	     td_tb    = tdead/tbin
         pi_N     = !pi / nft
         poisson           = 2.0-2.0*r0*tbin*(1.0-(float((nft-m))/float(nft))*      $
	                    (m+1-td_tb)^2.0*cos(m*2.0*indice*pi_N)                      $
                        +((float(nft-m-1.0))/float(nft))*(m-td_tb)^2*               $
			            cos((m+1)*indice*2.0*pi_N)                                  $ 
                        +2.0*cos(float(m+1)/2.0*2.0*indice*pi_N)*                   $
			            sin(float(m)/2.0*2.0*indice*pi_N)/sin(pi_N*float(indice))   $
                        -(float(m+1)/nft)*(sin(float(2*m+1.0)/2.0*2.*indice*pi_N)   $
			            /sin(indice*pi_N))                                          $
                        +1.0/float(nft)*((sin(float(m+1)/2.*2.0*indice*pi_N)^2/     $
			            (sin(indice*pi_N)^2))))
      endelse
	; **
	; ** Note that weighted sum of Poisson levels is calculated here 
	; ** The final averaging should take place in calling routine 
	; **
	      if(keyword_set(diff)) then begin
	; If required, don't integrate and return the current poissonian level
	          freq_dep_poi = (poisson + pvle)
	        endif else begin
	          freq_dep_poi = freq_dep_poi + float(poisson + pvle)
	      endelse
;
      end
