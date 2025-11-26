pro mu_zhang,cnts,time_per_spectrum,number_of_detectors, $
                 nft,tot_vle,i_vle,zpoi
;
; Computes poissonian levels
;
;--------------------------------------------------------------------------
; Parameters
;

tdead     = 10.5e-6  
tau_vle   = [20e-6,55e-6,150e-6,550e-6]
iwarnings = 0
;
; variables
;
poisson = fltarr(nft/2l+1l)    ; to avoid starting with 0
pvle    = fltarr(nft/2l+1l)    ; ditto
timedel = 0.0d0

fvle    = tau_vle(i_vle-1)

if (number_of_detectors ne 0) then begin
   vle = tot_vle /number_of_detectors
  endif else begin
   vle = 0.0
endelse
tbin = time_per_spectrum / float(nft)
ratio = tbin / tdead

if (number_of_detectors ne 0) then begin
   r0 = (cnts / float(nft) / tbin) / number_of_detectors 
  endif else begin
   r0 = 0.
endelse

indice = lindgen(nft/2l)+1

pvle                  = 2*vle*r0*fvle^2 *                               $
                      (sin(!pi*fvle*indice/time_per_spectrum)/          $
                      (!pi*fvle*indice/time_per_spectrum))^2

if(ratio ge 1) then begin
	;     *** FORMULA (24) of Zhang et al. 1995 ApJ 449, 930
	poisson             = 2*(1 -2*r0*tdead*(1-tdead/(2*tbin))) -         $
                          2*(nft-1)/nft*r0*tdead*(tdead/tbin)*           $
                          cos((2*!pi*indice)/nft) +pvle
endif else begin
	print,' Formula (27) of Zhang et al. not yet tested'
	;           *** FORMULA (27) of Zhang et. al. 1995 ApJ 449, 930
	m            = fix(tdead/tbin)
    td_tb        = tdead/tbin
    pi_N         = !pi / nft
    poisson           = 2 -2*r0*tbin*(1-(float((nft-m))/                      $
                        float(nft))*(m+1-td_tb)^2*cos(m*2.*indice*pi_N)       $
                        +((float(nft -m -1))/float(nft))*                     $
                        (m- td_tb)^2*cos((m-1)*indice*2.*pi_N)                $ 
                        +2*cos(float(m+1)/2.*2.*indice*pi_N)*sin(float(m)/    $
                        2.*2.*indice*pi_N)/sin(pi_N*float(indice))            $
                        -(float(m+1)/nft)*(sin(float(2*m+1)/2.                $
                        *2.*indice*pi_N)/sin(indice*pi_N))                    $
                        +1./float(nft)*((sin(float(m+1)/2.*2.*indice*pi_N)    $
                        ^2/(sin(indice*pi_N)^2)))) +pvle
endelse
zpoi = float(min(poisson))
xte_max = float(max(poisson))

end
      
