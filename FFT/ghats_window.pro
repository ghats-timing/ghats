function ghats_window,length,tipo,par=par,winn=winn
	
x = findgen(length)
t = x/length-0.5

SWITCH tipo OF
   'Boxcar':   BEGIN
	             finestra = x*0.0+1.0
	             winn=0
	             BREAK
	           END

   'Bartlett': BEGIN
	             finestra = 1.0-abs((x-length*0.5)/(length*0.5))
	             winn=1
	             BREAK
	           END
   'Hann':     BEGIN
		          finestra = 0.5*(1.0-cos(2.0*!PI*x/length))
		          winn=2
		          BREAK
			   END
   'Welch':    BEGIN
				  finestra = 1.0-((x-length*0.5)/(length*0.5))^2.0
                  winn=3
				  BREAK
			   END
   'Cosine':  BEGIN
	              finestra = cos(!PI*t)
	              winn=4
	              BREAK
	           END	
   'Hanning': BEGIN
	              finestra = (cos(!PI*t))^2.0
	              winn=5
	              BREAK
	          END
   'Hamming': BEGIN
	              IF(keyword_set(par)) THEN BEGIN
	                 finestra = par+(1.0-par)*cos(!PI*t)^2.0
	                 winn=6
	               ENDIF ELSE BEGIN
	                 PRINT,'Parameter needed for Hamming window!'
	                 finestra=-1
	              ENDELSE
	              BREAK
	          END
	'Triplet': BEGIN
		          IF(keyword_set(par)) THEN BEGIN
	                 finestra = exp(-par/abs(t))*(cos(!PI*t))^2.0
	                 winn=7
	               ENDIF ELSE BEGIN
	                 PRINT,'Parameter needed for Triplet window!'
	                 finestra=-1
	              ENDELSE
	              BREAK
	           END
	 'Gauss':  BEGIN
		          IF(keyword_set(par)) THEN BEGIN
	                 finestra = exp(-0.5*t^2.0/par^2)
	                 winn=8
	               ENDIF ELSE BEGIN
	                 PRINT,'Parameter needed for Gauss window!'
	                 finestra=-1
	              ENDELSE
	              BREAK
	           END
	  'Kaiser': BEGIN
		        	IF(keyword_set(par)) THEN BEGIN
			        xx = par*sqrt(1.0-(2.0*t)^2.0)
                    finestra = beseli(xx,0)/beseli(par,0)
                    winn=9
                  ENDIF ELSE BEGIN
                    PRINT,'Parameter needed for Kaiser-Bessel window!'
                    finestra=-1
                  ENDELSE
                BREAK
                END
  ELSE: BEGIN 
           PRINT,'Unknown window type!'
           finestra=-1
        END   
ENDSWITCH

return,finestra
END