pro gh_oplot_power,frequency,power,power_err
;+
; NAME: 
;      GH_OPLOT_POWER
; PURPOSE: 
;      Overplots a power spectrum to the current window
; EXPLANATION:
;      This procedure overplots a power spectrum to the current window.
;
; CALLING SEQUENCE: 
;       GH_OPLOT_POWER,FREQUENCY,POWER,POWER_ERR
; INPUTS:
;       FREQUENCY= Frequency array
;       POWER    = Power array
;       POWER_ERR= Array of errors on power
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       NONE
;
; EXAMPLE:
;       None
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  10 Jun 2002  removing negative points from plot
;	    T. Belloni  05 Apr 2009  gmu version adapted
;	    T. Belloni  05 May 2010  from mu
;-
;--------------------------------------------------------------------------
;
common sis,sistema

indici=where(power gt 0.0)
x=frequency(indici)
y=power(indici)
oplot,x,y,psym=10

if(sistema eq 'GDL') then begin
    muploterr,x,y,power_err(indici),/xlog,/ylog,psym=3
endif else begin
    muerrplot,x,y-power_err(indici),y+power_err(indici),width=0.0
endelse

end
