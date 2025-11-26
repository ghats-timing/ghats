pro gh_oplot_nupower,frequency,power,power_err
;+
; NAME: 
;      GH_OPLOT_NUPOWER
; PURPOSE: 
;      Overplots a power spectrum in nuPnu to the current window
; EXPLANATION:
;      This procedure overplots a power spectrum in nuPnu to the current window.
;
; CALLING SEQUENCE: 
;       MU_OPLOT_NUPOWER,FREQUENCY,POWER,POWER_ERR
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
;		T. Belloni  01 Dec 2010  from mu6
;-
;--------------------------------------------------------------------------
;
pwr     = power     * frequency
pwr_err = power_err * frequency
indici=where(pwr gt 0.0)
x=frequency(indici)
y=pwr(indici)

gh_oplot_power,frequency,y,pwr_err(indici)
end
