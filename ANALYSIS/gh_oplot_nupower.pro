pro gh_oplot_nupower,frequency,power,power_err,help=help
;+
; NAME: 
;      GH_OPLOT_NUPOWER
; PURPOSE: 
;      Overplots a power spectrum in nuPnu to the current window
; EXPLANATION:
;      This procedure overplots a power spectrum in nuPnu to the current window.
;
; CALLING SEQUENCE: 
;       GH_OPLOT_NUPOWER,FREQUENCY,POWER,POWER_ERR[,/HELP]
; INPUTS:
;       FREQUENCY= Frequency array
;       POWER    = Power array
;       POWER_ERR= Array of errors on power
;
; OUTPUTS:
;       HELP     = If set, print usage information and return
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
if(keyword_set(help)) then begin
   print,''
   print,'GH_OPLOT_NUPOWER'
   print,''
   print,'Overplot a PDS as f*P(f) on the current plot.'
   print,''
   print,'Usage:'
   print,'  GH_OPLOT_NUPOWER, frequency, power, power_err'
   print,'  GH_OPLOT_NUPOWER,/HELP'
   print,''
   return
endif
;
pwr     = power     * frequency
pwr_err = power_err * frequency
indici=where(pwr gt 0.0)
x=frequency(indici)
y=pwr(indici)

gh_oplot_power,frequency,y,pwr_err(indici)
end
