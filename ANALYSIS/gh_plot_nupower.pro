pro gh_plot_nupower,frequency,power,power_err,x1,x2,y1,y2,ps=psopt,help=help
;+
; NAME: 
;      GH_PLOT_POWER
; PURPOSE: 
;      Plots a power spectrum in nuPu in log-log to the current window
; EXPLANATION:
;      This procedure plots a power spectrum in nuPu in log-log 
;      to the current window.
;
; CALLING SEQUENCE: 
;       GH_PLOT_NUPOWER,FREQUENCY,POWER,POWER_ERR,X1,X2,Y1,Y2[,/PS][,/HELP]
; INPUTS:
;       FREQUENCY= Frequency array
;       POWER    = Power array
;       POWER_ERR= Array of errors on power
;       X1       = Optional minimum X value for plot
;       X2       = Optional maximum X value for plot (mandatory if X1 is set)
;       Y1       = Optional minimum Y value for plot (mandatory if 
;                                                     X1,X2 are set)
;       Y2       = Optional maximum Y value for plot (mandatory if Y1 is set)
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       PS       = If set, output goes to a PS file
;       HELP     = If set, print usage information and return
;
; EXAMPLE:
;       Plots a power spectrum in nuPnu with frequencies from 0.1 to 64.0 Hz:
;
;       MU> MU_NUPLOT_POWER,nu,pow,pow_e,0.1,64.0
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  21 Nov 2001  PS option
;       T. Belloni  10 Jun 2002  removing negative points from plot
;       T. Belloni  03 Apr 2009  plot_oo removed
;	    T. Belloni  05 Apr 2009  gmu version adapted
;		T. Belloni  01 Dec 2010  from mu
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_PLOT_NUPOWER'
   print,''
   print,'Plot a PDS as f*P(f) in log-log form.'
   print,''
   print,'Usage:'
   print,'  GH_PLOT_NUPOWER, frequency, power, power_err'
   print,'  GH_PLOT_NUPOWER, frequency, power, power_err, x1, x2'
   print,'  GH_PLOT_NUPOWER, frequency, power, power_err, x1, x2, y1, y2'
   print,''
   print,'Arguments:'
   print,'  frequency  Frequency array, in Hz'
   print,'  power      Power array; plotted as frequency*power'
   print,'  power_err  Error array for power; scaled by frequency'
   print,'  x1,x2      Optional frequency range to plot'
   print,'  y1,y2      Optional frequency*power range to plot'
   print,''
   print,'Keywords: /PS writes gh_nupower.ps; /HELP prints this message.'
   print,''
   return
endif
;
common sis,sistema

pwr     = power     * frequency
pwr_err = power_err * frequency
indici  = where(pwr gt 0.0)
x=frequency(indici)
y=pwr(indici)
;
; same as plot_power
;
n    = n_elements(frequency)
if(N_params() gt 3) then begin
   xmin = x1
   xmax = x2
  endif else begin
   xmin = frequency(0) * 0.8
   xmax = frequency(n-1) * 1.2
endelse
if(N_params() gt 3) then begin
   ymin = y1
   ymax = y2
  endif else begin
   ymin = min(y)   * 0.8
   ymax = max(y)   * 1.2
endelse
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_nupower.ps'
endif
plot,x,y,psym=10,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle='Power*Frequency', $
     /xlog,/ylog

if(sistema eq 'IDL') then begin
   muerrplot,x,y-pwr_err(indici),y+pwr_err(indici),width=0.0
endif else begin
	muploterr,x,y,pwr_err(indici),/xlog,/ylog,psym=3
endelse
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif

end
