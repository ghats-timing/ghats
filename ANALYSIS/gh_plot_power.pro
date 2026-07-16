pro gh_plot_power,frequency,power,power_err,x1,x2,y1,y2,ps=psopt
;+
; NAME: 
;      GH_PLOT_POWER
; PURPOSE: 
;      Plots a power spectrum in log-log to the current window
; EXPLANATION:
;      This procedure plots a power spectrum in log-log to the current window.
;
; CALLING SEQUENCE: 
;       GH_PLOT_POWER,FREQUENCY,POWER,POWER_ERR,X1,X2,Y1,Y2[,/PS]
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
;
; EXAMPLE:
;       Plots a power spectrum with frequencies from 0.1 to 64.0 Hz:
;
;       MU> GH_PLOT_POWER,nu,pow,pow_e,0.1,64.0
;
; COMMON BLOCKS: 
;      sistema    = IDL or GDL for plotting error bars 
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
;		T. Belloni  06 May 2010  from mu
;-
;--------------------------------------------------------------------------
;
common sis,sistema

indici=where(power gt 0.0)
x=frequency(indici)
y=power(indici)
n    = n_elements(frequency)
if (N_params() gt 3) then begin
   xmin = x1
   xmax = x2
  endif else begin
   xmin = frequency(0) * 0.8
   xmax = frequency(n-1) * 1.2
endelse
if (N_params() eq 7) then begin
   ymin = y1
   ymax = y2
  endif else begin
   ymin = min(y    )   * 0.8
   ymax = max(y    )   * 1.2
endelse
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_power.ps'
endif
;
plot,x,y,psym=10,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle='Power',/xlog,/ylog

if(sistema eq 'IDL') then begin
    	muerrplot,x,y-power_err(indici),y+power_err(indici),width=0.0
endif else begin
	muploterr,x,y,power_err(indici),/xlog,/ylog,psym=3
endelse
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif
end
