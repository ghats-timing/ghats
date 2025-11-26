pro gh_plot_lag,frequency,lag,lag_err,x1,x2,y1,y2,ps=psopt,time=tim,lin=lin,hyp=hyp
;+
; NAME: 
;      GH_PLOT_LAG
; PURPOSE: 
;      Plots a phase/time lag spectrum in log-lin to the current window
; EXPLANATION:
;      This procedure plots a phase/time lag spectrum in log-lin to 
;      the current window. By default it considers the input to be
;      a phase lag spectrum. Time lags can be done with the switch /time
;
; CALLING SEQUENCE: 
;       GH_PLOT_LAG,FREQUENCY,LAG,LAG_ERR,X1,X2,Y1,Y2[/TIME][,/PS]
; INPUTS:
;       FREQUENCY= Frequency array
;       LAG      = Phase/time lag array
;       LAG_ERR  = Array of errors on phase/time lag
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
;       TIME     = If set, time lags, if not, phase lags
;       PS       = If set, output goes to a PS file
;		LIN      = If set, the plot will be linear in Y
;		HYP      = If set, the hyperbolic sine of the lags is plotted, allowing 
;					a sort of log plot also for negative lags
;
; EXAMPLE:
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  01 Oct 2002  from mu_plot_power
;       T. Belloni  03 Apr 2009  plot_oi removed
;	    T. Belloni  05 Apr 2009  gmu version adapted
;		T. Belloni  07 Jun 2009  fixed plot symbol in GDL
;		T. Belloni  10 Aug 2009  added 0-lag line
;		T. Belloni  09 Dec 2011  /lin option added
;		T. Belloni  26 Jan 2012  from mu_plot_cross
;		T. Belloni  01 Mar 2012  changed default PS name
;		T. Belloni  16 Mar 2012  added hyp keyword
;-
;--------------------------------------------------------------------------
common sis,sistema
;
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
   ymin = min(lag    )   * 0.8
   ymax = max(lag    )   * 1.2
endelse
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_lag.ps'
endif
;
ystr = 'Phase lag (rad)'
if(keyword_set(tim)) then begin
	ystr='Time lag (s)'
endif

if(keyword_set(hyp)) then begin
	lag = sinh(lag)
endif
;plot,frequency,lag,psym=8,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
;     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle=ystr,/xlog

if(sistema eq 'IDL') then begin
	a=findgen(16)*(!PI*2/16.0)
	usersym,cos(a),sin(a),/fill
	if(keyword_set(lin)) then begin
	    plot,frequency,lag,psym=8,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
            yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle=ystr,/xlog
        muerrplot,frequency,lag-lag_err,lag+lag_err,width=0.0
    endif else begin
	    plot,frequency,lag,psym=8,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
            yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle=ystr,/xlog,/ylog
        muerrplot,frequency,lag-lag_err,lag+lag_err,width=0.0
    endelse
endif else begin
	if(keyword_set(lin)) then begin
    	plot,frequency,lag,psym=4,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
	        yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle=ystr,/xlog
    	muploterr,frequency,lag,lag_err,/xlog,psym=3
    endif else begin
		plot,frequency,lag,psym=4,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
	        yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle=ystr,/xlog,/ylog
    	muploterr,frequency,lag,lag_err,/xlog,psym=3
    endelse
endelse
oplot,[xmin,xmax],[0,0],linestyle=2
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif
end
