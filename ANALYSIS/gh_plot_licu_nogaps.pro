pro gh_plot_licu_nogaps,time,rate,x1,x2,y1,y2,ps=psopt
;+
; NAME: 
;      GH_PLOT_LICU_NOGAPS
; PURPOSE: 
;      Plots a light curve to the current window
; EXPLANATION:
;      This procedure plots a light curve to the current window.
;
; CALLING SEQUENCE: 
;       GH_PLOT_LICU_NOGAPS,RATE,X1,X2,Y1,Y2[,/PS]
; INPUTS:
;       RATE     = Rate array
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
;       Plots a light curve with times from 0 to 1000s:
;
;       MU> GH_PLOT_LICU,t,x,0.0,1000.0
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
;       T. Belloni  12 May 2001  Adapted to MUFFT
;		T. Belloni  07 May 2010  from mu_plot_power
;		T. Belloni  17 Oct 2019  from mu_plot_licu
;-
;--------------------------------------------------------------------------
;
; To plot in a nice way a light curve (no error bars in Y)
;
timei = findgen(n_elements(rate))   ; time is the index now, from 0
n     = n_elements(time)
if (N_params() gt 2) then begin
   xmin = x1
   xmax = x2
  endif else begin
   xmin = timei(0); * 0.8
   xmax = timei(n-1); * 1.2
endelse
if (N_params() eq 6) then begin
   ymin = y1
   ymax = y2
  endif else begin
   ymin = min(rate);   * 0.8
   ymax = max(rate);   * 1.2
endelse
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_licu.ps'
endif
plot,timei,rate,psym=10,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
     yrange=[ymin,ymax],xtitle='PDS index',ytitle='Rate (cts/s)'
;
; Identify the gaps
;
delta  = time[2:n-1]-time[1:n-2]
deltat = min(delta)

gaps = where(delta gt (deltat*1.1))+0.5
ng = n_elements(gaps)
!p.linestyle=1
FOR i=0,(ng-1) DO BEGIN
	oplot,[gaps[i],gaps[i]],[ymin,ymax]
ENDFOR
!p.linestyle=0
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif
end
