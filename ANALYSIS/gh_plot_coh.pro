pro gh_plot_coh,frequency,coh,coh_err,x1,x2,y1,y2,ps=psopt
;+
; NAME: 
;      GH_PLOT_COH
; PURPOSE: 
;      Plots a coherence spectrum in log-log to the current window
; EXPLANATION:
;      This procedure plots a cohrence spectrum in log-log to 
;      the current window. 
;
; CALLING SEQUENCE: 
;       GH_PLOT_COH,FREQUENCY,COH,COH_ERR,X1,X2,Y1,Y2[,/PS]
; INPUTS:
;       FREQUENCY= Frequency array
;       COH      = Coherence array
;       COH_ERR  = Array of errors on coherence
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
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       None
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  01 Oct 2002  from mu_plot_cross
;       T. Belloni  03 Apr 2009  plot_oo removed
;	    T. Belloni  05 Apr 2009  gmu version adapted
;		T. Belloni  07 Jun 2009  fixed usersym problem under GDL
;		T. Belloni  26 Jan 2012  from mu_plot_coherence
;		T. Belloni  01 Mar 2012  changed default PS name
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
   ymin = 0.0
   ymax = 1.1
endelse
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_coh.ps'
endif
;
;plot,frequency,coh,psym=8,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
;     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle='Coherence', $
;     /xlog,/ylog

if(sistema eq 'IDL') then begin
	a=findgen(16)*(!PI*2/16.0)
	usersym,cos(a),sin(a),/fill
	plot,frequency,coh,psym=8,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
	     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle='Coherence', $
	     /xlog;,/ylog
   muerrplot,frequency,coh-coh_err,coh+coh_err,width=0.0
endif else begin
	plot,frequency,coh,psym=4,xstyle=1,ystyle=1,xrange=[xmin,xmax], $
	     yrange=[ymin,ymax],xtitle='Frequency (Hz)',ytitle='Coherence', $
	     /xlog;,/ylog
   muploterr,frequency,coh,coh_err,/xlog,psym=3
endelse
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif
end
