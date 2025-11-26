pro gh_plot_hk_xte,filename,ps=psopt
;+
; NAME: 
;      GH_PLOT_HK_XTE
; PURPOSE: 
;      Plot HK time series
; EXPLANATION:
;      This procedure reads the XTE/PCA HK information from a PDS file
;      and plots the resulting series to the terminal. No output if
;      produced.
;
; CALLING SEQUENCE: 
;       GH_PLOT_HK_XTE,FILENAME
; INPUTS:
;       FILENAME = name of the input PDS file
;
; OUTPUTS:
;       NONE
;
; KEYWORDS:
;       PS       = If set, output goes to a PS file
;
; EXAMPLE:
;       None
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       MUXANA_HK: Get HK information from a XTE/PCA PDS file
; NOTES:
;       None
; MODIFICATION HISTORY: 
;       T. Belloni  10 Nov 2001  implementation, rough version
;       T. Belloni  19 Nov 2001  Better multi-window plot
;       T. Belloni  21 Nov 2001  PS option
;       T. Belloni  12 Mag 2002  adapted for MUFFT
;		T. Belloni	01 Dec 2010  from mu6
;		T. Belloni  22 May 2012  changed name from GH_PLOT_HK
;		T. Belloni  19 Dec 2013  gh_hk_xte call corrected
;-
;--------------------------------------------------------------------------
;
;  Get HK info from file
;
gh_hk_xte,filename,times,vle,ndet,poiss
notick=replicate(' ',30)
;
;  PS output setup
;
if(keyword_set(psopt)) then begin
   entry_device=!d.name
   set_plot,'PS'
   device,filename='gh_hk_xte.ps'
endif

plot,times,vle,charsize=1.0,ytitle='VLE rate',psym=10, $
                            ystyle=16,xtickname=notick,$
                            position=[0.15,0.68,0.95,0.95]
plot,times,ndet,charsize=1.0,ytitle='# det.', $
                            yrange=[-0.5,5.5],psym=10,ystyle=1,xtickname=notick,$
                            position=[0.15,0.41,0.95,0.68],/noerase
plot,times,poiss,charsize=1.0,xtitle='Time (s)',ytitle='Poiss. level',psym=10, $
                            ystyle=17,/noerase,position=[0.15,0.14,0.95,0.41]
;
;  PS output close
;
if(keyword_set(psopt)) then begin
   device,/close_file
   set_plot,entry_device
endif

end
