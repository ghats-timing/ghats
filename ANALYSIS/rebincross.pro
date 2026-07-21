pro rebincross,x,y,irf,xr,yr,nn
;+
; NAME: 
;      REBINCROSS
; PURPOSE: 
;      Rebins an array (for cross spectra) 
; EXPLANATION:
;      This procedure takes an array
;      and rebins it by a factor IRF. If IRF is negative, the rebinning
;      is logarythmic.
;
; CALLING SEQUENCE: 
;       REBINCROSS,X,Y,IRF,XR,YR,NN
; INPUTS:
;       X        = array of X values
;       Y        = array of Y values
;       IRF      = rebin factor. If IRF > 0, linear rebinning;
;                                if IRF < 0, log rebinning
;
; OUTPUTS:
;       XR       = array of rebinned X values
;       YR       = array of rebinned Y values
;       NN       = array with number of points rebinned per new bin
;
; KEYWORDS:
;       NONE
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       CROSSLOGREBXSYS: Log rebinning
; NOTES:
;       Logarythmic rebinning needs a negative rebin factor.
;       IRF = -20 means that each bin will have an amplitude
;                 increased by exp(1/20) with respect to the
;                 previous one.
;
;       The routine derives from murebin
;                 
; MODIFICATION HISTORY: 
;       T. Belloni  01 Oct 2002  from murebin
;		T. Belloni  03 Dec 2010  unchanged from Mu6
;		T. Belloni  05 Apr 2012  fixed nrd problem for irf=1
;-
;--------------------------------------------------------------------------
n = n_elements(x)
xr = fltarr(n)
yr = fltarr(n)

; **
; ** If rebin factor IRF equals 0 or 1 then transfer input directly to output *
; **
        if (irf eq 0 or irf eq 1) then begin
          ;;irf = 1
          xr = x
          yr = y
          nrd = n
          nn = intarr(n)*0+1
         endif else begin
          if (irf lt 0) then begin
; **
; ** Logarithmic rebin
; **
             crosslogrebxsys,x,y,xr,yr,-irf,nn,nrd
           endif else begin
             if (irf gt 1) then begin
                nrd = n / irf
                nn = intarr(nrd)*0+irf
; **
; ** Linear rebin
; **
                for i = 0l,nrd-1l do begin
                   xs  = 0.0
                   ys  = 0.0
                   sx2 = 0.0
                   sy2 = 0.0
                   is  = 0l
                   for j = 0l,irf-1l do begin
                      xs = xs + x((i)*irf + j)
                      ys = ys + y((i)*irf + j)
                      ;;is = is + 1l
                   ;;uno:
                   endfor
                  ;; if (is gt 0) then begin
                      xr(i)  = xs / float(irf)
                      yr(i)  = ys / float(irf)
                   ;;endif else begin
                      ;;xr(i) = 0.0
                      ;;yr(i) = 0.0
                   ;;endelse
                endfor
             endif
          endelse
       endelse

; trim output arrays
xr = xr(0:nrd-1)
yr = yr(0:nrd-1)
end
