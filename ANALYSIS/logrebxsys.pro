pro logrebxsys,x,sx,y,sy,n,xr,sxr,yr,syr,nrd,irf
;+
; NAME: 
;      LOGREBXSYS
; PURPOSE: 
;      Rebins a power spectrum (with errors) logarythmically
; EXPLANATION:
;      This procedure takes an array with (optional) errors in X and/or Y
;      and rebins it logarythmically by a factor abs(IRF).
;      This routine is intended for internal use only (LOGREBIN)
;
; CALLING SEQUENCE: 
;       LOGREBINXSYS,X,SX,Y,SY,XR,SXR,YR,SYR,NRD,IRF
; INPUTS:
;       X        = array of X values
;       SX       = array of errors on X; if no errors are available, 
;                  this variable is dummy (use X)
;       Y        = array of Y values
;       SY       = array of errors on Y; if no errors are available, 
;                  this variable is dummy (use Y)
;       IRF      = rebin factor. If IRF > 0, linear rebinning;
;                                if IRF < 0, log rebinning
;
; OUTPUTS:
;       XR       = array of rebinned X values
;       SXR      = array of errors on SXR
;       YR       = array of rebinned Y values
;       SYR      = array of errors on SYR
;       NRD      = number of elements of XR and YR
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
;       None
; NOTES:
;       IRF = -20 means that each bin will have an amplitude
;                 increased by exp(1/20) with respect to the
;                 previous one.
;                 
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  10 Nov 2001  switch to long integers
;-
;--------------------------------------------------------------------------

swsx = max(sx)
swsy = max(sy)

        i = 0l
        nres = 0l
        xend = x(0)
        factor = 10.0^(1.0/(abs(float(irf))))
tre:    if (xend eq 0.0) then begin
          xr(i) = x(i)
          yr(i) = y(i)
          if (swsx gt 0.0) then sxr(i) = sx(i)
          if (swsy gt 0.0) then syr(i) = sy(i)
          i     = i + 1l
          xend  = x(i)
          nres  = nres + 1l
          goto,tre
        endif
uno:    xend = x(i) * factor
        xs = x(i)
        ys = y(i)
        if (swsx) then sx2 = sx(i)^2
        if (swsy) then sy2 = sy(i)^2
        ns = 1l
        for j = 0l,n-1l do begin
          i = i + 1l
          if (x(i) gt xend or i ge (n-1)) then begin
            xr(nres) = xs / float(ns)
            yr(nres) = ys / float(ns)
            if (swsx) then sxr(nres) = sqrt(sx2) / float(ns)
            if (swsy) then syr(nres) = sqrt(sy2) / float(ns)
            nres = nres + 1l
            if (i ge (n-1)) then goto, due
            goto, uno
          endif else begin
            xs = xs + x(i)
            ys = ys + y(i)
            if (swsx) then sx2 = sx2 + sx(i)^2
            if (swsy) then sy2 = sy2 + sy(i)^2
            ns = ns + 1l
          endelse
        endfor
        xr(nres) = xs / float(ns)
        yr(nres) = ys / float(ns)
        if (swsx) then sxr(nres) = sqrt(sx2) / float(ns)
        if (swsy) then syr(nres) = sqrt(sy2) / float(ns)
        nres = nres + 1l
        goto,uno
due:    nrd = nres - 1l

end
