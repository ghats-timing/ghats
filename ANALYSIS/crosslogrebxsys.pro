pro crosslogrebxsys,x,y,xr,yr,irf,nn,nrd
;+
; NAME: 
;      CROSSLOGREBXSYS
; PURPOSE: 
;      Rebins an array x,y logarythmically
; EXPLANATION:
;      This procedure takes an array with X and y 
;      and rebins it logarythmically by a factor abs(IRF).
;      This routine is intended for internal use only (REBINCROSS)
;
; CALLING SEQUENCE: 
;       CROSSLOGREBINXSYS,X,Y,XR,YR,IRF
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
;       NRD      = length of output arrays
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
;       T. Belloni  01 Oct 2002  from logrebxsys
;-
;--------------------------------------------------------------------------
n = n_elements(x)
nn = intarr(n)

        i = 0l
        nres = 0l
        xend = x(0)
        factor = 10.0^(1.0/(abs(float(irf))))
tre:    if (xend eq 0.0) then begin
          xr(i) = x(i)
          yr(i) = y(i)
          i     = i + 1l
          xend  = x(i)
          nres  = nres + 1l
          goto,tre
        endif
uno:    xend = x(i) * factor
        xs = x(i)
        ys = y(i)
        ns = 1l
        for j = 0l,n-1l do begin
          i = i + 1l
          if (x(i) gt xend or i ge (n-1)) then begin
            xr(nres) = xs / float(ns)
            yr(nres) = ys / float(ns)
            nn(nres) = ns
            nres = nres + 1l
            if (i ge (n-1)) then goto, due
            goto, uno
          endif else begin
            xs = xs + x(i)
            ys = ys + y(i)
            ns = ns + 1l
          endelse
        endfor
        xr(nres) = xs / float(ns)
        yr(nres) = ys / float(ns)
        nn(nres) = ns
        nres = nres + 1l
        goto,uno
due:    nrd = nres - 1l

nn = nn(0:nrd-1)
end
