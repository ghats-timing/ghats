pro ghrebin,x,sx,y,sy,irf,      $
            xr,sxr,yr,syr,nrd
;+
; NAME: 
;      GHREBIN
; PURPOSE: 
;      Rebins a power spectrum (with errors) either linearly or logarythmically
; EXPLANATION:
;      This procedure takes an array with (optional) errors in X and/or Y
;      and rebins it by a factor NRD. If NRD is negative, the rebinning
;      is logarythmic.
;
; CALLING SEQUENCE: 
;       GHREBIN,X,SX,Y,SY,IRF,XR,SXR,YR,SYR,NRD
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
;       LOGREBXSYS: Log rebinning
; NOTES:
;       Logarythmic rebinning needs a negative rebin factor.
;       IRF = -20 means that each bin will have an amplitude
;                 increased by exp(1/20) with respect to the
;                 previous one.
;
;       The routine is a direct port of a F77 rebin routine  by Mariano
;       Mendez.
;                 
; MODIFICATION HISTORY: 
;       T. Belloni  20 Aug 2001  implementation
;       T. Belloni  10 Nov 2001  switch to long integers
;		T. Belloni  06 May 2010  from murebin  
;-
;--------------------------------------------------------------------------
;
; Logarithmic rebinning after Mariano's FORTRAN subroutine
;

; Determine whether the error arrays are filled

swsx = max(sx)
swsy = max(sy)

n = n_elements(x)
xr = dblarr(n)
sxr = dblarr(n)
yr = dblarr(n)
syr = dblarr(n)

; **
; ** If rebin factor IRF equals 0 or 1 then transfer input directly to output *
; **
        if (irf eq 0 or irf eq 1) then begin
          irf = 1
          nrd = -1l
          for i = 0l,n-1l do begin
             if (not(swsx gt 0.0 and sx(i) lt 0.0)    $
                    and  not (swsy gt 0.0 and sy(i) lt 0.0)) then begin
                nrd = nrd + 1l
                xr(nrd) = x(i)
                yr(nrd) = y(i)
                if (swsx gt 0.0) then begin
                  sxr(nrd) = sx(i)
                 endif else begin
                  sxr(nrd) = 0.0
                endelse
                if (swsy gt 0.0) then begin
                  syr(nrd) = sy(i)
                 endif else begin
                  syr(nrd) = 0.0
                endelse
             endif
          endfor
         endif else begin
          if (irf lt 0) then begin
; **
; ** Logarthmic rebin
; **
             logrebxsys,x,sx,y,sy,n,xr,sxr,yr,syr,nrd,-irf
           endif else begin
             if (irf gt 1) then begin
                nrd = n / irf
; **
; ** Linear rebin, ignoring points with negative sigma's **
; **
                for i = 0l,nrd-1l do begin
                   xs  = 0.0
                   ys  = 0.0
                   sx2 = 0.0
                   sy2 = 0.0
                   is  = 0l
                   for j = 0l,irf-1l do begin
                      if (swsx gt 0.0 and sx((i) * irf + j) lt 0.0)then goto,uno
                      if (swsy gt 0.0 and sy((i) * irf + j) lt 0.0)then goto,uno
                      xs = xs + x((i)*irf + j)
                      ys = ys + y((i)*irf + j)
                      if (swsx gt 0.0) then sx2 = sx2 + sx((i) * irf + j)^2
                      if (swsy gt 0.0) then sy2 = sy2 + sy((i) * irf + j)^2
                      is = is + 1l
                   uno:
                   endfor
                   if (is gt 0) then begin
                      xr(i)  = xs / float(is)
                      yr(i)  = ys / float(is)
                      if (swsx gt 0.0) then begin
                         sxr(i) = sqrt(sx2) / float(is)
                        endif else begin
                         sxr(i) = 0.0
                      endelse
                      if (swsy gt 0.0) then begin
                         syr(i) = sqrt(sy2) / float(is)
                        endif else begin
                         syr(i) = 0.0
                      endelse
                   endif else begin
                      xr(i) = 0.0
                      yr(i) = 0.0
                      if (swsx gt 0.0) then sxr(i) = -1.0
                      if (swsy gt 0.0) then syr(i) = -1.0
                   endelse
                endfor
             endif
          endelse
       endelse

; trim output arrays
xr = xr(0:nrd-1)
sxr = sxr(0:nrd-1)
yr = yr(0:nrd-1)
syr = syr(0:nrd-1)
end
