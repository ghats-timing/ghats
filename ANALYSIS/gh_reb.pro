pro gh_reb,x,y,sy,irf,      $
            xr,yr,syr,help=help
;+
; NAME: 
;      GH_REB
; PURPOSE: 
;      Rebins a power spectrum (with errors) either linearly or logarythmically
; EXPLANATION:
;      This procedure takes arrays with frequency, power and power error
;      and rebins them by a factor NRD. If NRD is negative, the rebinning
;      is logarythmic.
;
; CALLING SEQUENCE: 
;       GH_REB,X,Y,SY,IRF,XR,YR,SYR[,/HELP]
; INPUTS:
;       X        = array of X values
;       Y        = array of Y values
;       SY       = array of errors on Y
;       IRF      = rebin factor. If IRF > 0, linear rebinning;
;                                if IRF < 0, log rebinning
;
; OUTPUTS:
;       XR       = array of rebinned X values
;       YR       = array of rebinned Y values
;       SYR      = array of errors on SYR
;
; KEYWORDS:
;       HELP     = If set, print usage information and return
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS: 
;       None 
; ROUTINES USED: 
;       GHREBIN: general rebin routine (from M. Mendez)
; NOTES:
;       Logarythmic rebinning needs a negative rebin factor.
;       IRF = -20 means that each bin will have an amplitude
;                 increased by exp(1/20) with respect to the
;                 previous one.
;
; MODIFICATION HISTORY: 
;       T. Belloni  16 Jul 2002  from murebin
;		T. Belloni  06 May 2010  from mu_reb
;-
;--------------------------------------------------------------------------
if(keyword_set(help)) then begin
   print,''
   print,'GH_REB'
   print,''
   print,'Rebin arrays with errors, linearly or logarithmically.'
   print,''
   print,'Usage:'
   print,'  GH_REB, x, y, yerr, irf, xr, yr, yrerr'
   print,''
   print,'Arguments:'
   print,'  x,y,yerr   Input x values, y values, and errors'
   print,'  irf        Rebin factor: >0 linear, <0 logarithmic'
   print,'  xr,yr,yrerr Output rebinned arrays'
   print,''
   print,'Example:'
   print,'  GH_REB, freq, pow, pow_err, -100, fr, pr, pre'
   print,''
   return
endif
;
ghrebin,x,x,y,sy,irf,xr,dum,yr,syr,nrd

end
