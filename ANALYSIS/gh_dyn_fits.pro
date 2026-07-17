pro gh_dyn_fits, inputpds, outputfits,     $
                     index=index,frebin=reb,trebin=treb, arguments, $
                     show=show,help=help

; NAME:
;       GH_DYN_FITS (Federico Garcia)
; PURPOSE:
;       Runs GH_DYN to create a dynamical PDS from the INPUTPDS and 
;       store it into a FITS file together with TIME, FREQ and RATE.
;       You can use accompannying GH_dyn_fits.py to make the figures.
; UPDATES:
;       GH_DYN_FITS v2023/07/28 (Federico Garcia)
;             Added /show key to launch the python viewer GH_dyn_fits_tk.py
;       GH_DYN_FITS v2023/07/27 (Federico Garcia)
;             Original version
;
  compile_opt idl2  

  if keyword_set(help) then begin
       print,''
       print,'GH_DYN_FITS'
       print,''
       print,'Create a dynamic PDS with GH_DYN and save it to FITS.'
       print,''
       print,'Usage:'
       print,"  GH_DYN_FITS, 'input.pds', 'output.fits'"
       print,"  GH_DYN_FITS, 'input.pds', 'output.fits', INDEX=[i1,i2], FREBIN=-100, TREBIN=4"
       print,"  GH_DYN_FITS, 'input.pds', 'output.fits', /SHOW"
       print,''
       print,'Keywords:'
       print,'  INDEX=     Transform-index range passed to GH_DYN'
       print,'  FREBIN=    Frequency rebin factor passed to GH_DYN'
       print,'  TREBIN=    Time rebin factor passed to GH_DYN'
       print,'  ARGUMENTS= String of arguments passed to the Python viewer'
       print,'  /SHOW      Launch GH_dyn_fits_tk.py after writing the FITS file'
       print,'  /HELP      Print this message'
       print,''
       return
  endif

  if N_params() LT 2 then begin 
       print,'Syntax - gh_dyn_fits, inputpds, outputfits [,INDEX=INDEX][,FREBIN=FREBIN][,TREBIN=TREBIN][,ARGUMENTS][,/SHOW]'
       print,'ARGUMENTS=String of arguments passed to the Python plotter (must be in quotes)'
       print,'To see a list of possible arguments use xxx; the error message gives the valid options'
       return
  endif
  
  Catch, theError
  IF theError NE 0 then begin
	Catch,/Cancel
	void = cgErrorMsg(/quiet)
	RETURN
  ENDIF


if(keyword_set(index)) then begin
    starting = index[0]
    ending  = index[1]
endif else begin
    starting = 0L
    ending   = 100000L
endelse

if(keyword_set(reb)) then begin

endif else begin
    reb=1
endelse

if(keyword_set(treb)) then begin

endif else begin
    treb=1
endelse

strtreb=strtrim(string(treb),2)
if(arguments EQ !NULL) then arguments = ' '

; CALL GH_DYN
  gh_dyn,inputpds,time,rate,nu,dynima,index=index,frebin=reb,trebin=treb

; WRITE OUTPUT FITS FILE
  print,'Saving FITS file: ',outputfits
  writefits, outputfits, dynima
  writefits, outputfits, nu, /APPEND
  writefits, outputfits, time, /APPEND
  writefits, outputfits, rate, /APPEND

; PLOT THE FITS FILE USING GH_DYN_FITS_TK.PY
if(keyword_set(show)) then begin
  cmd = 'GH_dyn_fits_tk.py '+outputfits+' -trebin '+strtreb+' '+arguments+' '
  print,'Plotting with ',cmd,' ...'
  print,'Wait a moment; the plot is coming...'
  spawn,cmd
endif

  return
  end
