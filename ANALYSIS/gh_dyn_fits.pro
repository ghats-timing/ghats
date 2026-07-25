pro gh_dyn_fits, inputpds, outputfits,     $
                     index=index,frebin=reb,trebin=treb, arguments=arguments, $
                     show=show,help=help,cmhelp=cmhelp

; NAME:
;       GH_DYN_FITS (Federico Garcia)
; PURPOSE:
;       Runs GH_DYN to create a dynamical PDS from the INPUTPDS and 
;       store it into a FITS file together with TIME, FREQ and RATE.
;       You can use accompannying GH_dyn_fits.py to make the figures.
; UPDATES:
;       M. Mendez/Codex 2026/07/23
;             Moved the colormap catalogue from /HELP to /CMHELP.
;       M. Mendez/Codex 2026/07/23
;             Added compact colormap list to /HELP.
;       M. Mendez/Codex 2026/07/23
;             Expanded /HELP text for ARGUMENTS viewer options.
;       M. Mendez/Codex 2026/07/23
;             Made ARGUMENTS a real keyword and handled the unset case.
;       GH_DYN_FITS v2023/07/28 (Federico Garcia)
;             Added /show key to launch the python viewer GH_dyn_fits_tk.py
;       GH_DYN_FITS v2023/07/27 (Federico Garcia)
;             Original version
;
  compile_opt idl2  

  if keyword_set(cmhelp) then begin
       print,''
       print,'GH_DYN_FITS Colormaps'
       print,''
       print,'Use with ARGUMENTS, for example:'
       print,"  GH_DYN_FITS, file, 'dynpds.fits', /SHOW, ARGUMENTS='-cmap viridis'"
       print,"  GH_DYN_FITS, file, 'dynpds.fits', /SHOW, ARGUMENTS='-cmap viridis_r'"
       print,''
       print,'Append _r to reverse any colormap name.'
       print,''
       print,'Perceptual:'
       print,'  viridis, plasma, inferno, magma, cividis, turbo'
       print,''
       print,'Sequential:'
       print,'  Grays, Greys, Purples, Blues, Greens, Oranges, Reds'
       print,'  YlOrBr, YlOrRd, OrRd, PuRd, RdPu, BuPu, GnBu, PuBu'
       print,'  YlGnBu, PuBuGn, BuGn, YlGn, binary, gist_yarg, gist_gray'
       print,'  gray, bone, pink, spring, summer, autumn, winter, cool'
       print,'  Wistia, hot, afmhot, gist_heat, copper'
       print,''
       print,'Diverging:'
       print,'  PiYG, PRGn, BrBG, PuOr, RdGy, RdBu, RdYlBu, RdYlGn'
       print,'  Spectral, coolwarm, bwr, seismic, berlin, managua, vanimo'
       print,''
       print,'Cyclic:'
       print,'  twilight, twilight_shifted, hsv'
       print,''
       print,'Qualitative:'
       print,'  Pastel1, Pastel2, Paired, Accent, Dark2, Set1, Set2, Set3'
       print,'  tab10, tab20, tab20b, tab20c'
       print,''
       print,'Miscellaneous:'
       print,'  flag, prism, ocean, gist_earth, terrain, gist_stern'
       print,'  gnuplot, gnuplot2, CMRmap, cubehelix, brg, gist_rainbow'
       print,'  rainbow, jet, nipy_spectral, gist_ncar'
       print,''
       print,'Grey aliases: gray/grey, gist_gray/gist_grey, gist_yarg/gist_yerg'
       print,''
       return
  endif

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
       print,"  GH_DYN_FITS, 'input.pds', 'output.fits', /SHOW, $"
       print,"                ARGUMENTS='--linear -fmin 0.1 -fmax 50 -pmin 1.5 -pmax 4.0'"
       print,''
       print,'Keywords:'
       print,'  INDEX=     Transform-index range passed to GH_DYN'
       print,'  FREBIN=    Frequency rebin factor passed to GH_DYN'
       print,'  TREBIN=    Time rebin factor passed to GH_DYN'
       print,'  ARGUMENTS= String of arguments passed to the Python viewer'
       print,'  /SHOW      Launch GH_dyn_fits_tk.py after writing the FITS file'
       print,'  /HELP      Print this message'
       print,'  /CMHELP    Print the available colour maps for -cmap NAME'
       print,''
       print,'ARGUMENTS viewer options:'
       print,'  -title TEXT          Window title'
       print,'  -t0 MJD              Reference MJD for the time axis'
       print,'  -cmap NAME           Matplotlib colour map'
       print,'  -interp NAME         Image interpolation method'
       print,'  -pmin P, -pmax P     Colour-scale power limits'
       print,'  -fmin F, -fmax F     Frequency range to display'
       print,'  -ft F1 F2 ...        Frequency tick positions'
       print,'  -fh F1 F2 ...        Highlighted frequencies'
       print,'  -hc COLOR            Highlight colour'
       print,'  -tmin T, -tmax T     Time range to display'
       print,'  -tt T1 T2 ...        Time tick positions'
       print,'  -th T1 T2 ...        Highlighted times'
       print,'  --savepdf            Save figure as PDF'
       print,'  --savepng            Save figure as PNG'
       print,'  --noshow             Do not display the interactive window'
       print,'  --dynPDSonly         Hide the light curve panel'
       print,'  --nupnu              Plot nu*Pnu instead of Pnu'
       print,'  --linear             Use a linear colour scale'
       print,'  -cb, --colorbar      Show colour bar at startup'
       print,'  --fillGaps           Plot the full timespan including GTI gaps'
       print,'  --showGaps           Mark inferred GTI gaps with vertical lines'
       print,''
       print,'Use GH_DYN_FITS,/CMHELP to list available colour maps for -cmap NAME.'
       print,''
       print,'The -trebin option is added automatically from TREBIN=.'
       print,'For the Python viewer help from the shell: GH_dyn_fits_tk.py -h'
       print,''
       return
  endif

  if N_params() LT 2 then begin 
       print,'Syntax - gh_dyn_fits, inputpds, outputfits [,INDEX=INDEX][,FREBIN=FREBIN][,TREBIN=TREBIN][,ARGUMENTS][,/SHOW]'
       print,'ARGUMENTS=String of arguments passed to the Python plotter (must be in quotes)'
       print,'Use GH_DYN_FITS,/HELP to see the common ARGUMENTS options.'
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
if(n_elements(arguments) eq 0) then arguments = ' '

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
