pro gh_bispec_plot,filename,quantity,log=log,clip=clip,zrange=zrange,cmap=cmap, $
                    xlog=xlog,ylog=ylog, $
                    f1range=f1range,f2range=f2range, $
                    help=help
;+
; NAME:
;      GH_BISPEC_PLOT
;
; PURPOSE:
;      Displays one derived-product image extension from a GH_BISPEC_FITS or
;      GH_WRITE_BISPEC_FITS output file.
;
; CALLING SEQUENCE:
;      GH_BISPEC_PLOT,FILENAME,QUANTITY[,/LOG][,CLIP=[LO,HI]] $
;                     [,ZRANGE=[ZMIN,ZMAX]][,CMAP=CMAP][,/HELP]
;
; INPUTS:
;      FILENAME = multi-extension FITS file written by GH_BISPEC_FITS or
;                 GH_WRITE_BISPEC_FITS
;      QUANTITY = one of:
;                 'BREAL'
;                 'BIMAG'
;                 'BMOD'
;                 'BPHASE'
;                 'BICOH'
;                 'INT_BICOH'
;                 'RAW_DEN1_CORR'
;                 'RAW_DEN2_CORR'
;                 'NPROD_USED'
;
; KEYWORDS:
;      LOG    = if set, display log10(abs(image)).
;               Useful mostly for BMOD, BREAL and BIMAG.
;
;      CLIP   = two-element percentile range [lo,hi].
;               Default is [1,99].
;
;      F1RANGE = two-element range [f1min,f1max] to display.
;
;      F2RANGE = two-element range [f2min,f2max] to display.
;
;      ZRANGE = two-element range [zmin,zmax] for display scaling.
;               Overrides CLIP if supplied.
;
;      CMAP   = colour map name. Available names:
;               'gray'      = gray scale
;               'grey'      = gray scale
;               'rainbow'   = rainbow
;               'heat'      = red temperature
;               'bluewhite' = blue-white
;               'prism'     = prism
;               'redblue'   = blue-red, if available
;      XLOG   = force logarithmic x axis.
;
;      YLOG   = force logarithmic y axis.
;
;               If XLOG/YLOG are not given, the routine tries to read
;               AXLOG1/AXLOG2 from the FITS header.
;
;
;      HELP   = print usage message and return.
;
; EXAMPLE:
;      IDL> gh_bispec_plot,'out.fits','BICOH'
;      IDL> gh_bispec_plot,'out.fits','BMOD',/log,clip=[1,99.9],cmap='heat'
;      IDL> gh_bispec_plot,'out.fits','BICOH',/log,zrange=[-4,-2.5],cmap='heat'
;
; ROUTINES USED:
;      MRDFITS
;      SXPAR
;      LOADCT
;      CONTOUR
;      AXIS
;      XYOUTS
;
; MODIFICATION HISTORY:
;      M. Mendez  21 May 2026  first version, with help from ChatGPT
;-
;--------------------------------------------------------------------------

if(keyword_set(help)) then begin
   print,' '
   print,'GH_BISPEC_PLOT'
   print,'Display one derived-product image extension from a bispectrum FITS file.'
   print,' '
   print,'Usage:'
   print,"  gh_bispec_plot,'out.fits','BICOH',f1range=[0.1,20],f2range=[0.1,20]"
   print,"  gh_bispec_plot,'out.fits','BMOD',/log,clip=[1,99.9],cmap='heat'"
   print,"  gh_bispec_plot,'out.fits','BICOH',/log,zrange=[-4,-2.5],cmap='heat'"
   print,' '
   print,'Available quantities:'
   print,'  BREAL, BIMAG, BMOD, BPHASE, BICOH, INT_BICOH,'
   print,'  RAW_DEN1_CORR, RAW_DEN2_CORR, NPROD_USED'
   print,' '
   print,'  F1RANGE=[fmin,fmax]   frequency range to display on f1 axis'
   print,'  F2RANGE=[fmin,fmax]   frequency range to display on f2 axis'
   print,' '
   print,'Available colour maps:'
   print,"  'gray'       gray scale"
   print,"  'grey'       gray scale"
   print,"  'rainbow'    rainbow"
   print,"  'heat'       red temperature"
   print,"  'bluewhite'  blue-white"
   print,"  'prism'      prism"
   print,"  'redblue'    blue-red, if available"
   print,'  XLOG/YLOG   force logarithmic x/y axes'
   print,'              Otherwise AXLOG1/AXLOG2 are read from FITS header if present.'
   print,' '
   return
endif

if(n_params() lt 2) then begin
   massage,"Usage: gh_bispec_plot,filename,quantity"
   retall
endif

q = strupcase(strtrim(quantity,2))

case q of
   'BREAL': begin
      ext = 1
      qlabel = 'Real[B(f!D1!N,f!D2!N)]'
   end
   'BIMAG': begin
      ext = 2
      qlabel = 'Im[B(f!D1!N,f!D2!N)]'
   end
   'BMOD': begin
      ext = 3
      qlabel = '|B(f!D1!N,f!D2!N)|'
   end
   'BPHASE': begin
      ext = 4
      qlabel = 'Phase[B(f!D1!N,f!D2!N)]'
   end
   'BICOH': begin
      ext = 5
      qlabel = 'Bicoherence'
   end
   'INT_BICOH': begin
      ext = -1
      qlabel = 'Poisson-corrected intrinsic bicoherence'
      qfind = 'INT_BICOH'
   end
   'BICOH_INT_OBS': begin
      ext = -1
      qlabel = 'Poisson-corrected intrinsic bicoherence'
      qfind = 'BICOH_INT_OBS'
   end
   'RAW_DEN1_CORR': begin
      ext = -1
      qlabel = 'RAW DEN1 corrected'
      qfind = 'RAW_DEN1_CORR'
   end
   'RAW_DEN2_CORR': begin
      ext = -1
      qlabel = 'RAW DEN2 corrected'
      qfind = 'RAW_DEN2_CORR'
   end
   'NPROD_USED': begin
      ext = 6
      qlabel = 'N!Dprod!N'
   end
   else: begin
      massage,'Unknown quantity. Use BREAL, BIMAG, BMOD, BPHASE, BICOH, INT_BICOH, RAW_DEN1_CORR, RAW_DEN2_CORR, or NPROD_USED'
      retall
   end
endcase

;--------------------------------------------------------------------------
; Read image.
;--------------------------------------------------------------------------
if(ext lt 0) then begin
   primary_header = headfits(filename)
   nextend = long(fxpar(primary_header,'NEXTEND'))
   if(nextend le 0L) then nextend = 64L
   found_ext = 0B
   for iext=1L,nextend do begin
      test = mrdfits(filename,iext,thdr,/silent,status=status)
      if(status ne 0) then goto,done_find_ext
      extname = strupcase(strtrim(sxpar(thdr,'EXTNAME'),2))
      if((extname eq qfind) or ((qfind eq 'INT_BICOH') and (extname eq 'BICOH_INT_OBS'))) then begin
         img = test
         hdr = thdr
         found_ext = 1B
         goto,done_find_ext
      endif
   endfor
   done_find_ext:
   if(~found_ext) then begin
      massage,'Requested EXTNAME not found in FITS file: '+q
      retall
   endif
endif else begin
img = mrdfits(filename,ext,hdr)
endelse
disp = float(img)

;--------------------------------------------------------------------------
; Optional log10(abs(image)) display.
;--------------------------------------------------------------------------
if(keyword_set(log)) then begin
   w = where(finite(disp) and abs(disp) gt 0.0, nw)
   tmp = disp*0.0 + !VALUES.F_NAN
   if(nw gt 0) then tmp[w] = alog10(abs(disp[w]))
   disp = tmp
   qlabel = 'log!D10!N('+qlabel+')'
endif

show = disp

;--------------------------------------------------------------------------
; Reconstruct frequency axes.
;
; Prefer explicit FREQ1/FREQ2 extensions if present. These are needed for
; logarithmically rebinned files, because CRVAL/CDELT cannot fully describe
; non-uniform frequency spacing.
;
; Fall back to CRVAL/CDELT for older files.
;--------------------------------------------------------------------------
nx = (size(show))[1]
ny = (size(show))[2]

freq1 = dblarr(nx)
freq2 = dblarr(ny)

got_freq1 = 0
got_freq2 = 0

catch,error_status
if(error_status eq 0) then begin
   tmpfreq1 = mrdfits(filename,7,tmp_hdr1)
   if(n_elements(tmpfreq1) eq nx) then begin
      freq1 = double(tmpfreq1)
      got_freq1 = 1
   endif
endif
catch,/cancel

catch,error_status
if(error_status eq 0) then begin
   tmpfreq2 = mrdfits(filename,8,tmp_hdr2)
   if(n_elements(tmpfreq2) eq ny) then begin
      freq2 = double(tmpfreq2)
      got_freq2 = 1
   endif
endif
catch,/cancel

if(got_freq1 eq 0) then begin
   xmin = double(sxpar(hdr,'CRVAL1'))
   dx   = double(sxpar(hdr,'CDELT1'))
   freq1 = xmin + dx*dindgen(nx)
endif

if(got_freq2 eq 0) then begin
   ymin = double(sxpar(hdr,'CRVAL2'))
   dy   = double(sxpar(hdr,'CDELT2'))
   freq2 = ymin + dy*dindgen(ny)
endif

;--------------------------------------------------------------------------
; Determine whether to use log axes.
;
; Manual /XLOG and /YLOG override header defaults.
; Otherwise use AXLOG1/AXLOG2 from the FITS header.
;--------------------------------------------------------------------------
axlog1 = long(sxpar(hdr,'AXLOG1'))
axlog2 = long(sxpar(hdr,'AXLOG2'))

do_xlog = axlog1
do_ylog = axlog2

if(keyword_set(xlog)) then do_xlog = 1
if(keyword_set(ylog)) then do_ylog = 1

;--------------------------------------------------------------------------
; Optional frequency cropping for display.
;--------------------------------------------------------------------------
if(keyword_set(f1range)) then begin

   if(n_elements(f1range) ne 2) then begin
      massage,'F1RANGE must be [f1min,f1max]'
      retall
   endif

   w1 = where(freq1 ge f1range[0] and $
              freq1 le f1range[1], n1)

   if(n1 le 0) then begin
      massage,'No f1 bins found in F1RANGE'
      retall
   endif

endif else begin

   w1 = lindgen(n_elements(freq1))
   n1 = n_elements(freq1)

endelse

if(keyword_set(f2range)) then begin

   if(n_elements(f2range) ne 2) then begin
      massage,'F2RANGE must be [f2min,f2max]'
      retall
   endif

   w2 = where(freq2 ge f2range[0] and $
              freq2 le f2range[1], n2)

   if(n2 le 0) then begin
      massage,'No f2 bins found in F2RANGE'
      retall
   endif

endif else begin

   w2 = lindgen(n_elements(freq2))
   n2 = n_elements(freq2)

endelse

; Crop image and axes.
i1lo = w1[0]
i1hi = w1[n1-1]
i2lo = w2[0]
i2hi = w2[n2-1]

show  = show[i1lo:i1hi,i2lo:i2hi]
freq1 = freq1[i1lo:i1hi]
freq2 = freq2[i2lo:i2hi]

; Log axes require positive frequencies.
if(do_xlog ne 0) then begin
   if(min(freq1) le 0.0d0) then begin
      print,'Cannot use logarithmic x axis: FREQ1 contains values <= 0.'
      do_xlog = 0
   endif
endif

if(do_ylog ne 0) then begin
   if(min(freq2) le 0.0d0) then begin
      print,'Cannot use logarithmic y axis: FREQ2 contains values <= 0.'
      do_ylog = 0
   endif
endif

wfinite = where(finite(show), nfinite)
if(nfinite le 0) then begin
   massage,'No finite pixels to display'
   retall
endif

if(keyword_set(zrange)) then begin

   if(n_elements(zrange) ne 2) then begin
      massage,'ZRANGE must be [zmin,zmax]'
      retall
   endif

   dlo = float(zrange[0])
   dhi = float(zrange[1])

endif else begin

   if(keyword_set(clip)) then begin
      clo = clip[0]
      chi = clip[1]
   endif else begin
      clo = 1.0
      chi = 99.0
   endelse

   vals = show[wfinite]
   vals = vals[sort(vals)]
   nvals = n_elements(vals)

   ilo = long((clo/100.0)*(nvals-1)) > 0 < (nvals-1)
   ihi = long((chi/100.0)*(nvals-1)) > 0 < (nvals-1)

   dlo = vals[ilo]
   dhi = vals[ihi]

   if(dhi le dlo) then begin
      dlo = min(vals)
      dhi = max(vals)
   endif

endelse

if(dhi le dlo) then begin
   massage,'Image has no useful dynamic range'
   retall
endif

;--------------------------------------------------------------------------
; Treat NaN/Inf pixels as masked pixels.
;
; IDL CONTOUR,/FILL does not reliably leave NaNs blank. Therefore invalid
; pixels are replaced by a sentinel below the displayed data range. The first
; contour interval is assigned colour index 0, which we later force to white.
;--------------------------------------------------------------------------
masked_color_index = 0
sentinel = dlo - 1000.0*abs(dhi-dlo)

wbad = where(~finite(show), nbad)

if(nbad gt 0) then begin
   show[wbad] = sentinel
endif

; Clip only valid data. Leave sentinel pixels untouched.
wgood = where(show ne sentinel, ngood)

if(ngood gt 0) then begin
   show[wgood] = show[wgood] > dlo < dhi
endif

;--------------------------------------------------------------------------
; Levels and colours for filled contours.
; The first interval is reserved for masked pixels.
;--------------------------------------------------------------------------
nlev_data = 64
levels_data = dlo + (dhi-dlo)*findgen(nlev_data)/(nlev_data-1)

if(nbad gt 0) then begin

   levels = [sentinel, dlo, levels_data[1:*]]
   colors = [masked_color_index, $
             fix(1.0 + 254.0*findgen(nlev_data-1)/(nlev_data-2))]

endif else begin

   levels = levels_data
   colors = fix(1.0 + 254.0*findgen(nlev_data-1)/(nlev_data-2))

endelse


;--------------------------------------------------------------------------
; Choose colour table.
; This must be done after opening the graphics device/window on some systems.
;--------------------------------------------------------------------------
if(keyword_set(cmap)) then begin
   cm = strlowcase(strtrim(cmap,2))
endif else begin
   cm = 'rainbow'
endelse

;--------------------------------------------------------------------------
; Save current plotting parameters.
;--------------------------------------------------------------------------
old_position = !p.position
old_charsize = !p.charsize
old_font     = !p.font

!p.charsize = 1.2
!p.font     = -1
black = 255

;--------------------------------------------------------------------------
; Open window and force indexed colours so LOADCT works.
;--------------------------------------------------------------------------
if(!d.name eq 'X') then window,title='GH_BISPEC_PLOT: '+q,xsize=900,ysize=760
device,decomposed=0
erase

case cm of
   'gray':      loadct,0,/silent
   'grey':      loadct,0,/silent
   'rainbow':   loadct,13,/silent
   'heat':      loadct,3,/silent
   'bluewhite': loadct,1,/silent
   'prism':     loadct,6,/silent
   'redblue':   loadct,33,/silent
   else: begin
      print,'Unknown CMAP. Using rainbow.'
      cm = 'rainbow'
      loadct,13,/silent
   end
endcase

;--------------------------------------------------------------------------
; Reserve colour index 0 for masked pixels.
;--------------------------------------------------------------------------
tvlct,rct,gct,bct,/get
rct[0] = 255B
gct[0] = 255B
bct[0] = 255B
tvlct,rct,gct,bct

; Force black for axes/text
rct[255] = 0B
gct[255] = 0B
bct[255] = 0B
tvlct,rct,gct,bct
black = 255

;--------------------------------------------------------------------------
; Main image.
;--------------------------------------------------------------------------
!p.position = [0.10,0.11,0.80,0.90]

contour,show,freq1,freq2,/fill, $
        levels=levels, $
        c_colors=colors, $
        xtitle='f!D1!N (Hz)', $
        ytitle='f!D2!N (Hz)', $
        title=qlabel, $
        xstyle=1, $
        ystyle=1, $
        xlog=do_xlog, $
        ylog=do_ylog, $
        color=black

;--------------------------------------------------------------------------
; Colour bar using CONTOUR, not TV.
;--------------------------------------------------------------------------
!p.position = [0.84,0.11,0.88,0.90]

barx = [0.0,1.0]
bary = dlo + (dhi-dlo)*findgen(256)/255.0
barz = fltarr(2,256)

for i=0,255 do begin
   barz[0,i] = bary[i]
   barz[1,i] = bary[i]
endfor

bar_levels = dlo + (dhi-dlo)*findgen(64)/(64-1)
bar_colors = fix(1.0 + 254.0*findgen(63)/(63-1))

contour,barz,barx,bary,/fill,/noerase, $
        levels=bar_levels, $
        c_colors=bar_colors, $
        xstyle=4, $
        ystyle=4, $
        color=black

axis,yaxis=1,yrange=[dlo,dhi],ystyle=1,charsize=1.2,color=black

;--------------------------------------------------------------------------
; Colour-bar label on the right. Some IDL/X combinations stack rotated
; strings character by character. If that happens, comment this XYOUTS out
; and use a horizontal label above the bar instead.
;--------------------------------------------------------------------------
xyouts,0.93,0.50,qlabel, $
        /normal, $
        orientation=270.0, $
        charsize=1.2, $
        charthick=2, $
        alignment=0.5, $
        color=black

;--------------------------------------------------------------------------
; Restore plotting parameters.
;--------------------------------------------------------------------------
!p.position = old_position
!p.charsize = old_charsize
!p.font = old_font

print,'Displayed ',q
if(keyword_set(zrange)) then begin
   print,'Z range: ',dlo,' - ',dhi
endif else begin
   print,'Clip percentiles: ',clo,' - ',chi
endelse
print,'Display range: ',dlo,'  ',dhi
print,'X axis log: ',do_xlog
print,'Y axis log: ',do_ylog
if(keyword_set(f1range)) then $
   print,'Displayed F1 range: ',f1range[0],' - ',f1range[1]

if(keyword_set(f2range)) then $
   print,'Displayed F2 range: ',f2range[0],' - ',f2range[1]

print,'Colour map: ',cm

end
