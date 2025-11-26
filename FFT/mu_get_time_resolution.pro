pro mu_get_time_resolution,filenames,nfiles,nmetafiles,tress
;
; Procedure to extract time resolution info
;
;------------------------------------------------------------------------
; Parameters
;
; filenames                 I: array of input filenames
; nfiles                    I: X dimension of filenames
; nmetafiles                I: Y dimension of filenames
; tress                     O: output time resolution
;------------------------------------------------------------------------
type=''
unit = 10
for j=0,nmetafiles-1 do begin
   fxbopen,unit,filenames(0,j),1,header,errmsg=errmsg
   instrument=strtrim(fxpar(header,'INSTRUME'))
   extname   =strtrim(fxpar(header,'EXTNAME'))
   datamode  =strtrim(fxpar(header,'DATAMODE'))
;
   if(extname eq 'XTE_HK') then begin
      type='xte_hk'
      massage,'HK data are not supported!'
      retall
   endif
   if((extname eq 'XTE_SA') and (instrument eq 'PCA')) then begin
      type='pca_sa'
   endif
   if((extname eq 'XTE_SE') and (instrument eq 'PCA')) then begin
      type='pca_se'
      if(datamode eq 'GoodXenon') then type = 'pca_gx'
   endif
   if((extname eq 'XTE_SA') and (instrument eq 'HEXTE')) then begin
      type='hxt_sa'
   endif
   if((extname eq 'XTE_SE') and (instrument eq 'HEXTE')) then begin
      type='hxt_se'
   endif

   if(type eq '') then begin
      print,type
      massage,'Unrecognized data file!'
      retall
   endif

   timedel   =fxpar(header,'TIMEDEL')
   tdim	     =fxpar(header,'TDIM*')
   ntdim     = n_elements(tdim)

   tress(0,j) = timedel
   modo = strmid(type,4,5)
   if(modo eq 'sa') then begin
      for k=0,ntdim-1 do begin
	 tdim(k)=strtrim(tdim(k))
	 if(tdim(k) ne '') then begin
	    if(strmid(tdim(k),0,1) eq '(') then begin
	       temp_string = strmid(tdim(k),1,strlen(tdim(k))-2)
	       ipos = strpos(temp_string,',')
	       nchans = 1
	       nhisto = 0
	       if (ipos ge 0) then begin
;                 GDL-compatibility modification  (19 Nov 2008 TMB)
                  temp_string2=STRJOIN(STRSPLIT(temp_string, /EXTRACT,','), ' ')
		  reads,temp_string2,nhisto,nchans
	       endif ELSE begin
		  reads,temp_string,nhisto
	       endelse
;
	       if(datamode eq 'Standard2') then begin
		  ntemp	 = nhisto
		  nhisto = nchans
		  nchans = ntemp
		  if (nhisto ne 1) then begin
		     massage,'NHISTO is not 1!!'
		  endif
	       endif
	       tress(0,j) = timedel/double(nhisto)
	       goto,seisei
	    endif
	 endif
      endfor
   endif

   seisei: fxbclose,unit	; close FITS file (can give problems?)
endfor
end
