pro mu_sort_time_resolution,filenames,metafiles,tress,nfiles,nmetafiles
;
; Procedure to sort filenames on time resolution
;
;------------------------------------------------------------------------
; Parameters
;
; filenames                  I/O: input list of filenames
; metafiles                  I/O: input list of metafiles
; tress                      I/O: input list of time resolutions
; nfiles                     I  : X dimension of filenames
; nmetafiles                 I  : Y dim of filenames and dim of metafiles
;------------------------------------------------------------------------

filename = strarr(nfiles)

for j=1,nmetafiles-1 do begin
   tres = tress(0,j)
   filename = filenames(*,j)
   metafile = metafiles(j)
   for k=j-1,0,-1 do begin
      if(tress(0,k) ge tres) then goto,trecinque
      tress(0,k+1) = tress(0,k)
      for i=0,nfiles-1 do begin
	 filenames(i,k+1) = filenames(i,k)
      endfor
      metafiles(k+1) = metafiles(k)
   endfor
   k = -1
   trecinque: tress(0,k+1) = tres
   filenames(*,k+1) = filename
   metafiles(k+1) = metafile
endfor

end
