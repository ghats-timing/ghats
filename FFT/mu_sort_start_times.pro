pro mu_sort_start_times,     $
                 filenames,tstarts,nfiles,nmetafiles
;
; Sort filenames according to tstart
;
;-------------------------------------------------------------------------
; Parameters
;
; filenames                    I/O: array with input filenames
; nfiles                       I:   X dimension of filenames
; tstarts                      I/O: array of start times from keywords
; nmetafiles                   I:   Y dimension of filenames
;-------------------------------------------------------------------------
;
; A simple sorting routine
;
for j=0,nmetafiles-1 do begin
   for i=1,nfiles-1 do begin

      tstart   = tstarts(i,j)
      filename = filenames(i,j)
      for k=i-1,0,-1 do begin
	 if(tstarts(k,j) le tstart) then goto, salto
	 tstarts(k+1,j)	  = tstarts(k,j)
	 filenames(k+1,j) = filenames(k,j)
      endfor
      k=-1
      salto: tstarts(k+1,j) = tstart
      filenames(k+1,j)	    = filename
   endfor
endfor

end
