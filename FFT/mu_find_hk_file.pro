pro mu_find_hk_file,housefile,tstart,tstop,filelist,nfiles
;
; Procedure to find the Standard2 file that covers midpoint of data bin
;
;--------------------------------------------------------------------------
; Parameters
;
; housefile                   O: selected standard2 file
; tstart                      I: start time of bin
; tstop                       I: end time of bin
; filelist                    I: list of housekeeping (std2) files
; nfiles                      I: length of filelist
;--------------------------------------------------------------------------
;
housefile=''
tavg = 0.5d0 *(tstart + tstop)

for i=0,nfiles-1 do begin
   unit=11
   fxbopen,unit,filelist(i),1,header,errmsg=errmsg
   t1   = fxpar(header,'TSTART')
   t2   = fxpar(header,'TSTOP')
   fxbclose,unit
   if((t1 le tavg) and (t2 ge tavg)) then begin
      housefile=filelist(i)
      return
   endif
endfor

end
