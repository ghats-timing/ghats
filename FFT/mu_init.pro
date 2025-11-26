pro mu_init,filename,time_offset,sw_first_file, $
               fields_selected,jmeta,source,observatory, $
               instrument,mjdrefi,mjdreff,timezero
;
; Procedure to start the ball rolling for the accumulation. It is run
; on each input file, which is here opened then closed.
;
;-----------------------------------------------------------------------
;Parameters
; filename                 I: input filename
; time_offset              O: offset time (from header)
; sw_first_file            I: flag on whether it is the first file
; fields_selected          I/O: flag only to be set on file 1 (?)
; jmeta                    I: current index for metafile
; source                   O: output source name
; observatory              O: output satellite name
; instrument               O: output instrumant name
; mjdrefi:                 O: integer part of time_offset (minus t0)
; mjdreff                  O: fractional part of time_offset (minus t0)
; timezero                 O: timezero offset
;-----------------------------------------------------------------------
;
unit=10
fxbopen,unit,filename,1,header,errmsg=errmsg

datamode  =strtrim(fxpar(header,'DATAMODE'))
extname	  =strtrim(fxpar(header,'EXTNAME'))

if(extname eq 'XTE_HK') then begin
  massage,'HK data are not supported!'
   retall
endif

tstart=0.0
tstop=tstart

tstart  = fxpar(header,'TSTART')
tstop   = fxpar(header,'TSTOP')

timedel=double(0.0)
timedel = fxpar(header,'TIMEDEL')

mjdrefi = fxpar(header,'MJDREFI')
mjdreff = fxpar(header,'MJDREFF')
timezero= fxpar(header,'TIMEZERO')

time_offset = double(mjdrefi) + mjdreff + timezero/86400.0d0

tdim      = fxpar(header,'TDIM*')
nfields   = n_elements(tdim)
ttype     = fxpar(header,'TTYPE*')

observatory =strtrim(fxpar(header,'TELESCOP'))
;instrument  =strtrim(fxpar(header,'INSTRUME'))
source      =strtrim(fxpar(header,'OBJECT'))

if(sw_first_file eq 1) then begin
;   select_data,nfields,ttype,tdim,fields_selected(0,jmeta)  *****
;   here no user data selection is allowed. Only field number 2,
;                 the usual data field, is selected automatically

   fields_selected(0,jmeta) = 1
;   sw_first_file = 0   ; commented out in Michiel's version

endif
;
;  close input file
;
fxbclose,unit

end
