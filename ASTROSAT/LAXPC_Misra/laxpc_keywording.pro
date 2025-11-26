pro laxpc_keywording,namelist

nlines=file_lines(namelist)
OPENR, lun, namelist, /GET_LUN
  osservazioni = STRARR(1, nlines)
  READF, lun, osservazioni
CLOSE, lun
FREE_LUN, lun
nlines=file_lines(namelist)

FOR i=0,nlines-1 DO BEGIN
   filename = osservazioni[i]
   primary_header=headfits(filename)
   ;;fxbopen,u,filename,1,primary_header ; read in primary header
   ;;fxbclose,u
   ; extract keywords from primary header
   timezero = fxpar(primary_header,'TIMEZERO')
   tstart   = fxpar(primary_header,'TSTART')
   tstop    = fxpar(primary_header,'TSTOP')
   object   = fxpar(primary_header,'OBJECT')
   observer = fxpar(primary_header,'OBSERVER')
   ;timedel  = 1.0e-5
   
   fxhmodify,filename,'TIMEZERO',timezero ,EXTENSION=1
   fxhmodify,filename,'TSTART'  ,tstart   ,EXTENSION=1
   fxhmodify,filename,'TSTOP'   ,tstop    ,EXTENSION=1
   fxhmodify,filename,'OBJECT'  ,object   ,EXTENSION=1
   fxhmodify,filename,'OBSERVER',observer ,EXTENSION=1
   ;fxhmodify,filename,'TIMEDEL' ,timedel  ,EXTENSION=1
   
   ;fxbopen,u,filename,2,extension_header ; read in extension header
   ;fxaddpar,extension_header,'TIMEZERO',timezero
   ;fxaddpar,extension_header,'TSTART',tstart  
   ;fxaddpar,extension_header,'TSTOP', tstop   
   ;fxaddpar,extension_header,'OBJECT',object  
   ;fxaddpar,extension_header,'OBSERVER',observer

ENDFOR
END