pro ghx_read_pds_metafile,filename,pdsfiles
;+
; NAME:
;      GHX_READ_PDS_METAFILE
; PURPOSE:
;      Read a GHX PDS metafile.
; EXPLANATION:
;      A GHX metafile is passed as '@file' and contains one PDS filename per
;      non-empty line. Entries are first used as written, matching other GHATS
;      metafile readers. If a relative entry is not found, GHX also tries it
;      relative to the metafile location.
;-

metafile = strmid(filename,1)
if(metafile eq '') then begin
   massage,'Empty GHX metafile name!'
   retall
endif

openr,unit,metafile,/get_lun,error=err
if(err ne 0) then begin
   massage,'Cannot open GHX metafile '+metafile
   retall
endif

meta_dir = file_dirname(metafile)
pdsfiles = strarr(1)
nfiles = 0L
line = ''

while(~eof(unit)) do begin
   readf,unit,line
   entry = strtrim(line,2)
   if(entry ne '') then begin
      first_char = strmid(entry,0,1)
      if((first_char ne '#') and (first_char ne ';')) then begin
         if((strpos(entry,'/') ne 0) and (meta_dir ne '') and $
            (meta_dir ne '.') and (~file_test(entry))) then begin
            entry = meta_dir+'/'+entry
         endif
         if(nfiles eq 0L) then begin
            pdsfiles[0] = entry
         endif else begin
            pdsfiles = [pdsfiles,entry]
         endelse
         nfiles = nfiles + 1L
      endif
   endif
endwhile

free_lun,unit

if(nfiles eq 0L) then begin
   massage,'GHX metafile contains no PDS filenames!'
   retall
endif

end
