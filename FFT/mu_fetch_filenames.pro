pro mu_fetch_filenames,infilename,MAX_N_FILES,MAX_N_METAFILES, $
              filenames,metafiles,filenameroot,nfiles,nmetafiles
;
; Obtain input file names, dealing with metafiles and metametafiles
;   MU 3.0 version, system independent. Removed Linux spawn lines
;   MU 4.0 version: mu_get_lines replaces get_lines
;
;-------------------------------------------------------------------------
; Parameters
;
; infilename            I: name of input file name (or meta or meta2 file)
; MAX_N_FILES           I: input constant with maximum # of files
; MAX_N_METAFILES       I: input constant with maximum # of metafiles
;
; filenames             O: array with list of filenames
; metafiles             O: array with list of metafilenames
; filenameroot          O: root for output filenames
; nfiles                O: number of actual files
; nmetafiles            O: number of actual metafiles
;-------------------------------------------------------------------------

   sdum = ''
;
; Check whether the input file is a metafile (or metametafile)
;
   iposat = strpos(infilename,'@')
   if (iposat ge 0) then begin
;
;  filename starts with @: it's a (meta)metafile
;
      nmetafiles = 1
      metafiles(0) = infilename
      iposat2 = strpos(infilename,'@@')
      if(iposat2 ge 0) then begin
;
;    filename starts with @@: it's a metametafile
;
;
;        Open the file
;
     openr,1,infilename,ERROR=ierr
     if (ierr ne 0) then begin
        massage,'Input metametafile not found!'
        retall
     endif
;        Compute file length --> goes to nn
;        New system independent version in mu 3.0
;
       nn = mu_file_lines(infilename)
;
;        Read in metafile names within metametafile
;
     nmetafiles = 0
     for j=1,nn do begin
        nmetafiles = nmetafiles  +1
        readf,1,sdum
        metafiles(nmetafiles-1) = sdum
        iposat3 = strpos(sdum,'@')
        if (iposat3 lt 0) then begin
           massage,'File in metametafile is not a metafile!'
           retall
        endif
     endfor
     close,1
      endif
;
;    loop through the metafiles
;
      for j=0,nmetafiles-1 do begin
;
;        Compute number of lines in metafile
;        New system independent version
;
       nn = mu_file_lines(metafiles(j))
;
;        Open the metafile
;
     openr,1,metafiles(j),ERROR=ierr
     if (ierr ne 0) then begin
        massage,'Input metafile not found!'
        retall
     endif
;
;        Read in filenames from metafile
;
     nfiles = -1
     for k=1,nn do begin
        nfiles=nfiles+1
        readf,1,sdum
        filenames(nfiles,j) = sdum
     endfor
     if(j gt 0)  then begin
        if (nfiles ne nfiles_prev) then begin
           massage,'Metafiles not consistent!'
           retall
        endif
     endif
     nfiles_prev = nfiles
     nfiles = nfiles + 1  ; increases to have real n. of files
     close,1
      endfor
     endif ELSE begin
;
;   no metafile case (the simplest one)
;
      nfiles     = 1
      nmetafiles = 1
      metafiles(0) = 'NONE'
      filenames(0,0)=infilename
   endelse
;
;  extract root for output filename
;
if(iposat lt 0) then begin
;
;  No metafile case
;
   filenameroot = strtrim(infilename)
  endif ELSE begin
;
;  Metafile case
;
   ini=max([iposat+1,iposat2+2])
   filenameroot = strtrim(strmid(infilename,ini))
   punto=strpos(filenameroot,'.')
   if(punto ge 0) then begin
      filenameroot = strmid(filenameroot,0,punto)
   endif
endelse
;  remove directory info if there
; System-free call

filenameroot=file_basename(filenameroot,'.evt')
;filenameroot=fsc_base_filename(filenameroot)    ; gdl compatible version

end
