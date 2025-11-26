pro mu_get_channel_information,                      $
             filenames,                $
             channels1,channels2,nchannels,nhistos, $
             ioffset,ibit_chan,            $
         nfiles,nmetafiles
;
; Procedure for obtaining channel information from input files
;
;----------------------------------------------------------------------------
; Parameters
;
; filenames                   I: array of input filenames
; channels1                   O: start channel list for each file
; channels2                   O: end channel list for each file
; nchannels                   O: array of number of channels
; nhistos                     O: number of histograms (?)
; ioffset                     O: array of offsets (for event data)
; ibit_chan                   O: array of bits per chan (for event data)
; nfiles                      I: input number of files
; nmetafiles                  I: input number of metafiles
;----------------------------------------------------------------------------
;
ioffset=intarr(nfiles,nmetafiles)
nhistos=intarr(nfiles,nmetafiles)
type=''
nhisto=0
ichans = intarr(n_elements(channels1(*,0,0)))
ichane = ichans
unit = 10
for j=0,nmetafiles-1 do begin
   for i=0,nfiles-1 do begin

      fxbopen,unit,filenames(i,j),1,header,errmsg=errmsg
      instrument=strtrim(fxpar(header,'INSTRUME'))
      datamode  =strtrim(fxpar(header,'DATAMODE'))
      extname   =strtrim(fxpar(header,'EXTNAME'))

      if(extname eq 'XTE_HK') then begin
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
      timedel= fxpar(header,'TIMEDEL')
      tdim   = fxpar(header,'TDIM*')
      ntdim     = n_elements(tdim)

      modo = strmid(type,4,2)
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
;           GDL-compatibility modification  (19 Nov 2008 TMB)
            temp_string2=STRJOIN(STRSPLIT(temp_string, /EXTRACT,','), ' ')
            reads,temp_string2,nhisto,nchans
         endif ELSE begin
            reads,temp_string,nhisto
         endelse
;
         if(datamode eq 'Standard2') then begin
            ntemp  = nhisto
            nhisto = nchans
            nchans = ntemp
            if (nhisto ne 1) then begin
         massage,'NHISTO is not 1!!'
            endif
         endif
         goto,trecinque
           endif
        endif
     endfor
      endif
      trecinque:
      if(modo eq 'sa') then begin
         tddes=fxpar(header,'TDDES*')
         if(strmid(type,0,3) eq 'pca') then event_description = tddes(1)
         if(strmid(type,0,3) eq 'hxt') then event_description = tddes(3)
         mu_parse_channel_descriptor_se,event_description,nchans,ichans,ichane
        endif ELSE begin
         tevtb=fxpar(header,'TEVTB*')
         if(strmid(type,0,3) eq 'pca') then event_description = tevtb(1)
         if(strmid(type,0,3) eq 'hxt') then event_description = tevtb(3)
;
;        ioffset and ibit_chan elements must be passed by reference,
;        since they are supposed to be filled in by des2chan
;
         ioff = ioffset(i,j)
         ibi  = ibit_chan(i,j)
         mu_parse_channel_descriptor_sa,event_description,nchans,ichans,ichane, $
                  ioff,ibi
         ioffset(i,j)   = ioff
         ibit_chan(i,j) = ibi
      endelse

      nhistos(i,j) = nhisto
      nchannels(i,j) = nchans
;      channels1=channels1(0:nchans-1,*,*)     ; trim array
;      channels2=channels2(0:nchans-1,*,*)     ; trim array
;      ichans=ichans(0:nchans-1)               ; trim array
;      ichane=ichane(0:nchans-1)               ; trim array

      channels1(0:nchans-1,i,j) = ichans(0:nchans-1)
      channels2(0:nchans-1,i,j) = ichane(0:nchans-1)

      fxbclose,unit
   endfor
endfor
end
