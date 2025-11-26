pro mu_get_housekeeping,filename,tstart,tstop,housedata, $
                 vle_times,vle_counts,npcus,MAX_VLE,nvle,  $
		eband,rate1,rate2,rate3
;
; Procedure to get housekeeping information out of RXTE Standard2
; data. Mainly important for VLE rates, used for the Poissonian estimate.
; It uses only the Standard2 file that covers the mid-point of the data file
; This can fail when there are two Standard2 files for a given data file
;
; Mu 2.2 system-independent version
; Mu 3.0 std2 band calls
; Mu 6.0 vle_times zeroed at each call
;
;-----------------------------------------------------------------------------
; Parameters
; filename              I: input data filename
; tstart                I: start time for data file
; tstop                 O: end time for data file
; housedata             O: flag for housekeeping file found
; vle_times             O: array with times for VLE rate
; vle_counts            O: array with VLE rate
; npcus                 O: array with number of PCUs on
; MAX_VLE               I: maximum number of VLE rows allowed
; nvle                  O: actual number of VLE data lines
; eband                 I: array with energy bands for light curve
;                               production (e.g. [0,13,14,35,36,101])
; rate1                 O: array with rate in band 1
; rate2                 O: array with rate in band 2
; rate3                 O: array with rate in band 3
;-----------------------------------------------------------------------------
;
; Major constants
;
MAX_TYPES    =   5
C_TYPES      =   6
MAX_COL      = 256
MAX_FILES    =  40
;
; Prepare keyword names for VLE rates (per PCU)
;
datatypes    = strarr(MAX_TYPES)
datatypes(0) = 'VLECntPcu0'
datatypes(1) = 'VLECntPcu1'
datatypes(2) = 'VLECntPcu2'
datatypes(3) = 'VLECntPcu3'
datatypes(4) = 'VLECntPcu4'
;
; Prepare keyword names for PCU2  columns
;
datatypesc   = strarr(C_TYPES)
datatypesc(0) = 'X1LSpecPcu2'
datatypesc(1) = 'X1RSpecPcu2'
datatypesc(2) = 'X2LSpecPcu2'
datatypesc(3) = 'X2RSpecPcu2'
datatypesc(4) = 'X3LSpecPcu2'
datatypesc(5) = 'X3RSpecPcu2'
;
; Array to contain column numbers corresponding to VLE rates
;
vlecol       = intarr(MAX_TYPES)
;
; Array to contain column numbers corresponding to PCU2 rates
;
ccol         = intarr(C_TYPES)
;
; Look for FS4a filenames
;
; New call for IDL 6
dirname=file_dirname(filename)
; GDL compatible call 19-Nov-2008 TMB
;dum = fsc_base_filename(filename,Directory=dirname)
;islash = 0
;islash = strpos(filename,'/FS',/REVERSE_SEARCH)
;if(islash ge 0) then begin
;   dirname = strtrim(strmid(filename,0,islash),2)
;  endif ELSE begin
;   dirname = '.'
;endelse

;
; System independent call for MU 3.0
;
search_string = dirname+'/FS4a*'
filelist=file_search(search_string)
po              = strpos(filelist,'_bkg')
filelist        = filelist(where(po eq -1)) ; remove bkg files

nfiles=n_elements(filelist)
housefile=''
housedata=1
if(nfiles  le 0 or (strlen(filelist(0)) eq 0)) then begin
   print,'No Standard2 data found in directory ',dirname
   housedata=0
   return
endif
;
; Find the standard2 file that corresponds to the current data file times
;
mu_find_hk_file,housefile,tstart,tstop,filelist,nfiles

if(housefile eq '') then begin
   housedata=0
   print,'No Standard2 file found that covers midpoint of data file ',$
         strtrim(filename,2)
   return
endif
print,'Using Standard2 file ',strtrim(housefile,2),' for VLE rates'
;
; Now access the standard2 file
;
unit = 80
;
; Open the file
;
fxbopen,unit,housefile,1,header,errmsg=errmsg
;
; Read in basic header information
;
t1        = fxpar(header,'TSTART')
t2        = fxpar(header,'TSTOP')
timedel   = fxpar(header,'TIMEDEL')
datamode  = strtrim(fxpar(header,'DATAMODE'))
tdim      = fxpar(header,'TDIM*')
nfields   = strtrim(fxpar(header,'TFIELDS'))
ttype     = fxpar(header,'TTYPE*')
nrows     = fxpar(header,'NAXIS*')
idltype   = intarr(nfields)
;
; Identifies which column numbers correspond to the VLE and Time columns
;
timecol = -100  ; marker for column number corresponding to Time
fxbtform,header,tbcol,idltype,formato,numval,maxval
ttype       = fxpar(header,'TTYPE*')
for i=0,nfields-1 do begin
   for j=0,4 do begin
      if(ttype(i) eq datatypes(j)) then vlecol(j) = i+1
   endfor
   if(ttype(i) eq 'Time    ') then timecol = i+1
endfor
;
; Identifies which column numbers correspond to the PCU2 columns
;
for i=0,nfields-1 do begin
   for j=0,5 do begin
      if(ttype(i) eq datatypesc(j)) then ccol(j) = i+1
   endfor
endfor

nvle = nrows(1)   ; number of rows is second element, first is n. of bytes
if(nvle gt MAX_VLE) then begin
   massage,'Increase array size for vle counts'
   retall
endif
;
; Now loop over the STD2 rows
;
vle_times = vle_times * 0  ; nulling for new use: otherwise problems in mu_read_std2_rates
for irow=0,nvle-1 do begin
   npcus(irow)      = 0
   vle_counts(irow) = 0.0
;
; For each row, read in time and vle rates
; The assumption is that if the VLE counts for a PCU are > 1, the PCU is on
   fxbread,unit,timel,timecol,irow+1
   vle_times(irow) = timel
   for j=0,4 do begin
      fxbread,unit,icounts,vlecol(j),irow+1
      if (icounts gt 1) then begin
         npcus(irow) = npcus(irow)+1
      endif
      vle_counts(irow) = vle_counts(irow) + icounts
   endfor
;
;  Here accumulate the three PCU2 light curves
;
   spettro = intarr(129)
   for j=0,5 do begin
;       Add up all layers, left and right
      fxbread,unit,spe,ccol(j),irow+1
      spettro =  spettro + spe
   endfor
   rate1(irow) = float(total(spettro(eband(0):eband(1))))
   rate2(irow) = float(total(spettro(eband(2):eband(3))))
   rate3(irow) = float(total(spettro(eband(4):eband(5))))
;
;  Convert VLE counts into counts per second (VLE rate)
;
   vle_counts(irow) = vle_counts(irow)/timedel
endfor
;
; Close standard2 file
;
fxbclose,unit

end
