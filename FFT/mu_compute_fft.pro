pro mu_compute_fft,rdata,pwr,np,T,start_time,source,observatory,           $
           instrument,out_file,chout,ndet,i_vle,                   $
           current_vle_rate,tstart_fft,                       $
       new_output_flag,output_unit,ntotal_ffts,canali,        $
	   std21,std22,std23, $
	   BANDE,turbo,gti_flag,proliferation     ; MU6
;
; Perform and write FFT/PDS to disk
; The heart of MUFFT
;
;-----------------------------------------------------------------------
; Parameters
;
; new_subfile       I/O: flag for first call (needed by write_*)
;-----------------------------------------------------------------------

common finestre,finestra,winn

EPS         =  1.0e-6      ; precision tolerance

if(chout eq 'POWER') then begin
   out_string='PowerSpectrum'
  endif else begin
   if(chout eq 'FOURIER') then out_string='FourierSpectrum'
endelse

time_now    = start_time
nsel_chan   = 1
;
; Prepare info string to write to output header
;
osserv     ='                '
sorgente   ='                '
strumento  ='                '
strput,osserv,observatory,0
strput,sorgente,source,0
strput,strumento,instrument,0
;
cnts = total(rdata)
;
; computation of Poissonian level
;
if(observatory eq 'XTE') then begin
   mu_zhang,cnts,T,ndet,np,current_vle_rate,i_vle,poi
endif else begin
    poi = 2.0
endelse
;
; Data windowing
rdata = rdata * finestra
cnts = total(rdata)    ; recalculation
;
;  Do the FFT
;  Here we assume that the normalization of the IDL FFT is the same
;  as that of Ramach's routine
;
rdata=fft(rdata,1,/OVERWRITE)   ; overwritten on rdata (now complex)
;
;  This is a call to the INVERSE FFT from IDL, which according to
;  Ramachandran corresponds to the FORWARD transform in his routine.
;  The normalization should be the same (fingers crossed)
;
if(cnts ge EPS) then begin
   ; non-zero total counts
   ; computer power spectrum
   a0=abs(rdata(0))*2.0   ; a0 should be real anyway
   pwr = abs(rdata(1:np/2))^2 *2.0/cnts ; compute power spectrum (Leahy)
  endif else begin
   a0 = 0.0
   pwr = pwr*0.0
endelse

time_now = start_time
current_channel = 1
chou2=strlowcase(strmid(chout,0,1))

case chou2 of
   'p':     begin
              mu_write_single_power,a0,pwr,np,poi,cnts,T,start_time,out_file,   $
                         ndet,             $
              current_vle_rate,new_output_flag,             $
              output_unit,ntotal_ffts,canali,osserv,sorgente, $
              strumento,std21,std22,std23,   $
              BANDE,turbo,gti_flag,proliferation,i_vle     ; MU6
            end
   'f':     begin
                 fftdata = rdata(1:np/2)
                 mu_write_single_fft,fftdata,np,poi,cnts,T,start_time,out_file,   $
                       ndet,             $
              current_vle_rate,new_output_flag,             $
              output_unit,ntotal_ffts,canali,osserv,sorgente, $
              strumento,std21,std22,std23,   $
              BANDE,turbo,gti_flag,proliferation,i_vle     ; MU6
              end

   else: begin
           print,chout
            massage,'I do not know what to write this: please choose either POWER or FFT'
            retall
     end
endcase
;
; Register that first call to write_* has been done
;
new_subfile = 0
;
; Rdata redefined to get back to real
;
rdata = fltarr(np)*0.0

end
