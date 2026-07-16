pro mu_read_czti,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,quadrante,rsour
;
;  Read in CZTI photon files
;  jmeta is the metafile number but ALSO the quadrant number!
;
;-------------------------------------------------------------------------
; Parameters
;
; unit                 I: unit for file reading (already opened)
; tag                  I: file tag identifier
; rdata                O: array with extracted light curve
; np                   I: number of points to accumulate
; channels             I: array of selected/non-selected channels
; selected             I: still a mystery to me
; jmeta                I: index of current metafile
; tres                 I: time resolution of current file
; tres_fft             I: time resolution needed for light curve
; tstart_fft           I: fft start time
; tend_fft             I: fft end time
; irowl                I: current row
;-------------------------------------------------------------------------
common dati,tempi0,tempi1,tempi2,tempi3,energy0,energy1,energy2,energy3,peso0,peso1,peso2,peso3, $
	   canale_corrente0,canale_corrente1,canale_corrente2,canale_corrente3, $
		   std21,std22,std23,  $
			   nrows0,nrows1,nrows2,nrows3   ; common block for keeping the datacommon barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)
common barycentered,baryflag               ; flag for barycentered photons (11-Oct-2009)

CASE quadrante of

	0: mu_read_czti0,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,rsour

	1: mu_read_czti1,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,rsour

	2: mu_read_czti2,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,rsour

	3: mu_read_czti3,unit,tag,rdata,np,channels,selected,jmeta,tres,tres_fft, $
                  tstart_fft,tend_fft,irowl, $
		          channels_std21,channels_std22,channels_std23,nffts,energy,fraction,rsour
ENDCASE
end
