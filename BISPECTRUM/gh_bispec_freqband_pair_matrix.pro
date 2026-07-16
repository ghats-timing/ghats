pro gh_bispec_freqband_pair_matrix,f1,f2,bsr,bsi,den1,den2,nprod,band_edges,bicoh_pair,bphase_pair,npix_pair,nprod_pair,ceff,neff,vphase,order_index,exclude_centers=exclude_centers,exclude_widths=exclude_widths,exclude_factor=exclude_factor,help=help

if keyword_set(help) then begin
   print,'GH_BISPEC_FREQBAND_PAIR_MATRIX'
   print,'Build frequency-band pair bispectral matrix from raw 2D sums.'
   return
endif

if n_elements(exclude_factor) eq 0 then exclude_factor=1.0d0

nband=n_elements(band_edges)-1
nf1=n_elements(f1)
nf2=n_elements(f2)

bicoh_pair=fltarr(nband,nband)*!values.f_nan
bphase_pair=fltarr(nband,nband)*!values.f_nan
npix_pair=lonarr(nband,nband)
nprod_pair=dblarr(nband,nband)

bsum_pair=dcomplexarr(nband,nband)
den1_pair=dblarr(nband,nband)
den2_pair=dblarr(nband,nband)

for i=0L,nband-1L do begin
   for j=0L,nband-1L do begin

      f1lo=double(band_edges[i])
      f1hi=double(band_edges[i+1])
      f2lo=double(band_edges[j])
      f2hi=double(band_edges[j+1])

      for ii=0L,nf1-1L do begin
         ff1=double(f1[ii])
         if ff1 lt f1lo or ff1 ge f1hi then continue

         for jj=0L,nf2-1L do begin
            ff2=double(f2[jj])
            if ff2 lt f2lo or ff2 ge f2hi then continue

            if nprod[ii,jj] le 0 then continue
            if ~finite(bsr[ii,jj]) or ~finite(bsi[ii,jj]) then continue
            if ~finite(den1[ii,jj]) or ~finite(den2[ii,jj]) then continue
            if den1[ii,jj] le 0.0d0 or den2[ii,jj] le 0.0d0 then continue

            fsum=ff1+ff2
            excluded=0

            if n_elements(exclude_centers) gt 0 then begin
               for k=0L,n_elements(exclude_centers)-1L do begin
                  ew=0.5d0*exclude_factor*double(exclude_widths[k])
                  ec=double(exclude_centers[k])
                  if abs(ff1-ec) le ew then excluded=1
                  if abs(ff2-ec) le ew then excluded=1
                  if abs(fsum-ec) le ew then excluded=1
               endfor
            endif

            if excluded then continue

            z=dcomplex(double(bsr[ii,jj]),double(bsi[ii,jj]))
            bsum_pair[i,j]=bsum_pair[i,j]+z
            den1_pair[i,j]=den1_pair[i,j]+double(den1[ii,jj])
            den2_pair[i,j]=den2_pair[i,j]+double(den2[ii,jj])
            npix_pair[i,j]=npix_pair[i,j]+1L
            nprod_pair[i,j]=nprod_pair[i,j]+double(nprod[ii,jj])
         endfor
      endfor

      if den1_pair[i,j] gt 0.0d0 and den2_pair[i,j] gt 0.0d0 then begin
         bicoh_pair[i,j]=abs(bsum_pair[i,j])^2/(den1_pair[i,j]*den2_pair[i,j])
         bphase_pair[i,j]=atan(imaginary(bsum_pair[i,j]),float(bsum_pair[i,j]))
      endif

   endfor
endfor

w=where(finite(bicoh_pair) and bicoh_pair gt 0.0,nw)
ceff=!values.f_nan
neff=!values.f_nan
vphase=!values.f_nan
order_index=!values.f_nan

if nw gt 0 then begin
   b=bicoh_pair[w]
   psum=total(double(b))
   if psum gt 0.0d0 then begin
      p=double(b)/psum
      ceff=max(p)
      neff=1.0d0/total(p^2)

      vx=0.0d0
      vy=0.0d0
      for kk=0L,nw-1L do begin
         idx=w[kk]
         ph=bphase_pair[idx]
         vx=vx+double(bicoh_pair[idx])*cos(ph)
         vy=vy+double(bicoh_pair[idx])*sin(ph)
      endfor
      vphase=sqrt(vx^2+vy^2)/psum

      sum_order=0.0d0
      sum_offdiag=0.0d0
      for i=0L,nband-1L do begin
         for j=0L,nband-1L do begin
            if i ne j and finite(bicoh_pair[i,j]) then begin
               sum_offdiag=sum_offdiag+double(bicoh_pair[i,j])
               if i lt j then sum_order=sum_order+double(bicoh_pair[i,j])
            endif
         endfor
      endfor
      if sum_offdiag gt 0.0d0 then order_index=sum_order/sum_offdiag
   endif
endif

print,'GH_BISPEC_FREQBAND_PAIR_MATRIX'
print,'  nband: ',nband
print,'  C_pair max/sum: ',ceff
print,'  N_eff pairs  : ',neff
print,'  V_phase      : ',vphase
print,'  Order index  : ',order_index
print,' '
print,'  band_edges:'
print,band_edges
print,' '
print,'  bicoh_pair:'
print,bicoh_pair
print,' '
print,'  bphase_pair:'
print,bphase_pair
print,' '
print,'  npix_pair:'
print,npix_pair

end
