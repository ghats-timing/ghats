pro gh_bispec_target_pair_matrix,f1,f2,bsr,bsi,den1,den2,nprod,band_edges,target_c,target_w,bicoh_pair,bphase_pair,npix_pair,nprod_pair,btot,ceff,neff,vphase,help=help

if keyword_set(help) then begin
   print,'GH_BISPEC_TARGET_PAIR_MATRIX'
   print,'Low-frequency band-pair matrix, restricted by f1+f2 near target frequencies.'
   return
endif

nband=n_elements(band_edges)-1
ntarg=n_elements(target_c)
nf1=n_elements(f1)
nf2=n_elements(f2)

bicoh_pair=fltarr(nband,nband,ntarg)*!values.f_nan
bphase_pair=fltarr(nband,nband,ntarg)*!values.f_nan
npix_pair=lonarr(nband,nband,ntarg)
nprod_pair=dblarr(nband,nband,ntarg)

bsum_pair=dcomplexarr(nband,nband,ntarg)
den1_pair=dblarr(nband,nband,ntarg)
den2_pair=dblarr(nband,nband,ntarg)

btot=dblarr(ntarg)*!values.f_nan
ceff=dblarr(ntarg)*!values.f_nan
neff=dblarr(ntarg)*!values.f_nan
vphase=dblarr(ntarg)*!values.f_nan

for it=0L,ntarg-1L do begin
   tc=double(target_c[it])
   tw=0.5d0*double(target_w[it])

   for i=0L,nband-1L do begin
      f1lo=double(band_edges[i])
      f1hi=double(band_edges[i+1])

      for j=0L,nband-1L do begin
         f2lo=double(band_edges[j])
         f2hi=double(band_edges[j+1])

         for ii=0L,nf1-1L do begin
            ff1=double(f1[ii])
            if ff1 lt f1lo or ff1 ge f1hi then continue

            for jj=0L,nf2-1L do begin
               ff2=double(f2[jj])
               if ff2 lt f2lo or ff2 ge f2hi then continue

               fsum=ff1+ff2
               if abs(fsum-tc) gt tw then continue

               if nprod[ii,jj] le 0 then continue
               if ~finite(bsr[ii,jj]) or ~finite(bsi[ii,jj]) then continue
               if ~finite(den1[ii,jj]) or ~finite(den2[ii,jj]) then continue
               if den1[ii,jj] le 0.0d0 or den2[ii,jj] le 0.0d0 then continue

               z=dcomplex(double(bsr[ii,jj]),double(bsi[ii,jj]))
               bsum_pair[i,j,it]=bsum_pair[i,j,it]+z
               den1_pair[i,j,it]=den1_pair[i,j,it]+double(den1[ii,jj])
               den2_pair[i,j,it]=den2_pair[i,j,it]+double(den2[ii,jj])
               npix_pair[i,j,it]=npix_pair[i,j,it]+1L
               nprod_pair[i,j,it]=nprod_pair[i,j,it]+double(nprod[ii,jj])
            endfor
         endfor

         if den1_pair[i,j,it] gt 0.0d0 and den2_pair[i,j,it] gt 0.0d0 then begin
            bicoh_pair[i,j,it]=abs(bsum_pair[i,j,it])^2/(den1_pair[i,j,it]*den2_pair[i,j,it])
            bphase_pair[i,j,it]=atan(imaginary(bsum_pair[i,j,it]),float(bsum_pair[i,j,it]))
         endif

      endfor
   endfor

   tmp=bicoh_pair[*,*,it]
   w=where(finite(tmp) and tmp gt 0.0,nw)

   if nw gt 0 then begin
      b=double(tmp[w])
      psum=total(b)
      btot[it]=psum

      if psum gt 0.0d0 then begin
         p=b/psum
         ceff[it]=max(p)
         neff[it]=1.0d0/total(p^2)

         ph=bphase_pair[*,*,it]
         vx=0.0d0
         vy=0.0d0
         for kk=0L,nw-1L do begin
            vx=vx+b[kk]*cos(ph[w[kk]])
            vy=vy+b[kk]*sin(ph[w[kk]])
         endfor
         vphase[it]=sqrt(vx^2+vy^2)/psum
      endif
   endif

endfor

print,'GH_BISPEC_TARGET_PAIR_MATRIX'
print,'  nband: ',nband
print,'  ntarg: ',ntarg
print,'  band_edges:'
print,band_edges

for it=0L,ntarg-1L do begin
   print,' '
   print,'TARGET ',it,'  center=',target_c[it],'  width=',target_w[it]
   print,'  Btot  = ',btot[it]
   print,'  Ceff  = ',ceff[it]
   print,'  Neff  = ',neff[it]
   print,'  Vphase= ',vphase[it]
   print,'  bicoh_pair:'
   print,bicoh_pair[*,*,it]
   print,'  bphase_pair:'
   print,bphase_pair[*,*,it]
   print,'  npix_pair:'
   print,npix_pair[*,*,it]
endfor

end
