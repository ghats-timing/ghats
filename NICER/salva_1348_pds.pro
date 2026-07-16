PRO salva_1348_pds,filename

   inputfile = filename+'.pds'
   buoni_1348,buoni,/flares
   ghx,inputfile,nu,pow,pow_e,sel=buoni
   gh_reb,nu,pow,pow_e,-100,x,y,ye
   gh_xspec,x,y,ye,filename
   gh_licu,inputfile,time,licu
   print,mean(licu[buoni])

END