PRO laxpc_level1_split,dirname
;
; Input filename: gti file
;
;gtifile=file_search('*.gti') 
gtifile = 'usergti.fits'
gti=mrdfits(gtifile,1)

gti1 = gti.start
gti2 = gti.stop

; divide into intervals of 2ks exposure

diff = gti2-gti1
n    = n_elements(diff)
inizio = gti1[0]
fine   = gti2[n-1]

cumsum = total(diff,/cumulative)

slice  = 2000
quanti = fix(cumsum[n-1]/slice)+1

eventfile1=file_search(dirname+'/lxp1/modeEA/*fits') 
eventfile2=file_search(dirname+'/lxp2/modeEA/*fits') 
eventfile3=file_search(dirname+'/lxp3/modeEA/*fits') 

t1=gti1[0]

xxx
FOR i=1,quanti DO BEGIN
    outfilename1 = dirname+'/lxp1/modeEA/pippo'+strtrim(string(i),1)
	outfilename2 = dirname+'/lxp2/modeEA/pippo'+strtrim(string(i),1)
	outfilename3 = dirname+'/lxp3/modeEA/pippo'+strtrim(string(i),1)
	limite  = i*slice
	iii     = max(where(cumsum le limite))
	t2      = gti2[iii]
    stringa = 'fselect '+eventfile1+' '+outfilename1+' "(TIME > '+string(t1)+') && (TIME <= '+string(t2)+')"' 
	print,stringa
    stringa = 'fselect '+eventfile2+' '+outfilename2+' "(TIME > '+string(t1)+') && (TIME <= '+string(t2)+')"' 
	print,stringa
    stringa = 'fselect '+eventfile3+' '+outfilename3+' "(TIME > '+string(t1)+') && (TIME <= '+string(t2)+')"' 
	print,stringa
	t1      = t2
ENDFOR

END