function mu_bit_extractor,inte,nume
;
; Check whether a certain bit of a short integer variable is set
;  or not. Input bit number runs from 1 to 16.
;
;      Author: TMB       Date: 14-APR-2009   Vectorialized
;
mask = [  1,  2,  4,   8 ,  16,  32,  64,  128,   $
        256,512,1024,2048,4096,8192,16384,-32767-1 ]

inume=17-nume
aa = inte or mask(inume-1)
bb = (aa eq inte)
return,bb
;if(aa eq inte) then begin
;   return,1
; endif
;
;return,0
end
