pro gh_read_colors,unit,nft,rmjd,r1,r2,r3

rmjd             = 0.0d0
cnts             = 0.0
poisson          = 0.0
current_vle_rate = 0.0
fndet            = 0.0
a0               = 0.0
r1               = 0.0
r2               = 0.0
r3               = 0.0
pwr              = fltarr(nft)
;
; actual reading
;
readu,unit,rmjd,cnts,fndet,a0,poisspn,current_vle_rate, $
            r1,r2,r3,pwr

end
