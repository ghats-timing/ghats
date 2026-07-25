pro gh_read_colors,unit,nft,rmjd,r1,r2,r3,help=help
;+
; NAME:
;      GH_READ_COLORS
;
; PURPOSE:
;      Read one GHATS colour/light-curve record from an opened PDS file.
;
; CALLING SEQUENCE:
;      GH_READ_COLORS, unit, nft, rmjd, r1, r2, r3
;      GH_READ_COLORS, /HELP
;
; INPUTS:
;      UNIT    Open file unit positioned at a colour/PDS record.
;      NFT     Number of Fourier/time bins used to size the ignored PDS array.
;
; OUTPUTS:
;      RMJD    Record time.
;      R1      First colour/rate channel.
;      R2      Second colour/rate channel.
;      R3      Third colour/rate channel.
;
; KEYWORDS:
;      HELP    Print help and return.
;
; NOTES:
;      This is a helper routine normally called by GH_COLORS, not a standalone
;      analysis command.
;
; MODIFICATION HISTORY:
;      2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.
;-

if(keyword_set(help)) then begin
   print,'GH_READ_COLORS'
   print,'Helper routine: normally called by GH_COLORS, not run standalone.'
   print,'Calling sequence: gh_read_colors, unit, nft, rmjd, r1, r2, r3'
   print,'Reads one colour/light-curve record from an opened GHATS file.'
   return
endif

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
