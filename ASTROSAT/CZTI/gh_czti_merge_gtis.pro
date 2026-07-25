PRO gh_czti_merge_gtis,gti1,gti2,ngti,deltagap,help=help
;+
; NAME:
;      GH_CZTI_MERGE_GTIS
; PURPOSE:
;      Merging of CZTI GTIs shorter than a threshold
; EXPLANATION:
;      The CZTI data come with a large number of microgaps. This procedure allows
;      to ignore the microgaps shorter than deltagap
;
; CALLING SEQUENCE:
;       GH_CZTI_MERGE_GTIS,gti1,gti2,ngti,deltagap
;       H_CZTI_MERGE_GTIS,gti1,gti2,ngti,deltagap
; INPUTS:
;       GTI1        = start times of GTIs
;       GTI2        = end times of GTIs
;       NGTI        = number of GTIs
;       deltagap    = maximum gap to remove
;		
; KEYWORDS:
;		HELP        = Print help and return
;
; OUTPUTS:
;       NONE
;
; EXAMPLE:
;       NONE
;
; COMMON BLOCKS:
;       NONE
; ROUTINES USED:
;       NONE
; NOTES
;       This is a helper routine for CZTI GTI processing, not a standalone
;       analysis command.
; MODIFICATION HISTORY:
;
;   T. Belloni  27 Jun 2017    from scratch
;   2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.
;-
;-------------------------------------------------------------

if(keyword_set(help)) then begin
   print,'GH_CZTI_MERGE_GTIS'
   print,'Helper routine: normally called by CZTI processing code.'
   print,'Calling sequence: gh_czti_merge_gtis, gti1, gti2, ngti, deltagap'
   print,'Merges/handles CZTI GTI gaps shorter than DELTAGAP.'
   return
endif

gaps = gti1(1:ngti-1)-gti2(0:ngti-1)

END
