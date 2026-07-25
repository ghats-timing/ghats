;+
; NAME: 
;      GH_LOWCOH_SSQ_CI
; PURPOSE: 
;      Calculate the necesarry input to obtain Coherence error estimates
;      in the High Powers, Low Measured Coherence case of VN97.
; CALLING SEQUENCE:
;      GH_LOWCOH_SSQ_CI, CL, asq, nsq
; KEYWORDS:
;      HELP    Print help and return.
; NOTES:
;      This legacy procedure is not the callable helper used by
;      GH_CROSS_RE_IM_NEW. GH_CROSS_RE_IM_NEW uses the FUNCTION
;      GH_LOWCOH_SSQ_CI defined in gh_lowcoh_ssq.pro. On case-sensitive
;      systems, this mixed-case filename may also not autocompile from a
;      lower-case IDL command.
; MODIFICATION HISTORY: 
;      F. Garcia  22 Mar 2024  from scratch
;      2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.


pro gh_lowcoh_ssq_CI,CL,asq,nsq,help=help
    if keyword_set(help) then begin
       print,'GH_LOWCOH_SSQ_CI'
       print,'Legacy procedure. GH_CROSS_RE_IM_NEW uses the function in gh_lowcoh_ssq.pro.'
       print,'For callable function help, compile gh_lowcoh_ssq.pro and use:'
       print,'  dummy = gh_lowcoh_ssq_CI(/help)'
       print,'Calling sequence: gh_lowcoh_ssq_CI, CL, asq, nsq'
       print,'Computes confidence limits for SSQ in the VN97 low-coherence case.'
       return
    endif
    ; Calculate the Confidence Limits of SSQ (VN97 eqn 9)
    if (asq/nsq GT 300) then begin
    ; if S/N is large, use the Gaussian approximation (VN97 eqn 7)
       lower_bound = asq-sqrt(2*nsq*asq)
       upper_bound = asq+sqrt(2*nsq*asq)
    endif else begin
    ; Confidence levels as roots of the CDF (VN97 eqn 6)
        lower_bound = ZBRENT(1e-6,asq+100*nsq,FUNC='INVERSE_CDF', asq=asq, nsq=nsq, LVL=0.5*(1-CL))
        upper_bound = ZBRENT(1e-6,asq+100*nsq,FUNC='INVERSE_CDF', asq=asq, nsq=nsq, LVL=0.5*(1+CL))
    endelse
    return, [lower_bound, upper_bound]
END
