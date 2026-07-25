;+
; NAME: 
;      GH_LOWCOH_SSQ
; PURPOSE: 
;      Calculate the necesarry input to obtain Coherence error estimates
;      in the High Powers, Low Measured Coherence case of VN97.
; CALLING SEQUENCE:
;      ssq = GH_LOWCOH_SSQ(asq, nsq)
;      ci  = GH_LOWCOH_SSQ_CI(CL, asq, nsq)
; KEYWORDS:
;      HELP    Print help and return a NaN sentinel.
; NOTES:
;      These functions are helper routines used by GH_CROSS_RE_IM_NEW, not
;      standalone analysis commands.
;      Because these are IDL FUNCTIONs, /HELP must be called in function form:
;      dummy = GH_LOWCOH_SSQ(/HELP) or dummy = GH_LOWCOH_SSQ_CI(/HELP).
; MODIFICATION HISTORY: 
;      F. Garcia  22 Mar 2024  from scratch
;      2026 Jul 25  M. Mendez/Codex  Added helper /HELP text.


FUNCTION PROBDIST_CROSS, ssq, asq=asq, nsq=nsq
    ; Probability distribution of SSQ (VN97 eq 6), typo corrected (s*n)^-1.
    a = sqrt(asq)
    n = sqrt(nsq)
    s = sqrt(ssq)
    axn = 2.0*a*s/nsq < 700
    eexp = EXP(-0.5*(asq+2*ssq)/nsq)
    b1 = BESELI(axn, 0)
    maann = 0.5*asq/nsq < 700
    b2 = BESELI(maann, 0)
    ret = eexp*b1/b2/sqrt(!pi)/(n*s)
    return, ret
END

FUNCTION INVERSE_CDF, ssq, asq=asq, nsq=nsq, LVL=LVL
    ; Calculate CDF confidence interval by integrating a PDS.    
    S = 0
    FOR i = 0, 33 DO BEGIN
        TRAPZD,'PROBDIST_CROSS', 1e-6, ssq, S, 1000, asq=asq, nsq=nsq
    ENDFOR
    return, S - LVL
END

FUNCTION gh_lowcoh_ssq_CI,CL,asq,nsq,help=help
    if keyword_set(help) then begin
       print,'GH_LOWCOH_SSQ_CI'
       print,'Helper function: normally called by GH_CROSS_RE_IM_NEW.'
       print,'This is an IDL FUNCTION; use dummy = gh_lowcoh_ssq_CI(/help).'
       print,'Calling sequence: ci = gh_lowcoh_ssq_CI(CL, asq, nsq)'
       print,'Returns confidence limits for SSQ in the VN97 low-coherence case.'
       return, [!values.d_nan,!values.d_nan]
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

FUNCTION gh_lowcoh_ssq,asq,nsq,help=help
     if keyword_set(help) then begin
        print,'GH_LOWCOH_SSQ'
        print,'Helper function: normally called by GH_CROSS_RE_IM_NEW.'
        print,'This is an IDL FUNCTION; use dummy = gh_lowcoh_ssq(/help).'
        print,'Calling sequence: ssq = gh_lowcoh_ssq(asq, nsq)'
        print,'Returns the expected SSQ value for the VN97 low-coherence case.'
        return, !values.d_nan
     endif
     ; Analytic formula for average SSQ (eqn 3.20 Chakrabarty 1995)
     xsq = asq/nsq    
     mxsq = 0.5*xsq < 700
     return,0.5*nsq*(1+xsq+xsq*BESELI(mxsq,1)/BESELI(mxsq,0))  
END
