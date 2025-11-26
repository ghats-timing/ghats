;+
; NAME: 
;      GH_LOWCOH_SSQ
; PURPOSE: 
;      Calculate the necesarry input to obtain Coherence error estimates
;      in the High Powers, Low Measured Coherence case of VN97.
; NOTES:
;      This functions are used by GH_CROSS_RE_IM_NEW.
; MODIFICATION HISTORY: 
;      F. Garcia  22 Mar 2024  from scratch


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

FUNCTION gh_lowcoh_ssq_CI,CL,asq,nsq
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

FUNCTION gh_lowcoh_ssq,asq,nsq
     ; Analytic formula for average SSQ (eqn 3.20 Chakrabarty 1995)
     xsq = asq/nsq    
     mxsq = 0.5*xsq < 700
     return,0.5*nsq*(1+xsq+xsq*BESELI(mxsq,1)/BESELI(mxsq,0))  
END
