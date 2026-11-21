*! version 0.2.1  16oct2026
*! svytmle: survey-weighted, PSU-level cross-fitted TMLE (or AIPW) of the population
*! ATE with a design-based (Taylor linearization) variance.
*! Stata implementation of the R package svytmle (Karim 2026, arXiv:2606.30918).
*! Copyright (C) 2026 M. Ehsan Karim. License: GPL-3.0-or-later (LICENSE.md in the repository below).
*! https://github.com/ehsanx/svytmle-stata

program define svytmle, eclass properties(mi)
    version 15
    if replay() {
        if "`e(cmd)'" != "svytmle" error 301
        Display `0'
        exit
    }
    Estimate `0'
end

program define Estimate, eclass
    version 15
    syntax varname(numeric) [if] [in] , Treatment(varname numeric)      ///
        COVariates(varlist numeric fv) [ SUBpop(string asis)            ///
        Folds(integer 5) FOLDVar(varname numeric) SEED(string)          ///
        GBound(real 0.05) NUISance(string) ESTimator(string)            ///
        STRICTfolds Level(cilevel) RR OR GENerate(name) REPLACE         ///
        noHEADer * ]
    _get_diopts diopts, `options'

    local y `varlist'
    local a `treatment'
    local cmdline `"svytmle `0'"'

    * ---------------------------------------------------------------- options
    if `"`if'`in'"' != "" {
        di as err "{p 0 4 2}{bf:if} and {bf:in} are not allowed: select a "  ///
            "survey domain with {bf:subpop()}, so that the variance is "      ///
            "computed on the full design (subset the design, not the data){p_end}"
        exit 101
    }
    if "`y'" == "`a'" {
        di as err "the outcome and {bf:treatment()} must be different variables"
        exit 198
    }
    if `folds' < 2 {
        di as err "{bf:folds()} must be at least 2"
        exit 198
    }
    if `gbound' <= 0 | `gbound' >= 0.5 {
        di as err "{bf:gbound()} must be in (0, 0.5)"
        exit 198
    }
    if inlist(`"`nuisance'"', "", "unweighted") local nuisw 0
    else if `"`nuisance'"' == "survey" local nuisw 1
    else {
        di as err "{bf:nuisance()} must be {bf:unweighted} or {bf:survey}"
        exit 198
    }
    local nuisance = cond(`nuisw', "survey", "unweighted")
    if inlist(`"`estimator'"', "", "tmle") local est 0
    else if `"`estimator'"' == "aipw" local est 1
    else {
        di as err "{bf:estimator()} must be {bf:tmle} or {bf:aipw}"
        exit 198
    }
    local estimator = cond(`est', "aipw", "tmle")
    local strict = ("`strictfolds'" != "")
    quietly fvrevar `covariates', list
    local covbase `r(varlist)'
    if `: list y in covbase' | `: list a in covbase' {
        di as err "{bf:covariates()} must not contain the outcome or the treatment"
        exit 198
    }
    * generate(): validate the names now, create the variables only after a
    * successful fit (estimates go to temporary variables that are renamed)
    local gennames
    local gentmp
    if "`generate'" != "" {
        local gsuf fold g Q0 Q1 Q0star Q1star eif
        if `est' local gsuf fold g Q0 Q1 eif
        foreach s of local gsuf {
            local v `generate'`s'
            capture confirm variable `v', exact
            if !_rc & "`replace'" == "" {
                di as err "variable {bf:`v'} already exists; specify {bf:replace} to overwrite it"
                exit 110
            }
            if _rc confirm new variable `v'
            tempvar t
            local gennames `gennames' `v'
            local gentmp `gentmp' `t'
        }
    }
    else if "`replace'" != "" {
        di as err "{bf:replace} requires {bf:generate()}"
        exit 198
    }

    * ------------------------------------------------------- survey design
    quietly svyset
    if inlist(`"`r(settings)'"', "", ", clear") {
        di as err "{p 0 4 2}data not {bf:svyset}; declare the design first, e.g. "  ///
            "{bf:svyset psu [pweight=wt], strata(stratum)}{p_end}"
        exit 119
    }
    local wtype `"`r(wtype)'"'
    if `"`wtype'"' == "" {
        di as err "{p 0 4 2}svytmle requires a sampling weight: svyset the design "  ///
            "with {bf:[pweight=}{it:varname}{bf:]}{p_end}"
        exit 101
    }
    if "`wtype'" != "pweight" {
        di as err "the survey design must use a {bf:pweight} (found `wtype')"
        exit 101
    }
    local wvar   `r(wvar)'
    local psu    `r(su1)'
    local strat  `r(strata1)'
    local fpc    `r(fpc1)'
    local vce    `r(vce)'
    local stages = r(stages)
    local sdof   = r(dof)
    if `"`r(poststrata)'`r(postweight)'`r(calmethod)'`r(rake)'`r(regress)'"' != "" {
        di as err "{p 0 4 2}post-stratified or calibrated designs are outside the "  ///
            "supported scope of the estimator; svyset the uncalibrated design{p_end}"
        exit 459
    }
    if "`psu'" == "_n" local psu
    if "`vce'" != "" & "`vce'" != "linearized" {
        if "`psu'" == "" & "`strat'" == "" {
            di as err "{p 0 4 2}replicate-weight designs without first-stage PSUs/strata "  ///
                "are not supported; svyset the PSU and strata variables{p_end}"
            exit 459
        }
        di as txt "{p 0 6 2}Note: svytmle uses the linearization variance on the "  ///
            "first-stage PSUs and strata; the svyset vce(`vce') is ignored.{p_end}"
    }
    if "`fpc'" != "" {
        di as txt "{p 0 6 2}Note: the svyset FPC (`fpc') is ignored: the estimator uses "  ///
            "the with-replacement first-stage (no-FPC) variance, which is calibrated "    ///
            "for the super-population ATE.{p_end}"
    }
    if `stages' > 1 & `stages' < . {
        di as txt "{p 0 6 2}Note: only the first sampling stage is used for the variance.{p_end}"
    }
    if `sdof' < . {
        di as txt "{p 0 6 2}Note: the svyset dof(`sdof') is ignored; svytmle computes "  ///
            "the design degrees of freedom from the PSUs and strata.{p_end}"
    }

    tempvar D S H P
    quietly gen byte `D' = !missing(`wvar')
    if "`psu'"   != "" quietly replace `D' = 0 if missing(`psu')
    if "`strat'" != "" quietly replace `D' = 0 if missing(`strat')
    quietly count if `D' & `wvar' < 0
    if r(N) {
        di as err "the sampling weight `wvar' has negative values"
        exit 402
    }

    * ------------------------------------------------------------- domain
    quietly gen byte `S' = 1
    local subpopdesc
    if `"`subpop'"' != "" {
        local 0 `subpop'
        capture noisily syntax [varname(default=none numeric)] [if/]
        if _rc {
            di as err "(in {bf:subpop()}: specify {it:varname} and/or {bf:if} {it:exp})"
            exit _rc
        }
        if "`varlist'" != "" {
            quietly count if `D' & missing(`varlist')
            if r(N) {
                di as err "{p 0 4 2}the subpop() variable `varlist' has " r(N)        ///
                    " missing value(s) in the design; code domain membership 0/1 "    ///
                    "(svytmle does not guess whether those observations are in the "  ///
                    "domain){p_end}"
                exit 416
            }
            quietly replace `S' = 0 if `varlist' == 0
        }
        if `"`if'"' != "" quietly replace `S' = 0 if !(`if')
        local subpopdesc `subpop'
    }
    quietly count if `S' & !`D'
    if r(N) {
        di as err "{p 0 4 2}" r(N) " observation(s) in the analytic domain have a missing "  ///
            "weight, PSU, or stratum; svytmle does not drop observations within a "          ///
            "domain. Exclude them explicitly with {bf:subpop()}.{p_end}"
        exit 416
    }
    quietly count if `D'
    local N = r(N)
    if `N' == 0 {
        di as err "no observations with nonmissing design variables"
        exit 2000
    }
    if `N' < _N {
        di as txt "{p 0 6 2}Note: " _N - `N' " observation(s) with missing "  ///
            "design variables are outside the design.{p_end}"
    }
    if "`strat'" == "" quietly gen long `H' = 1 if `D'
    else quietly egen long `H' = group(`strat') if `D'
    if "`psu'" == "" quietly gen long `P' = sum(`D') if `D'
    else quietly egen long `P' = group(`H' `psu') if `D'

    quietly count if `S'
    local Nsub = r(N)
    if `Nsub' == 0 {
        di as err "the subpopulation contains no observations"
        exit 2000
    }
    foreach v in y a {
        local vn ``v''
        quietly count if `S' & missing(`vn')
        if r(N) {
            di as err "{p 0 4 2}`vn' has " r(N) " missing value(s) in the analytic domain; "  ///
                "svytmle does not impute or drop observations within a domain{p_end}"
            exit 416
        }
        quietly count if `S' & !inlist(`vn', 0, 1)
        if r(N) {
            di as err "`vn' must be coded 0/1 in the analytic domain"
            exit 450
        }
        quietly summarize `vn' if `S', meanonly
        if r(min) == r(max) {
            di as err "`vn' is constant in the analytic domain"
            exit 459
        }
    }
    foreach v of local covbase {
        quietly count if `S' & missing(`v')
        if r(N) {
            di as err "{p 0 4 2}covariate `v' has " r(N) " missing value(s) in the analytic "  ///
                "domain; svytmle does not impute or drop observations within a domain{p_end}"
            exit 416
        }
    }
    quietly summarize `wvar' if `S', meanonly
    if r(sum) <= 0 {
        di as err "all sampling weights are zero in the analytic domain"
        exit 2000
    }
    forvalues t = 0/1 {
        quietly summarize `wvar' if `S' & `a' == `t', meanonly
        if r(sum) <= 0 {
            di as err "{p 0 4 2}the `a' = `t' group has zero total sampling weight in the "  ///
                "analytic domain{p_end}"
            exit 459
        }
    }
    if "`foldvar'" != "" {
        quietly count if `S' & (missing(`foldvar') | `foldvar' != int(`foldvar') | `foldvar' < 1)
        if r(N) {
            di as err "{bf:foldvar()} must hold positive integers in the analytic domain"
            exit 459
        }
    }

    * ------------------------------------------- covariate design matrix
    * fvexpand on the domain (base = lowest level present, as R's droplevels).
    * The FULL marked term list goes to fvrevar: listing levels without their
    * base (2.race 3.race) makes Stata take the lowest LISTED level as the base
    * and silently zero 2.race. Base/omitted terms are then left out by position.
    * Collinearity is NOT screened here: as in R's glm, aliased columns are
    * resolved within each training-fold fit, in column order.
    quietly fvexpand `covariates' if `S'
    local xfull `r(varlist)'
    local xfull : list uniq xfull
    quietly fvrevar `xfull'
    local allv `r(varlist)'
    local xterms
    local xvars
    local k 0
    foreach t of local xfull {
        local ++k
        _ms_parse_parts `t'
        if !r(omit) {
            local xterms `xterms' `t'
            local xvars `xvars' `: word `k' of `allv''
        }
    }
    if "`xterms'" == "" {
        di as err "the covariates have no usable columns in the analytic domain"
        exit 459
    }

    * ---------------------------------------------------------- estimation
    if `"`seed'"' != "" set seed `seed'
    local rngstate `"`c(rngstate)'"'
    local rngname `c(rng_current)'
    tempname b V FP
    mata: _svytmle_main("`y'", "`a'", "`xvars'", "`H'", "`P'", "`wvar'",   ///
        "`D'", "`S'", "`foldvar'", `folds', `gbound', `nuisw', `est',        ///
        `strict', "`b'", "`V'", "`FP'", "`gentmp'")

    if `n_single' > 0 {
        di as err "{p 0 4 2}`n_single' stratum/strata contain a single PSU in the design; "  ///
            "the design-based variance needs >= 2 PSUs per stratum. Collapse such "          ///
            "strata before svyset (see {help svydescribe}).{p_end}"
        exit 460
    }
    if `n_noncon' > 0 {
        di as err "the covariates are constant in the analytic domain"
        exit 459
    }
    if `n_notrain' > 0 {
        di as txt "{p 0 6 2}Warning: in `n_notrain' fold-by-stratum case(s) a stratum had no "  ///
            "PSU in the training fold (it has a single PSU in the analytic domain); its "        ///
            "held-out rows are predicted from a fit that excludes the stratum. Consider "        ///
            "collapsing strata.{p_end}"
    }
    if `n_notrain_fv' > 0 {
        di as txt "{p 0 6 2}Warning: in `n_notrain_fv' fold-by-stratum case(s) foldvar() put "  ///
            "all of a stratum's domain PSUs in the held-out fold, so its held-out rows are "     ///
            "predicted from a fit that excludes the stratum. Assign folds to PSUs within "       ///
            "strata.{p_end}"
    }
    if `rr_clamped' {
        di as txt "{p 0 6 2}Note: a potential-outcome mean lies outside (0, 1) (or within "     ///
            "1e-6 of a bound), so the risk and odds ratios are not reported (stored as "         ///
            "missing).{p_end}"
    }
    if `n_nonconv' > 0 {
        di as txt "{p 0 6 2}Warning: `n_nonconv' logistic fit(s) did not converge in 25 "  ///
            "IRLS iterations (as R's glm.fit; possible separation).{p_end}"
    }
    if `n_dropfits' > 0 {
        di as txt "{p 0 6 2}Note: in `n_dropfits' nuisance fit(s), covariate columns that were "  ///
            "constant or collinear in the training fold were dropped (aliased, as R's glm "       ///
            "does).{p_end}"
    }
    if `df' < 1 {
        di as txt "{p 0 6 2}Warning: the design df are `df'; inference uses t with 1 df "  ///
            "(as the R package).{p_end}"
    }

    * ------------------------------------------- generated variables
    if "`generate'" != "" {
        local labs `" "cross-fitting fold" "out-of-fold propensity score (truncated)" "out-of-fold Q(0,W) (truncated)" "out-of-fold Q(1,W) (truncated)" "targeted Q*(0,W)" "targeted Q*(1,W)" "influence function of the ATE" "'
        if `est' local labs `" "cross-fitting fold" "out-of-fold propensity score (truncated)" "out-of-fold Q(0,W) (truncated)" "out-of-fold Q(1,W) (truncated)" "influence function of the ATE (AIPW)" "'
        local k 0
        foreach v of local gennames {
            local ++k
            local t : word `k' of `gentmp'
            if "`replace'" != "" {
                capture confirm variable `v', exact
                if !_rc drop `v'
            }
            rename `t' `v'
            label variable `v' `"svytmle: `: word `k' of `labs''"'
        }
        if `est' & "`replace'" != "" {
            foreach s in Q0star Q1star {
                capture confirm variable `generate'`s', exact
                if !_rc {
                    local vl : variable label `generate'`s'
                    if substr(`"`vl'"', 1, 8) == "svytmle:" drop `generate'`s'
                }
            }
        }
    }

    matrix colnames `b' = ATE POmean1 POmean0
    matrix colnames `V' = ATE POmean1 POmean0
    matrix rownames `V' = ATE POmean1 POmean0
    tempvar es
    quietly gen byte `es' = `D'
    ereturn post `b' `V', depname(`y') obs(`N') dof(`=max(1, `df')') esample(`es')
    ereturn scalar N_sub        = `Nsub'
    ereturn scalar N_strata     = `n_strata'
    ereturn scalar N_clust      = `n_psu'
    ereturn scalar N_strata_sub = `n_strata_sub'
    ereturn scalar N_clust_sub  = `n_psu_sub'
    ereturn scalar N_strata_df  = `n_strata_df'
    ereturn scalar N_clust_df   = `n_psu_df'
    ereturn scalar df_design    = `df'
    ereturn scalar N_pop        = `N_pop'
    ereturn scalar N_subpop     = `N_subpop'
    ereturn scalar folds        = `folds'
    ereturn scalar V_eff        = `V_eff'
    ereturn scalar fold_exact   = `fold_exact'
    ereturn scalar n_uneven     = `n_uneven'
    ereturn scalar strictfolds  = `strict'
    ereturn scalar gbound       = `gbound'
    if !`est' {
        ereturn scalar gbound_tmle  = `gb_tmle'
        ereturn scalar n_gbound_tmle = `n_gbtmle'
        ereturn scalar eps0         = `eps0'
        ereturn scalar eps1         = `eps1'
    }
    ereturn scalar ps_min       = `g_min'
    ereturn scalar ps_max       = `g_max'
    ereturn scalar n_nonconv    = `n_nonconv'
    ereturn scalar rr           = `rr_est'
    ereturn scalar se_logrr     = `se_logrr'
    ereturn scalar or           = `or_est'
    ereturn scalar se_logor     = `se_logor'
    ereturn scalar rr_clamped   = `rr_clamped'
    ereturn scalar level        = `level'
    ereturn local treatment  `a'
    ereturn local covariates `covariates'
    ereturn local xterms     `xterms'
    ereturn local estimator  `estimator'
    ereturn local learner    logit
    ereturn local nuisance   `nuisance'
    ereturn local subpop     `"`subpopdesc'"'
    ereturn local foldvar    `foldvar'
    ereturn local seed       `"`seed'"'
    if "`foldvar'" == "" {
        ereturn local rngstate `"`rngstate'"'
        ereturn local rng      `rngname'
    }
    ereturn local wtype      pweight
    ereturn local wexp       `"= `wvar'"'
    ereturn local wvar       `wvar'
    ereturn local su1        `psu'
    ereturn local strata1    `strat'
    ereturn local generate   `gennames'
    ereturn local vce        linearized
    ereturn local vcetype    Linearized
    ereturn local predict      _no_predict
    ereturn local marginsnotok _ALL
    ereturn local estat_cmd    estat_vce_only
    ereturn local title      = cond(`est', "Survey-weighted cross-fitted AIPW", ///
        "Survey-weighted cross-fitted TMLE")
    ereturn matrix fold_psu = `FP'
    ereturn local cmdline    `"`cmdline'"'
    ereturn local cmdversion 0.2.1
    ereturn local cmd        svytmle

    Display, level(`level') `rr' `or' `header' `diopts'
end

program define Display
    version 15
    syntax [, Level(cilevel) RR OR noHEADer *]
    _get_diopts diopts, `options'
    if "`header'" == "" {
        * right-hand column width adapts to the largest count (svy-style)
        local w = max(11, min(17, strlen(string(max(e(N_pop), e(N)), "%24.0fc"))))
        local c3 = 78 - `w'
        local c2 = `c3' - 18
        local f %`w'.0fc
        local dv = abbrev("`e(depvar)'", 20)
        local tv = abbrev("`e(treatment)'", 20)
        local L1 "`e(title)'"
        local L2 "Estimand:  ATE (risk difference)"
        local L3 "Outcome:   `dv'"
        local L4 "Treatment: `tv'"
        local L5 "Nuisance:  logit, `e(nuisance)'"
        local L6 "Cross-fit: `e(V_eff)' folds of whole PSUs"
        local ws = cond(e(N_strata) > 1, "within strata", "")
        if missing(e(fold_exact)) local bal
        else if e(fold_exact) == 1 local bal "balanced"
        else local bal = "`e(n_uneven)' " + cond(e(N_strata) == 1, "uneven",            ///
            cond(e(n_uneven) == 1, "stratum uneven", "strata uneven"))
        if "`ws'" != "" & "`bal'" != "" local L7 "           `ws', `bal'"
        else local L7 "           `ws'`bal'"
        local R1 "Number of strata"
        local R2 "Number of PSUs"
        local R3 "Number of obs"
        local R4 "Population size"
        local V1 = e(N_strata)
        local V2 = e(N_clust)
        local V3 = e(N)
        local V4 = e(N_pop)
        if `"`e(subpop)'"' != "" {
            local R5 "Subpop. no. obs"
            local R6 "Subpop. size"
            local R7 "Design df"
            local V5 = e(N_sub)
            local V6 = e(N_subpop)
            local V7 = cond(missing(e(df_design)), e(df_r), e(df_design))
        }
        else {
            local R5 "Design df"
            local V5 = cond(missing(e(df_design)), e(df_r), e(df_design))
        }
        di
        forvalues i = 1/7 {
            if `"`R`i''"' != "" {
                di as txt `"`L`i''"' _col(`c2') "`R`i''" _col(`c3') "= " as res `f' `V`i''
            }
            else di as txt `"`L`i''"'
        }
        if `"`e(subpop)'"' != "" {
            if udstrlen(`"`e(subpop)'"') <= 66 {
                di as txt "Subpop.:" _col(12) as res `"`e(subpop)'"'
            }
            else di as txt "{p 0 11 2}Subpop.: " as res `"`e(subpop)'"' "{p_end}"
        }
        di
    }
    ereturn display, level(`level') `diopts'
    if e(df_design) != e(N_clust) - e(N_strata) {
        di as txt "{p 0 6 2}Note: design df = " e(N_clust_df) " - " e(N_strata_df)       ///
            " = " e(df_design) ": only PSUs and strata containing subpopulation "         ///
            "observations with nonzero weight are counted (R survey::degf; Stata 19's "   ///
            "{bf:dofsubpop}).{p_end}"
    }
    if "`rr'`or'" != "" {
        if e(rr_clamped) == 1 {
            di as txt _n "{p 0 6 2}Relative effects are not defined for these results: a "  ///
                "potential-outcome mean lies outside (0, 1).{p_end}"
        }
        else if missing(e(rr)) | missing(e(or)) {
            di as txt _n "{p 0 6 2}Relative effects are not stored in these results "        ///
                "(for example after {bf:mi estimate}); obtain them with {bf:nlcom} or "       ///
                "{bf:mi estimate (logRR: ln(_b[POmean1]/_b[POmean0])): svytmle ...}{p_end}"
        }
        else {
            tempname crit
            scalar `crit' = invttail(e(df_r), (1 - `level'/100)/2)
            di as txt _n "Relative effects (delta method on the log scale, t with design df):"
            if "`rr'" != "" {
                di as txt "  Risk ratio  " as res %9.4f e(rr) as txt "   [`level'% CI "    ///
                    as res %9.4f exp(ln(e(rr)) - `crit'*e(se_logrr)) as txt ", "        ///
                    as res %9.4f exp(ln(e(rr)) + `crit'*e(se_logrr)) as txt "]"
            }
            if "`or'" != "" {
                di as txt "  Odds ratio  " as res %9.4f e(or) as txt "   [`level'% CI "    ///
                    as res %9.4f exp(ln(e(or)) - `crit'*e(se_logor)) as txt ", "        ///
                    as res %9.4f exp(ln(e(or)) + `crit'*e(se_logor)) as txt "]"
            }
        }
    }
end

version 15
mata:

// ---------------------------------------------------------------------------
// R glm.fit's logit link helpers (C code in R's family.c): linkinv and mu.eta
// clamp |eta| > 30 to DBL_EPSILON-based values.
// ---------------------------------------------------------------------------
real colvector _svytmle_linkinv(real colvector eta)
{
    real colvector inr, t
    inr = (eta :>= -30) :& (eta :<= 30)
    t = exp(eta :* inr)
    t = inr :* t + (eta :< -30) :* epsilon(1) + (eta :> 30) :* (1 / epsilon(1))
    return(t :/ (1 :+ t))
}

real colvector _svytmle_mueta(real colvector eta)
{
    real colvector inr, t
    inr = (eta :>= -30) :& (eta :<= 30)
    t = exp(eta :* inr)
    return(inr :* (t :/ ((1 :+ t) :^ 2)) + (1 :- inr) :* epsilon(1))
}

real scalar _svytmle_dev(real colvector y, real colvector mu, real colvector wt)
{
    // binomial deviance for y in [0,1] (R binomial_dev_resids, y_log_y)
    real colvector t1, t0
    t1 = y :* ln(mu) :* (y :> 0)
    t0 = (1 :- y) :* ln(1 :- mu) :* (y :< 1)
    t1 = editmissing(t1, 0)
    t0 = editmissing(t0, 0)
    return(-2 * sum(wt :* (t1 + t0)) + 2 * sum(wt :* _svytmle_ylogy(y)))
}

real colvector _svytmle_ylogy(real colvector y)
{
    // y*log(y) + (1-y)*log(1-y), 0 at the boundaries (0 for binary y)
    real colvector r
    r = y :* ln(y) + (1 :- y) :* ln(1 :- y)
    return(editmissing(r, 0))
}

// ---------------------------------------------------------------------------
// Logistic regression by IRLS, replicating R's glm.fit (binomial, logit):
// mustart = (w*y + 0.5)/(w + 1), eta = logit(mustart) (offset NOT added at the
// start), convergence |dev - devold|/(|dev| + 0.1) < 1e-8, maxit = 25.
// Each step is a QR least-squares solve of the weighted working response on
// sqrt(W) X (as R's Cdqrls); X must have full column rank on the positive-
// weight rows, which the callers ensure.
// ---------------------------------------------------------------------------
real colvector _svytmle_irls(real colvector y, real matrix X, real colvector wt,
                             real colvector off, real scalar conv)
{
    real colvector mu, eta, z, wz, beta, mueta
    real scalar dev, devold, iter

    mu = (wt :* y :+ 0.5) :/ (wt :+ 1)
    eta = logit(mu)
    devold = _svytmle_dev(y, mu, wt)
    conv = 0
    beta = J(cols(X), 1, 0)
    for (iter = 1; iter <= 25; iter++) {
        mueta = _svytmle_mueta(eta)
        z = (eta - off) + (y - mu) :/ mueta
        wz = wt :* (mueta :^ 2) :/ (mu :* (1 :- mu))
        beta = qrsolve(sqrt(wz) :* X, sqrt(wz) :* z)
        eta = X * beta + off
        mu = _svytmle_linkinv(eta)
        dev = _svytmle_dev(y, mu, wt)
        if (abs(dev - devold) / (abs(dev) + 0.1) < 1e-8) {
            conv = 1
            break
        }
        devold = dev
    }
    return(beta)
}

// clip x to [lo, hi]
real colvector _svytmle_clip(real colvector x, real scalar lo, real scalar hi)
{
    return(rowmin((rowmax((x, J(rows(x), 1, lo))), J(rows(x), 1, hi))))
}

// ---------------------------------------------------------------------------
// Logistic nuisance fit on training rows, prediction on new rows.
// Aliasing as R's glm (dqrdc2 limited pivoting): columns constant in the
// training rows are dropped (collinear with the intercept); any remaining
// exact collinearity is resolved ONCE, before IRLS, by sweeping the weighted
// cross-product in natural column order, so the LATER of collinear columns is
// dropped, as R does. Kept columns are centred and scaled (an affine
// reparametrisation: the same IRLS iterates and predictions in exact
// arithmetic). ndrop = number of covariate columns dropped.
// ---------------------------------------------------------------------------
real matrix _svytmle_logit_fitpred(real colvector ytr, real matrix Xtr,
                                   real colvector wt, real matrix Xnew,
                                   real scalar conv, real scalar ndrop)
{
    real rowvector keep, cen, sc, nz
    real matrix Z, Zn, P
    real colvector beta, d
    real scalar ntr, j, k, nn

    ntr = rows(Xtr)
    keep = selectindex(colmax(Xtr) :> colmin(Xtr))
    if (cols(keep) > 0) {
        cen = mean(Xtr[., keep])
        sc  = sqrt(colsum((Xtr[., keep] :- cen) :^ 2) :/ (ntr - 1))
        Z   = J(ntr, 1, 1), ((Xtr[., keep] :- cen) :/ sc)
    }
    else Z = J(ntr, 1, 1)
    d  = diagonal(invsym(cross(Z, wt, Z), (1..cols(Z))))
    nz = selectindex(d' :!= 0)
    Z  = Z[., nz]
    ndrop = cols(Xtr) - (cols(nz) - 1)
    beta = _svytmle_irls(ytr, Z, wt, J(ntr, 1, 0), conv)
    // Xnew holds one or more prediction matrices side by side, each with
    // cols(Xtr) columns
    k = cols(Xnew) / cols(Xtr)
    nn = rows(Xnew)
    P = J(nn, k, .)
    for (j = 1; j <= k; j++) {
        Zn = Xnew[., ((j - 1) * cols(Xtr) + 1)..(j * cols(Xtr))]
        if (cols(keep) > 0) Zn = J(nn, 1, 1), ((Zn[., keep] :- cen) :/ sc)
        else Zn = J(nn, 1, 1)
        P[., j] = _svytmle_linkinv(Zn[., nz] * beta)
    }
    return(P)
}

// ---------------------------------------------------------------------------
// Cross-fitting folds (R make_cf_folds): whole PSUs within strata,
// V_eff = max(2, min(V, min_h m_h)) over the domain's strata; within each
// stratum a random permutation of rep_len(1..V_eff, m_h).
// ---------------------------------------------------------------------------
real colvector _svytmle_folds(real colvector hd, real colvector pd,
                              real scalar V, real scalar Veff)
{
    real matrix up, info
    real colvector mh, fpsu, lookup, pat
    real scalar k, m

    up = uniqrows((hd, pd))
    info = panelsetup(up, 1)
    mh = info[., 2] - info[., 1] :+ 1
    Veff = max((2, min((V, min(mh)))))
    fpsu = J(rows(up), 1, .)
    for (k = 1; k <= rows(info); k++) {
        m = mh[k]
        pat = mod((0::(m - 1)), Veff) :+ 1
        fpsu[|info[k, 1] \ info[k, 2]|] = jumble(pat)
    }
    lookup = J(max(up[., 2]), 1, .)
    lookup[up[., 2], 1] = fpsu
    return(lookup[pd, 1])
}

// number of domain strata with no PSU in the training rows of fold v;
// nmulti = how many of those have >= 2 domain PSUs (absent because of the
// supplied folds, not because the stratum has a single domain PSU)
real scalar _svytmle_nmiss(real colvector hd, real colvector pd,
                           real colvector fold, real scalar v, real scalar nmulti)
{
    real colvector tr, inh, intr, mh, miss
    real matrix up
    real scalar k
    tr = selectindex(fold :!= v)
    inh = J(max(hd), 1, 0)
    intr = J(max(hd), 1, 0)
    inh[hd, 1] = J(rows(hd), 1, 1)
    intr[hd[tr], 1] = J(rows(tr), 1, 1)
    mh = J(max(hd), 1, 0)
    up = uniqrows((hd, pd))
    for (k = 1; k <= rows(up); k++) mh[up[k, 1]] = mh[up[k, 1]] + 1
    miss = inh :& !intr
    nmulti = sum(miss :& (mh :>= 2))
    return(sum(miss))
}

// training-fold survey weights (R .cf_train_weights): w * m_h / m_train,h,
// normalised to mean 1 over the training rows
real colvector _svytmle_trainwt(real colvector wd, real colvector hd,
                                real colvector pd, real colvector fold,
                                real scalar v)
{
    real colvector tr, mh, mtr, wt
    real matrix up, ut
    real scalar k

    tr = selectindex(fold :!= v)
    mh = J(max(hd), 1, 0)
    mtr = J(max(hd), 1, 0)
    up = uniqrows((hd, pd))
    for (k = 1; k <= rows(up); k++) mh[up[k, 1]] = mh[up[k, 1]] + 1
    ut = uniqrows((hd[tr], pd[tr]))
    for (k = 1; k <= rows(ut); k++) mtr[ut[k, 1]] = mtr[ut[k, 1]] + 1
    wt = wd[tr] :* (mh[hd[tr], 1] :/ mtr[hd[tr], 1])
    return(wt :/ mean(wt))
}

// greatest common divisor of two positive integers
real scalar _svytmle_gcd(real scalar a, real scalar b)
{
    real scalar x, y, t
    x = a
    y = b
    while (y > 0) {
        t = mod(x, y)
        x = y
        y = t
    }
    return(x)
}

// ---------------------------------------------------------------------------
// Fold balance within strata: a stratum is evenly split when every fold holds
// the same number of its domain PSUs (for random folds: V divides m_h).
// exact = 1 when every stratum is evenly split; cnt = PSUs per fold (summed
// over strata), in the order of the sorted fold labels.
// ---------------------------------------------------------------------------
void _svytmle_foldbal(real colvector hd, real colvector pd, real colvector fold,
                      real colvector folds, real scalar exact,
                      real scalar nuneven, real rowvector cnt)
{
    real matrix up, info
    real rowvector c
    real scalar k, j, V, q

    up = uniqrows((hd, pd, fold))
    info = panelsetup(up, 1)
    V = rows(folds)
    nuneven = 0
    cnt = J(1, V, 0)
    for (k = 1; k <= rows(info); k++) {
        c = J(1, V, 0)
        for (j = info[k, 1]; j <= info[k, 2]; j++) {
            q = selectindex(folds :== up[j, 3])
            c[q] = c[q] + 1
        }
        cnt = cnt + c
        if (max(c) != min(c)) nuneven = nuneven + 1
    }
    exact = (nuneven == 0)
}

// ---------------------------------------------------------------------------
// Stratified first-stage with-replacement linearization covariance of the
// domain Hajek means of the columns of E (zero outside the domain): PSU totals
// of w*s*(e - mean)/sum(w*s), centred within strata over ALL design PSUs of the
// stratum, scaled by n_h/(n_h - 1) (R svymean on a subset() design).
// ---------------------------------------------------------------------------
real matrix _svytmle_linvar(real matrix E, real colvector h, real colvector p,
                            real colvector w, real colvector s)
{
    real colvector ws, psuh, o, hs
    real rowvector m
    real matrix U, T, Ts, Tj, C, info, Sig
    real scalar Nhat, i, j, nh

    ws = w :* s
    Nhat = sum(ws)
    m = cross(ws, E) :/ Nhat
    U = (ws :* (E :- m)) :/ Nhat
    T = J(max(p), cols(E), 0)
    psuh = J(max(p), 1, .)
    for (i = 1; i <= rows(p); i++) {
        T[p[i], .] = T[p[i], .] + U[i, .]
        psuh[p[i]] = h[i]
    }
    o = order(psuh, 1)
    Ts = T[o, .]
    hs = psuh[o]
    info = panelsetup(hs, 1)
    Sig = J(cols(E), cols(E), 0)
    for (j = 1; j <= rows(info); j++) {
        Tj = Ts[|info[j, 1], 1 \ info[j, 2], cols(E)|]
        nh = rows(Tj)
        if (nh < 2) continue
        C = Tj :- mean(Tj)
        Sig = Sig + (nh / (nh - 1)) :* cross(C, C)
    }
    return(Sig)
}

void _svytmle_putlocal(string scalar name, real scalar x)
{
    st_local(name, strofreal(x, "%21.0g"))
}

// ---------------------------------------------------------------------------
// main driver
// ---------------------------------------------------------------------------
void _svytmle_main(string scalar yv, string scalar av, string scalar xv,
                   string scalar hv, string scalar pv, string scalar wv,
                   string scalar dv, string scalar sv, string scalar foldv,
                   real scalar V, real scalar gbound, real scalar nuisw,
                   real scalar est, real scalar strict,
                   string scalar bname, string scalar Vname,
                   string scalar FPname, string scalar genlist)
{
    real colvector h, p, w, s, y, a, hd, pd, wd, fold, folds, tr, ho, wt
    real colvector Q0, Q1, g1, Q0b, Q1b, Q0t, Q1t, QAt, off, g1den, g0den
    real colvector H0, H1, omega, eps, Q0s, Q1s, eif1, eif0, eif, idx, dom
    real matrix X, XA, Xho, PQ, Pg, hp, info, E, Sig, A, Vm, vals
    real scalar n, Veff, k, v, conv, nnc, nmt, nd1, nd2, ndf, gb, sw
    real scalar psi, psi1, psi0, df, p1, p0, lrr, lor, selrr, selor, j
    real scalar exact, nuneven, gg, nbad, dmin, nm, nmm, clamp
    real colvector grr, gor, D1, D0, mhd
    real rowvector fcnt
    real matrix upd, infod
    string rowvector gn

    // ---- full design sample
    h = st_data(., hv, dv)
    p = st_data(., pv, dv)
    w = st_data(., wv, dv)
    s = st_data(., sv, dv)

    // ---- every stratum of the design needs >= 2 PSUs
    hp = uniqrows((h, p))
    info = panelsetup(hp, 1)
    _svytmle_putlocal("n_single", sum((info[., 2] - info[., 1]) :== 0))
    _svytmle_putlocal("n_strata", rows(info))
    _svytmle_putlocal("n_psu", rows(hp))
    _svytmle_putlocal("n_noncon", 0)
    if (sum((info[., 2] - info[., 1]) :== 0) > 0) return

    // ---- domain data
    y  = st_data(., yv, sv)
    a  = st_data(., av, sv)
    X  = st_data(., tokens(xv), sv)
    hd = st_data(., hv, sv)
    pd = st_data(., pv, sv)
    wd = st_data(., wv, sv)
    n  = rows(y)
    _svytmle_putlocal("n_strata_sub", rows(uniqrows(hd)))
    _svytmle_putlocal("n_psu_sub", rows(uniqrows(pd)))
    _svytmle_putlocal("N_pop", sum(w))
    _svytmle_putlocal("N_subpop", sum(wd))
    if (cols(selectindex(colmax(X) :> colmin(X))) == 0) {
        _svytmle_putlocal("n_noncon", 1)
        return
    }

    // ---- folds
    if (foldv != "") {
        fold = st_data(., foldv, sv)
        if (rows(uniqrows((pd, fold))) != rows(uniqrows(pd))) {
            errprintf("foldvar() must be constant within each PSU\n")
            exit(459)
        }
        Veff = rows(uniqrows(fold))
    }
    else {
        if (strict) {
            // exact balance: V must divide every domain stratum's PSU count
            upd = uniqrows((hd, pd))
            infod = panelsetup(upd, 1)
            mhd = infod[., 2] - infod[., 1] :+ 1
            nbad = sum(mod(mhd, V) :!= 0)
            if (nbad > 0) {
                gg = mhd[1]
                for (j = 2; j <= rows(mhd); j++) gg = _svytmle_gcd(gg, mhd[j])
                errprintf("strictfolds: folds(%g) does not divide the number of domain PSUs\n", V)
                errprintf("in %g of %g strata, so the folds cannot be exactly balanced.\n",
                          nbad, rows(mhd))
                if (gg >= 2) {
                    dmin = 2
                    while (mod(gg, dmin) != 0) dmin = dmin + 1
                    errprintf("The domain PSU counts have greatest common divisor %g:\n", gg)
                    errprintf("folds(%g) (or any divisor of %g of at least 2) is exact.\n", dmin, gg)
                }
                else errprintf("No number of folds splits every stratum evenly.\n")
                exit(459)
            }
        }
        fold = _svytmle_folds(hd, pd, V, Veff)
    }
    folds = uniqrows(fold)
    if (rows(folds) < 2 & foldv != "") {
        errprintf("foldvar() defines only one fold in the analytic domain; cross-fitting\n")
        errprintf("needs at least 2 folds\n")
        exit(459)
    }
    if (rows(folds) < 2) {
        errprintf("cross-fitting needs at least 2 folds of PSUs, but every stratum of\n")
        errprintf("the analytic domain contains a single domain PSU: PSU-level\n")
        errprintf("cross-fitting is not possible for this domain\n")
        exit(459)
    }
    _svytmle_foldbal(hd, pd, fold, folds, exact, nuneven, fcnt)
    if (strict & !exact) {
        errprintf("strictfolds: the folds in foldvar() are not exactly balanced: in %g\n", nuneven)
        errprintf("strata the domain PSUs are not split evenly across the %g folds\n", rows(folds))
        exit(459)
    }

    // ---- out-of-fold nuisances (logistic regression)
    XA = (a, X)
    Q0 = J(n, 1, .)
    Q1 = J(n, 1, .)
    g1 = J(n, 1, .)
    nnc = 0
    nmt = 0
    nmm = 0
    ndf = 0
    for (k = 1; k <= rows(folds); k++) {
        v  = folds[k]
        tr = selectindex(fold :!= v)
        ho = selectindex(fold :== v)
        nmt = nmt + _svytmle_nmiss(hd, pd, fold, v, nm)
        nmm = nmm + nm
        if (nuisw) {
            if (sum(wd[tr]) <= 0) {
                errprintf("nuisance(survey): the training rows of fold %g have zero total weight\n", v)
                exit(459)
            }
            wt = _svytmle_trainwt(wd, hd, pd, fold, v)
        }
        else wt = J(rows(tr), 1, 1)
        Xho = X[ho, .]
        PQ = _svytmle_logit_fitpred(y[tr], XA[tr, .], wt,
                 ((J(rows(ho), 1, 1), Xho), (J(rows(ho), 1, 0), Xho)), conv, nd1)
        nnc = nnc + (conv == 0)
        Pg = _svytmle_logit_fitpred(a[tr], X[tr, .], wt, Xho, conv, nd2)
        nnc = nnc + (conv == 0)
        ndf = ndf + (nd1 > 0) + (nd2 > 0)
        Q1[ho] = PQ[., 1]
        Q0[ho] = PQ[., 2]
        g1[ho] = Pg[., 1]
    }

    // ---- truncation of the out-of-fold predictions (R run_estimators)
    g1  = _svytmle_clip(g1, gbound, 1 - gbound)
    Q0b = _svytmle_clip(Q0, 1e-3, 1 - 1e-3)
    Q1b = _svytmle_clip(Q1, 1e-3, 1 - 1e-3)

    sw = sum(wd)
    if (est) {
        // ---- cross-fitted AIPW (R svyaipw AIPW-CF): Hajek means of the
        // pseudo-outcomes; no targeting and no tmle lower bound
        D1 = (a :/ g1) :* (y - Q1b) + Q1b
        D0 = ((1 :- a) :/ (1 :- g1)) :* (y - Q0b) + Q0b
        psi1 = sum(wd :* D1) / sw
        psi0 = sum(wd :* D0) / sw
        psi  = sum(wd :* (D1 - D0)) / sw
        eif1 = D1 :- psi1
        eif0 = D0 :- psi0
        eif  = (D1 - D0) :- psi
        gb = .
        eps = (. \ .)
    }
    else {
        // ---- targeting (replicates tmle 2.1.1 with user-supplied Q and g1W)
        gb = 5 / sqrt(n) / ln(n)                  // tmle's default lower-only g bound
        g1den = rowmax((g1, J(n, 1, gb)))
        g0den = rowmax((1 :- g1, J(n, 1, gb)))
        Q0t = _svytmle_clip(Q0b, 1 - 0.9995, 0.9995)
        Q1t = _svytmle_clip(Q1b, 1 - 0.9995, 0.9995)
        QAt = (1 :- a) :* Q0t + a :* Q1t
        off = logit(QAt)
        H0 = (1 :- a) :/ g0den
        H1 = a :/ g1den
        omega = wd :/ sum(wd) :* n
        eps = _svytmle_irls(y, (H0, H1), omega, off, conv)
        nnc = nnc + (conv == 0)
        eps = editmissing(eps, 0)
        Q0s = invlogit(logit(Q0t) + eps[1] :/ g0den)
        Q1s = invlogit(logit(Q1t) + eps[2] :/ g1den)

        // ---- estimates and per-unit (unweighted) influence functions
        psi  = sum(wd :* (Q1s - Q0s)) / sw
        psi1 = sum(wd :* Q1s) / sw
        psi0 = sum(wd :* Q0s) / sw
        eif1 = (a :/ g1den) :* (y - Q1s) + Q1s :- psi1
        eif0 = ((1 :- a) :/ g0den) :* (y - Q0s) + Q0s :- psi0
        eif  = (a :/ g1den - (1 :- a) :/ g0den) :* (y - ((1 :- a) :* Q0s + a :* Q1s)) +
               (Q1s - Q0s) :- psi
    }

    // ---- design covariance of (psi1, psi0) and df
    E = J(rows(h), 2, 0)
    idx = selectindex(s)
    E[idx, .] = (eif1, eif0)
    Sig = _svytmle_linvar(E, h, p, w, s)
    A = (1, -1 \ 1, 0 \ 0, 1)
    Vm = A * Sig * A'
    Vm = (Vm + Vm') :/ 2
    dom = selectindex(s :& (w :!= 0))
    df = rows(uniqrows(p[dom])) - rows(uniqrows(h[dom]))

    // ---- relative effects (R delta_rr_or, eps = 1e-6 clamp); a PO mean
    // outside the clamp range makes the ratios meaningless: report missing
    clamp = (psi1 < 1e-6) | (psi1 > 1 - 1e-6) | (psi0 < 1e-6) | (psi0 > 1 - 1e-6)
    p1 = min((max((psi1, 1e-6)), 1 - 1e-6))
    p0 = min((max((psi0, 1e-6)), 1 - 1e-6))
    lrr = ln(p1) - ln(p0)
    grr = (1 / p1 \ -1 / p0)
    selrr = sqrt(max((grr' * Sig * grr, 0)))
    lor = logit(p1) - logit(p0)
    gor = (1 / (p1 * (1 - p1)) \ -1 / (p0 * (1 - p0)))
    selor = sqrt(max((gor' * Sig * gor, 0)))

    // ---- results back to Stata
    st_matrix(bname, (psi, psi1, psi0))
    st_matrix(Vname, Vm)
    _svytmle_putlocal("df", df)
    _svytmle_putlocal("n_psu_df", rows(uniqrows(p[dom])))
    _svytmle_putlocal("n_strata_df", rows(uniqrows(h[dom])))
    _svytmle_putlocal("V_eff", Veff)
    _svytmle_putlocal("fold_exact", exact)
    _svytmle_putlocal("n_uneven", nuneven)
    st_matrix(FPname, fcnt)
    st_matrixcolstripe(FPname, (J(rows(folds), 1, ""), "fold" :+ strofreal(folds, "%21.0g")))
    st_matrixrowstripe(FPname, ("", "PSUs"))
    _svytmle_putlocal("gb_tmle", gb)
    _svytmle_putlocal("n_gbtmle", est ? . : sum((g1 :< gb) :| ((1 :- g1) :< gb)))
    _svytmle_putlocal("eps0", eps[1])
    _svytmle_putlocal("eps1", eps[2])
    _svytmle_putlocal("g_min", min(g1))
    _svytmle_putlocal("g_max", max(g1))
    _svytmle_putlocal("n_nonconv", nnc)
    _svytmle_putlocal("n_notrain", nmt - nmm)
    _svytmle_putlocal("n_notrain_fv", nmm)
    _svytmle_putlocal("n_dropfits", ndf)
    _svytmle_putlocal("rr_clamped", clamp)
    _svytmle_putlocal("rr_est", clamp ? . : exp(lrr))
    _svytmle_putlocal("se_logrr", clamp ? . : selrr)
    _svytmle_putlocal("or_est", clamp ? . : exp(lor))
    _svytmle_putlocal("se_logor", clamp ? . : selor)

    if (genlist != "") {
        gn = tokens(genlist)
        if (est) vals = (fold, g1, Q0b, Q1b, eif)
        else vals = (fold, g1, Q0b, Q1b, Q0s, Q1s, eif)
        for (j = 1; j <= cols(gn); j++) {
            (void) st_addvar("double", gn[j])
            st_store(., gn[j], sv, vals[., j])
        }
    }
}

end
