{smcl}
{* *! version 0.2.1  16oct2026}{...}
{vieweralsosee "[SVY] svyset" "help svyset"}{...}
{vieweralsosee "[SVY] svy estimation" "help svy_estimation"}{...}
{vieweralsosee "[SVY] subpopulation estimation" "mansection SVY subpopulationestimation"}{...}
{vieweralsosee "[MI] mi estimate" "help mi_estimate"}{...}
{viewerjumpto "Syntax" "svytmle##syntax"}{...}
{viewerjumpto "Description" "svytmle##description"}{...}
{viewerjumpto "Options" "svytmle##options"}{...}
{viewerjumpto "Method" "svytmle##method"}{...}
{viewerjumpto "Fold balance" "svytmle##foldbalance"}{...}
{viewerjumpto "Postestimation" "svytmle##postestimation"}{...}
{viewerjumpto "Examples" "svytmle##examples"}{...}
{viewerjumpto "Stored results" "svytmle##results"}{...}
{viewerjumpto "Validation" "svytmle##validation"}{...}
{viewerjumpto "References" "svytmle##references"}{...}
{title:Title}

{phang}
{bf:svytmle} {hline 2} Survey-weighted cross-fitted TMLE (or AIPW) of the
population average treatment effect with a design-based variance


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:svytmle} {depvar}{cmd:,}
{opt t:reatment(varname)}
{opt cov:ariates(varlist)}
[{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt t:reatment(varname)}}binary (0/1) treatment or exposure{p_end}
{synopt:{opt cov:ariates(varlist)}}adjustment covariates; factor variables
allowed{p_end}

{syntab:Domain}
{synopt:{cmd:subpop(}[{varname}] [{it:{help if}}]{cmd:)}}analytic domain
(subpopulation){p_end}

{syntab:Estimator}
{synopt:{opt est:imator(type)}}{cmd:tmle} (default) or {cmd:aipw}{p_end}

{syntab:Cross-fitting}
{synopt:{opt f:olds(#)}}target number of cross-fitting folds; default
{cmd:5}{p_end}
{synopt:{opt strict:folds}}require exactly balanced folds{p_end}
{synopt:{opt seed(#)}}set the random-number seed before the folds are
drawn{p_end}
{synopt:{opt foldv:ar(varname)}}use these fold numbers instead of random
folds{p_end}

{syntab:Nuisance models}
{synopt:{opt nuis:ance(type)}}{cmd:unweighted} (default) or {cmd:survey}{p_end}
{synopt:{opt gb:ound(#)}}truncate out-of-fold propensity scores to
[{it:#}, 1-{it:#}]; default {cmd:0.05}{p_end}

{syntab:Reporting}
{synopt:{opt l:evel(#)}}confidence level; default is {cmd:level(95)}{p_end}
{synopt:{opt rr}}also report the risk ratio{p_end}
{synopt:{opt or}}also report the odds ratio{p_end}
{synopt:{opt gen:erate(stub)}}save folds, nuisance predictions, and the
influence function{p_end}
{synopt:{opt replace}}overwrite the {opt generate()} variables{p_end}
{synopt:{opt nohead:er}}suppress the output header{p_end}
{synopt:{it:display_options}}{opt coefl:egend}, {opt cformat(%fmt)}, and other
{help estimation options##display_options:display options}{p_end}
{synoptline}
{p 4 6 2}
The data must be {helpb svyset} with a {cmd:pweight}.
{cmd:if} and {cmd:in} are not allowed; use {cmd:subpop()}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:svytmle} estimates the population average treatment effect (risk
difference) of a binary treatment on a binary outcome from complex survey
data. It uses a targeted maximum likelihood estimator (TMLE) whose nuisance
models are cross-fitted at the level of the first-stage sampling unit (PSU),
with survey-weighted targeting and a design-based (Taylor-series
linearization) variance that treats the PSU as the independent unit.
Confidence intervals use a t reference with the survey design degrees of
freedom.

{pstd}
The command reads the design from {helpb svyset}: the first-stage PSU
variable, the strata variable, and the {cmd:pweight}. It does not change the
svyset settings. PSU identifiers that repeat across strata (for example
NHANES {cmd:SDMVPSU}) are handled as nested within strata.

{pstd}
{depvar} and {opt treatment()} must be coded 0/1 and be nonmissing in the
analytic domain, with both values present. Unlike {cmd:svy} estimation
commands, {cmd:svytmle} never drops observations from a domain. If the
outcome, the treatment, a covariate, or a design variable is missing for an
observation in the domain, it stops with an error. Define the domain as
complete cases (for example {cmd:subpop(inpop)}) or impute first (see
{help svytmle##examples:Examples}).

{pstd}
With {cmd:estimator(aipw)}, the same cross-fitted nuisance models and the same
design-based variance give the cross-fitted augmented inverse-probability
weighted (AIPW) estimator instead. In the double machine learning (DML)
literature this is the cross-fitted estimator based on the efficient
(AIPW) score for the ATE.

{pstd}
{cmd:svytmle} is the Stata implementation of the R package {bf:svytmle}
(Karim 2026). With the same folds, it reproduces the R package's primary
estimator ({cmd:svytmle(..., learners = "SL.glm", arms = "cf")}) and the ATE,
SE, df, and CI of its cross-fitted AIPW comparator ({cmd:svyaipw(..., learners = "SL.glm",}
{cmd:arms = "cf")}). {cmd:svyaipw()} reports only the ATE; the AIPW
potential-outcome means and relative effects are checked against a reference
computed from the R package's internal functions. See
{help svytmle##validation:Validation}.

{pstd}
The file
{browse "https://github.com/ehsanx/svytmle-stata/blob/main/R_VS_STATA.md":R_VS_STATA.md}
in the package's source repository,
{browse "https://github.com/ehsanx/svytmle-stata"}, lists every
difference between the R and Stata implementations: what is identical, how to
reproduce an R analysis, and the features present on only one side.

{pstd}
This version fits the nuisance models (the outcome regression and the
propensity score) by main-terms logistic regression. Flexibility comes from
the covariate specification: factor variables, interactions, and polynomial
terms are all allowed in {opt covariates()}.


{marker options}{...}
{title:Options}

{phang}
{opt treatment(varname)} specifies the treatment, coded 0/1.

{phang}
{opt estimator(type)} selects the estimator. {cmd:tmle}, the default, is the
cross-fitted TMLE. {cmd:aipw} is the cross-fitted AIPW (one-step) estimator,
which uses the same folds and out-of-fold nuisance models but no targeting
step. With the same {opt seed()} or {opt foldvar()}, the two estimators use
identical folds, so their results can be compared directly. (In R, by
contrast, {cmd:svytmle()} and {cmd:svyaipw()} draw different folds under the
same seed; see {opt foldvar()}.)

{phang}
{opt covariates(varlist)} specifies the adjustment covariates. Factor-variable
notation ({cmd:i.}, {cmd:c.}, {cmd:#}, {cmd:##}) is allowed. Factor variables
are indicator-coded with the lowest level present in the domain as the base.
Within each training fold, as in R's {cmd:glm}, columns that are constant are
dropped, and of a set of exactly collinear columns the later ones are dropped
(their coefficients are set to zero); a note reports how many fits this
affected. Nearly collinear columns, whose variance not explained by the other
columns is below about 1e-9 of their own, are also dropped. R's {cmd:glm}
keeps such columns down to a much smaller tolerance and fits a numerically
unstable direction. Covariates must be nonmissing in the analytic domain.

{phang}
{cmd:subpop(}[{varname}] [{it:if}]{cmd:)} restricts estimation to an analytic
domain. Observations with {it:varname} equal to 0, or for which {it:if} is
false, are outside the domain. The nuisance models and the targeting use the
domain, while the variance is computed on the full design: PSUs that only
partly intersect the domain contribute correctly ("subset the design, not the
data"; see {manlink SVY subpopulation estimation}). Two differences from
{cmd:svy, subpop()}:

{phang2}
o A missing value of {it:varname} within the design is an error. {cmd:svy}
would drop those observations from the design. (An {cmd:if} expression, by
contrast, follows Stata's usual rules, exactly as in {cmd:svy, subpop()}.
Missing values count as larger than any number, so {cmd:subpop(if bmi >= 25)}
includes observations with missing {cmd:bmi}. Write
{cmd:subpop(if bmi >= 25 & !missing(bmi))}. R's {cmd:subset} argument stops
with an error instead.)

{phang2}
o The design degrees of freedom count only the PSUs and strata that contain
domain observations with nonzero weight, as R's {cmd:survey::degf} does and as
Stata 19's {cmd:dofsubpop} option does. Stata's default would count all PSUs
of the strata that contain domain observations. The variance itself is the
same as {cmd:svy, subpop(): mean} of the influence function.

{phang}
{opt folds(#)} sets the target number of cross-fitting folds V; the default
is 5. Folds consist of whole PSUs within strata. V is reduced to the smallest
number of domain PSUs in any stratum (but never below 2), so every training
fold contains every stratum that has at least two domain PSUs. With 2 PSUs per
stratum (NHANES), this is leave-one-PSU-per-stratum-out cross-fitting. A
stratum with a single domain PSU is missing from one training fold, and a
warning says so. If every stratum of the domain has a single domain PSU,
cross-fitting is impossible and {cmd:svytmle} stops with an error.

{phang}
{opt strictfolds} requires exactly balanced folds. V must divide the number of
domain PSUs in every stratum, and V is never reduced. If this fails,
{cmd:svytmle} stops with an error. The error suggests a number of folds that
would work when one exists, that is, when the PSU counts share a divisor of at
least 2. In pooled NHANES, strata have 2 or 3 PSUs, so none does. With
{opt foldvar()}, the supplied folds must split every stratum's domain PSUs
evenly. See {help svytmle##foldbalance:Fold balance}.

{phang}
{opt seed(#)} sets the random-number seed (see {helpb set seed}) before the
folds are drawn. As for any Stata command, the generator depends on the
{helpb version} in effect (see {helpb set rng}), so a seed reproduces the
folds under the same version and {cmd:set rng} setting. The state of the
generator before the folds are drawn is stored in {cmd:e(rngstate)}, so any
run, seeded or not, can be reproduced with
{cmd:set rngstate `e(rngstate)'}. In observation-level designs
({cmd:svyset _n}) the folds also depend on the sort order of the data. Stata's
and R's generators differ, so the same seed does not give the same folds in
Stata and in R; use {opt foldvar()} to reproduce an R analysis exactly.

{phang}
{opt foldvar(varname)} supplies the fold of each observation in the analytic
domain (positive integers, constant within each PSU, at least 2 folds). It
overrides {opt folds()} and {opt seed()}. To reproduce an R analysis, export
the folds R used.

{phang2}
o {cmd:svytmle()} draws its folds first, so
{cmd:set.seed(s); svytmle:::make_cf_folds(strata, psu, V)} on the domain rows
reproduces them.

{phang2}
o {cmd:svyaipw()} draws its folds only after its single-fit models, and
returns neither set. Its folds must therefore be captured, as the
certification script does by tracing {cmd:svytmle:::make_cf_folds}.

{phang}
{opt nuisance(type)} chooses how the out-of-fold nuisance models are fit.
{cmd:unweighted}, the default, fits them without survey weights. This is the
estimator of Karim (2026). It is consistent when sampling is non-informative
given the covariates, which is what the theory assumes, and it avoids the
separation that weighted fits can produce. {cmd:survey} fits them with the
survey weights, rescaled within each training fold by m_h/m_train,h (the
number of domain PSUs in stratum h over the number in the training fold) and
normalised to mean 1, so that each training fold represents the population.
{cmd:survey} is an option for sensitivity analysis when sampling might be
informative given the covariates. It is not covered by the theory, and
weighted fits can become unstable when the weights vary a lot. Targeting and
the variance are the same in both modes.

{phang}
{opt gbound(#)} truncates the out-of-fold propensity scores to [{it:#},
1-{it:#}] before targeting (TMLE) or before the pseudo-outcomes are formed
(AIPW); the default is 0.05, and {it:#} must be in (0, 0.5). The range of the
out-of-fold propensity score is stored in {cmd:e(ps_min)} and
{cmd:e(ps_max)}. Many values at the bound indicate limited overlap, and then
TMLE and AIPW estimates can differ noticeably.

{phang}
{opt level(#)}; see {helpb estimation options}.

{phang}
{opt rr} and {opt or} display the risk ratio and the odds ratio with
delta-method confidence intervals, computed on the log scale with the design
covariance of the two potential-outcome means and the design degrees of
freedom. They are always stored in {cmd:e()}; these options only control the
display. They can be requested again on replay: {cmd:svytmle, rr or}.

{pmore}
If a potential-outcome mean lies outside (0, 1) (within 1e-6 of a bound), the
ratios are not defined. This can happen with AIPW, whose means are not
bounded, or with a degenerate fit. The ratios are then stored as missing and
{cmd:e(rr_clamped)} is 1. {opt rr} and {opt or} cannot be used under
{cmd:mi estimate}, which treats them as its own options; see
{help svytmle##postestimation:Postestimation}.

{phang}
{opt generate(stub)} creates, for the analytic domain, the variables
{it:stub}{cmd:fold} (fold), {it:stub}{cmd:g} (out-of-fold propensity score,
truncated to [{it:gbound}, 1-{it:gbound}]), {it:stub}{cmd:Q0} and
{it:stub}{cmd:Q1} (out-of-fold outcome predictions under control and
treatment, truncated to [0.001, 0.999]), {it:stub}{cmd:Q0star} and
{it:stub}{cmd:Q1star} (targeted predictions; TMLE only), and
{it:stub}{cmd:eif} (the per-observation influence function of the ATE). The
variables are created only after a successful fit. {opt replace} overwrites
existing variables of exactly these names. With {cmd:estimator(aipw)},
{opt replace} also removes {it:stub}{cmd:Q0star} and {it:stub}{cmd:Q1star}
left behind by an earlier {cmd:svytmle} TMLE run.

{phang}
{opt noheader} suppresses the header above the coefficient table.


{marker method}{...}
{title:Method}

{pstd}
Let Y be the outcome, A the treatment, W the covariates, and w the survey
weight. In the analytic domain (n observations):

{phang2}
1. {it:Cross-fitting.} PSUs are assigned to V folds within each stratum (a
random permutation of the balanced pattern 1, ..., V, 1, ... within each
stratum). For each fold v, logistic regressions Q(A,W) = Pr(Y=1|A,W) and
g(W) = Pr(A=1|W) are fit on the other folds and predicted for fold v, giving
out-of-fold Q0, Q1, and g. The fits replicate R's {cmd:glm.fit}: IRLS with QR
least-squares steps, the same starting values and convergence rule, and a
25-iteration limit. Then g is truncated to [{it:gbound}, 1-{it:gbound}] and Q
to [0.001, 0.999].

{phang2}
2. {it:Targeting.} A two-parameter logistic fluctuation of logit Q_A on
H0 = (1-A)/g0 and H1 = A/g1 is fit with the survey weights as observation
weights. Here g1 = max(g, b) and g0 = max(1-g, b), with
b = 5/(sqrt(n) log n) (the R {cmd:tmle} package's default bound). This gives
targeted predictions Q0* and Q1*.

{phang2}
3. {it:Estimates.} The potential-outcome means and the ATE are the
survey-weighted (Hajek) means of Q1*, Q0*, and Q1*-Q0*.

{phang2}
4. {it:Variance.} The per-observation influence functions
D1 = A/g1 (Y-Q1*) + Q1* - psi1 and D0 = (1-A)/g0 (Y-Q0*) + Q0* - psi0 are
formed without weights. Their 2x2 design covariance is the stratified
first-stage with-replacement linearization variance of the domain Hajek
means. Any FPC and later sampling stages are ignored, which is the variance
that is calibrated for the super-population ATE. Inference uses t with the
design degrees of freedom (at least 1).

{pstd}
{it:AIPW.} With {cmd:estimator(aipw)}, step 2 is skipped. The
potential-outcome means are the survey-weighted (Hajek) means of the
pseudo-outcomes D1 = A/g (Y-Q1) + Q1 and D0 = (1-A)/(1-g) (Y-Q0) + Q0, which
use the truncated out-of-fold nuisances. Their centred values are the
influence functions in step 4.

{pstd}
Every stratum of the design must contain at least 2 PSUs. Post-stratification
or calibration declared in {cmd:svyset} ({cmd:poststrata()},
{cmd:postweight()}) is not supported. Weights that were calibrated before
release, as NHANES weights are, are used as fixed design weights, as
{cmd:svy} commands use them. The variance then ignores the calibration step.

{marker foldbalance}{...}
{title:Fold balance}

{pstd}
The large-sample theory of Karim (2026) assumes that the folds are balanced
within every stratum up to remainders, in two senses:

{phang2}
o each stratum's domain PSUs are split as evenly as possible;

{phang2}
o the total imbalance D = max_v sum_h |m_hv - m_h/V| is small relative to
sqrt(m), where m is the number of domain PSUs (formally, o(sqrt(m))).

{pstd}
This holds when only a few strata have a PSU count that V does not divide, or
when the number of strata stays fixed. It is not guaranteed when a fixed share
of many small strata is uneven. Exact balance (V divides every m_h, so D = 0)
is the clean case.

{pstd}
By default {cmd:svytmle} follows the R package in two ways:

{phang2}
o it reduces V to the smallest number of domain PSUs in any stratum;

{phang2}
o when V does not divide a stratum's PSU count, the leftover PSUs of that
stratum go to the lowest-numbered folds. For example, with {cmd:folds(2)} and
3 PSUs in every stratum, fold 1 holds two PSUs of every stratum.

{pstd}
{cmd:e(fold_exact)} records whether the folds were exactly balanced, and
{cmd:e(n_uneven)} counts the strata that were not split evenly.
{cmd:e(fold_psu)} gives the number of PSUs in each fold, and the header shows
the balance. With random folds and V = 2, D = {cmd:e(n_uneven)}/2. Pooled
NHANES 2007-2018 has 4 uneven strata out of 90, so D = 2. If a large share of
strata is uneven, choose a {opt folds()} value that divides the PSU counts. To
make exact balance a requirement, use {opt strictfolds}.

{pstd}
The variance treats the first stage as sampled with replacement (the
ultimate-cluster approximation). The theory of Karim (2026) is stated for that
idealized design. For first stages that select PSUs without replacement with
unequal probabilities, as NHANES does, the variance is an approximation, not a
proven bound, and exact fold balance does not change that.


{marker postestimation}{...}
{title:Postestimation}

{pstd}
{cmd:e(b)} holds {cmd:ATE}, {cmd:POmean1}, and {cmd:POmean0}, with their
design covariance in {cmd:e(V)} and the design degrees of freedom in
{cmd:e(df_r)}. {helpb lincom}, {helpb nlcom} (with {cmd:df(`e(df_r)')}),
{helpb test}, and {helpb estimates} work as usual. {cmd:test} reports an
unadjusted F(k, df) statistic, not {cmd:svy}'s adjusted Wald test.
{cmd:predict} and {cmd:margins} are not available; use {opt generate()}.
{helpb mi estimate} works directly ({cmd:svytmle} has the {cmd:mi}
property). With {opt seed()}, {cmd:mi estimate} resets the same seed in every
imputation, so all imputations use the same folds when the design and domain
are the same. The between-imputation variance then conditions on one split.
Omit {opt seed()} to draw different folds in each imputation, as the R
package's {cmd:svytmle_mi()} does. Do not specify {opt rr} or {opt or} under
{cmd:mi estimate}:
{cmd:mi estimate} treats them as its own options ({opt rr} is read as
{opt rrr}) and stops. Pool the relative effects as transformations instead,
as shown in the examples.


{marker examples}{...}
{title:Examples}

{pstd}{it:Public data.} NHANES II, one of Stata's example datasets ({cmd:webuse}
downloads it, so it needs internet), has 2 PSUs in each of 31 strata, so 2 folds
split every stratum evenly:{p_end}
{phang2}{cmd:. webuse nhanes2}{p_end}
{phang2}{cmd:. svytmle highbp, treatment(rural) covariates(c.age##c.age i.female i.race bmi i.region) folds(2) strictfolds seed(12345)}{p_end}

{pstd}The cross-fitted AIPW estimator on the same folds:{p_end}
{phang2}{cmd:. svytmle highbp, treatment(rural) covariates(c.age##c.age i.female i.race bmi i.region) folds(2) strictfolds seed(12345) estimator(aipw)}{p_end}

{pstd}A domain (adults aged 40 or older), with risk and odds ratios. The
variance uses the full design:{p_end}
{phang2}{cmd:. svytmle highbp, treatment(rural) covariates(c.age##c.age i.female i.race bmi i.region) subpop(if age >= 40) seed(1) rr or}{p_end}

{pstd}Fold balance, contrasts, and relative effects from the stored
potential-outcome means:{p_end}
{phang2}{cmd:. display e(fold_exact), e(n_uneven)}{p_end}
{phang2}{cmd:. matrix list e(fold_psu)}{p_end}
{phang2}{cmd:. lincom POmean1 - POmean0}{p_end}
{phang2}{cmd:. nlcom (logRR: ln(_b[POmean1]/_b[POmean0])), df(`e(df_r)')}{p_end}

{pstd}Survey-weighted nuisance models as a sensitivity analysis, saving the
per-observation values. The SE is the {cmd:svy} linearized SE of the mean of the
saved influence function:{p_end}
{phang2}{cmd:. svytmle highbp, treatment(rural) covariates(c.age##c.age i.female i.race bmi i.region) nuisance(survey) seed(1) generate(t_)}{p_end}
{phang2}{cmd:. svy: mean t_eif}{p_end}

{pstd}{it:NHANES 2007-2018.} The analytic files of Karim (2026) are not shipped
with the package. The public pipeline at
{browse "https://github.com/ehsanx/Survey-TMLE-codes"} builds them as R files
({it:E1_imputed.rds}, ...), and the repository's {cmd:cert/make_fixtures.R}
converts them to Stata ({it:cert/fixtures/N1_nhanes_E1.dta}, ...). Adults in the
analytic domain {cmd:inpop}, pooled MEC weights:{p_end}
{phang2}{cmd:. svyset SDMVPSU [pweight=WTMEC_POOLED], strata(SDMVSTRA)}{p_end}
{phang2}{cmd:. svytmle Y, treatment(A) covariates(age i.sex i.race i.educ pir i.cycle i.smk i.married i.pa phq) subpop(inpop) seed(20260607) rr or}{p_end}

{pstd}A richer outcome and propensity specification (quadratics and an
interaction):{p_end}
{phang2}{cmd:. svytmle Y, treatment(A) covariates(c.age##c.age i.sex##i.race i.educ c.pir##c.pir i.cycle i.smk i.married i.pa c.phq##c.phq) subpop(inpop) seed(20260607)}{p_end}

{pstd}In NHANES E1 a few strata have 3 PSUs, so 2 folds cannot be exactly
balanced. {cmd:strictfolds} reports this and stops:{p_end}
{phang2}{cmd:. svytmle Y, treatment(A) covariates(age i.sex pir) subpop(inpop) folds(2) strictfolds}{p_end}

{pstd}{it:Multiple imputation} (Rubin's rules), with a pooled log risk ratio. Set
the seed once and omit {opt seed()}, so that each imputation draws its own
folds and the run is still reproducible:{p_end}
{phang2}{cmd:. set seed 1}{p_end}
{phang2}{cmd:. mi estimate: svytmle Y, treatment(A) covariates(age i.sex pir) subpop(inpop)}{p_end}
{phang2}{cmd:. mi estimate (logRR: ln(_b[POmean1]/_b[POmean0])): svytmle Y, treatment(A) covariates(age i.sex pir) subpop(inpop)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:svytmle} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations in the design{p_end}
{synopt:{cmd:e(N_sub)}}number of observations in the analytic domain{p_end}
{synopt:{cmd:e(N_strata)}, {cmd:e(N_clust)}}strata and PSUs in the
design{p_end}
{synopt:{cmd:e(N_strata_sub)}, {cmd:e(N_clust_sub)}}strata and PSUs containing
domain observations{p_end}
{synopt:{cmd:e(N_strata_df)}, {cmd:e(N_clust_df)}}strata and PSUs counted in
the design df{p_end}
{synopt:{cmd:e(N_pop)}, {cmd:e(N_subpop)}}estimated population and domain
sizes{p_end}
{synopt:{cmd:e(df_r)}}degrees of freedom used for inference,
max(1, {cmd:e(df_design)}){p_end}
{synopt:{cmd:e(df_design)}}design degrees of freedom,
{cmd:e(N_clust_df)} - {cmd:e(N_strata_df)}{p_end}
{synopt:{cmd:e(folds)}, {cmd:e(V_eff)}}requested and effective number of
folds{p_end}
{synopt:{cmd:e(fold_exact)}}1 if the folds were exactly balanced within
strata{p_end}
{synopt:{cmd:e(n_uneven)}}strata whose domain PSUs were not split
evenly{p_end}
{synopt:{cmd:e(strictfolds)}}1 if {opt strictfolds} was specified{p_end}
{synopt:{cmd:e(gbound)}}propensity truncation bound{p_end}
{synopt:{cmd:e(gbound_tmle)}}targeting lower bound 5/(sqrt(n) log n) (TMLE){p_end}
{synopt:{cmd:e(n_gbound_tmle)}}observations where that bound binds (TMLE){p_end}
{synopt:{cmd:e(eps0)}, {cmd:e(eps1)}}fluctuation parameters (TMLE){p_end}
{synopt:{cmd:e(ps_min)}, {cmd:e(ps_max)}}range of the out-of-fold propensity
score{p_end}
{synopt:{cmd:e(n_nonconv)}}logistic fits that did not converge{p_end}
{synopt:{cmd:e(rr)}, {cmd:e(se_logrr)}}risk ratio and the SE of its log{p_end}
{synopt:{cmd:e(or)}, {cmd:e(se_logor)}}odds ratio and the SE of its log{p_end}
{synopt:{cmd:e(rr_clamped)}}1 if a potential-outcome mean was outside (0, 1),
so that the ratios are missing{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:svytmle}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(cmdversion)}}version of {cmd:svytmle} that produced the results{p_end}
{synopt:{cmd:e(depvar)}, {cmd:e(treatment)}}outcome and treatment{p_end}
{synopt:{cmd:e(covariates)}}covariates as specified{p_end}
{synopt:{cmd:e(xterms)}}expanded covariate terms (base levels excluded){p_end}
{synopt:{cmd:e(estimator)}}{cmd:tmle} or {cmd:aipw}{p_end}
{synopt:{cmd:e(learner)}}{cmd:logit}{p_end}
{synopt:{cmd:e(nuisance)}}{cmd:unweighted} or {cmd:survey}{p_end}
{synopt:{cmd:e(subpop)}}domain specification{p_end}
{synopt:{cmd:e(wvar)}, {cmd:e(su1)}, {cmd:e(strata1)}}design variables
used{p_end}
{synopt:{cmd:e(seed)}, {cmd:e(foldvar)}}fold source{p_end}
{synopt:{cmd:e(rngstate)}, {cmd:e(rng)}}random-number state and generator
before the folds were drawn (not stored with {opt foldvar()}){p_end}
{synopt:{cmd:e(generate)}}variables created by {opt generate()}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Linearized}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}{cmd:ATE}, {cmd:POmean1}, {cmd:POmean0}{p_end}
{synopt:{cmd:e(V)}}design covariance (singular by construction:
ATE = POmean1 - POmean0){p_end}
{synopt:{cmd:e(fold_psu)}}number of PSUs in each fold{p_end}

{p2col 5 22 26 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the observations in the design{p_end}


{marker validation}{...}
{title:Validation}

{pstd}
The package's source repository, {browse "https://github.com/ehsanx/svytmle-stata"}, contains a
certification suite in its {cmd:cert/} folder. An R script produces reference
results with the R package {bf:svytmle} 0.1.2 ({cmd:svytmle()} and
{cmd:svyaipw()}), including the exact folds R drew for each case. The cases
are simulated designs and the four NHANES analyses of Karim (2026). A Stata
do-file re-estimates every case with {opt foldvar()} and compares the
estimates, standard errors, df, confidence intervals, relative effects,
targeting parameters, and the per-observation propensity scores and influence
functions. It also checks invariances and error handling.

{pstd}
The repository ships the simulated cases, so {cmd:do cert/run_cert.do} reruns
their Stata-side checks without R. The NHANES analytic files are built by the
public pipeline at {browse "https://github.com/ehsanx/Survey-TMLE-codes"}.

{pstd}
The simulated designs include:

{phang2}o nested PSU ids;{p_end}
{phang2}o domains that leave PSUs empty;{p_end}
{phang2}o unstratified and observation-level designs;{p_end}
{phang2}o n = 80;{p_end}
{phang2}o 2 PSUs per stratum;{p_end}
{phang2}o a training fold that is rank-deficient;{p_end}
{phang2}o zero weights;{p_end}
{phang2}o complete separation in the outcome and in the propensity model;{p_end}
{phang2}o survey-weighted nuisances.{p_end}

{pstd}
Stata and R agree to about 1e-12 for both estimators, with two exceptions.
Under complete separation, where R's {cmd:glm.fit} itself does not converge,
they agree to about 1e-10. For a covariate whose mean is a million times its
standard deviation, they agree to about 2.5e-9.


{marker references}{...}
{title:References}

{phang}
Karim, M. E. 2026. Cross-fitted survey-weighted TMLE with design-based
variance for causal machine learning. arXiv:2606.30918.
{browse "https://arxiv.org/abs/2606.30918"}. R package {bf:svytmle}:
{browse "https://github.com/ehsanx/svytmle"} (version 0.1.2, commit 5534234).

{phang}
Binder, D. A. 1983. On the variances of asymptotically normal estimators from
complex surveys. {it:International Statistical Review} 51: 279-292.

{phang}
Chernozhukov, V., D. Chetverikov, M. Demirer, E. Duflo, C. Hansen, W. Newey, and
J. Robins. 2018. Double/debiased machine learning for treatment and structural
parameters. {it:Econometrics Journal} 21: C1-C68.

{phang}
van der Laan, M. J., and S. Rose. 2011. {it:Targeted Learning: Causal Inference}
{it:for Observational and Experimental Data}. New York: Springer.

{phang}
Data: National Center for Health Statistics. National Health and Nutrition
Examination Survey, 2007-2018. Hyattsville, MD: Centers for Disease Control and
Prevention. {browse "https://www.cdc.gov/nchs/nhanes/"}.


{title:Author}

{pstd}
M. Ehsan Karim, School of Population and Public Health, University of British
Columbia, Vancouver, Canada. {browse "mailto:ehsan.karim@ubc.ca":ehsan.karim@ubc.ca}
{p_end}

{pstd}
Source, certification suite, and issue tracker:
{browse "https://github.com/ehsanx/svytmle-stata"}. Please report problems at
{browse "https://github.com/ehsanx/svytmle-stata/issues"}, with the output of {cmd:which svytmle}
and your Stata version.
{p_end}

{pstd}
Requires Stata 15 or newer; tested on Stata 15.1. License: GPL-3.0-or-later
({browse "https://github.com/ehsanx/svytmle-stata/blob/main/LICENSE.md":LICENSE.md}).
{p_end}
