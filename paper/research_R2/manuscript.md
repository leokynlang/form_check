# Silent Acceptance: A Benchmark of Specification Gates and Extrapolation Failures in Physical Model Assessment

**Revision status - R2, 2026-10-18. Working manuscript.**
Research presentation revised using the R1 audited results.
Appendix G records the correction history and the scope of this revision.

## Abstract

Specification tests are used to challenge candidate equations, but their
value depends on which departures they detect and what errors remain
after a model passes them. We benchmark residual-permutation curvature,
RESET variants, and nested-quadratic comparisons on prespecified
power-law candidates, measuring false alarms, detection, declined
evaluations, and extrapolation-related joint release.

On a common pool of 18 non-monomial physics tasks, regressor-power
HC3 RESET has the highest mean detection fraction among the four
specification tests at each evaluated sample size under the
H0-reference 5% calibration target. Its paired task-level advantage
over permutation curvature is supported at n=200 and 1000, while
constructed odd-symmetric, opposing-curvature, and interaction
alternatives reveal distinct diagnostic blind spots. Native false-positive
rates also vary with design, sample size, and numerical implementation.

In the 17-task pool with finite extrapolation outcomes, the complete
workflow combines R² ≥ 0.99 with four valid, non-rejecting diagnostics.
It releases 481 wrong models whose relative error exceeds 10% at the
registered upper-end probe: 7.5% of screened dangerous wrong models
(481/6,403), or 0.73% of all screened wrong models (481/65,992).
The releases are concentrated in two of the seventeen tasks; excluding
the threshold-borderline task bonus.15, whose median probe error is
9.2% against the 10% cutoff, the release fraction is 0% at n=200 and
1000. The release fraction among all dangerous eligible fits decreases
from 2.94% at n=50 to 0.21% at n=1000. Thus the gates remove many
dangerous fits while leaving a measurable residual failure channel
whose size depends on task composition and threshold placement.

The benchmark connects diagnostic performance to complete-workflow
outcomes and identifies which mismatch directions warrant additional
checks. Its results apply to the stated candidate family, task pools,
and decision rules; independent joint-policy calibration and evaluation
over a fixed extrapolation domain are the next validation steps.

## 1 Introduction

A candidate physical equation is useful only to the extent that the
available evidence supports its intended use. Prediction error describes
agreement on sampled inputs; specification tests ask whether selected
departures from the assumed form remain detectable. An engineering
workflow often combines the two: fit a candidate, require a high fit
score, run diagnostics, and retain models that trigger no alarm. The
practical question is how much protection that decision rule provides,
and which errors it still releases.

The distinction is concrete in a controlled friction example. A power
law fitted to Prandtl-Karman data reaches R² = 0.99, passes
cross-validation and an information-criterion comparison, yet incurs
about -20% error one decade beyond training. In the corresponding
synthetic scenario, the four specification tests studied here detect the
departure at power ≥ 0.99. This motivates evaluating the diagnostics
alongside fit quality: a fit-only failure can be detected, but other
mismatch directions may survive the same checks.

We study this behavior for prespecified power-law candidates.
Permutation curvature compresses departure from log-linearity into one
quadratic direction; RESET adds powers of fitted values or regressors;
nested-quadratic comparisons evaluate a richer rival using
cross-validation or AICc. These choices create testable differences in
sensitivity. False alarms and numerical stability also matter: an alarm
that detects more wrong forms may reject more correct forms, and an
unstable inference implementation can distort both sides of the
comparison.

The benchmark combines controlled synthetic alternatives, physics
equations sampled from SRSD configurations, family-stratified tasks,
a family-viability battery, nested-rival comparisons, and observational
compatibility probes. Ground truth supplies form labels in the synthetic
and physics-derived suites. Repeated datasets measure sampling
variability within those tasks; equations provide the units for the
task-level comparisons. Decision thresholds, eligible populations, and
declined evaluations are reported with each result.

Three findings organize the study. First, diagnostic sensitivity follows
the mismatch direction. Opposing curvature and odd-symmetric departures
can be weak in a mean-square projection, while interactions expose
differences between fitted-value and regressor-power augmentation.
On the common 18-task physics pool, SVD-based regressor-power HC3
has the largest mean detection fraction among A-D at all three sample
sizes, with positive paired D-A intervals at n=200 and 1000.
Second, detecting residual structure and controlling false alarms are
separate empirical requirements. Native size differs across designs;
a stable numerical path is necessary for a meaningful comparison.
Third, complete-workflow release is more informative than gate silence
alone. In the 17-task extrapolation pool, 481 dangerous wrong models
pass both R² ≥ 0.99 and every valid diagnostic: 7.5% of the
screened-dangerous subset. These releases are concentrated in two
tasks; excluding the threshold-borderline task bonus.15, the release
fraction is 0% at n=200 and 1000, so the pooled 7.5% is not a
task-uniform failure rate. The resulting failure channel can be
quantified alongside the errors the gates remove.

We make three contributions:

1. **A diagnostic benchmark with explicit decision accounting.**
   The six core suites provide ground-truth form labels, specified
   sampling configurations, threshold policies, and typed failure
   registration. We distinguish simulation precision within tasks
   from variation across tasks.

2. **A map of detection, false alarms, and mismatch directions.**
   Controlled alternatives connect method rankings to specific blind
   spots. Family viability and nested-rival comparisons show how a
   diagnostic can support subsequent checks, while implementation
   cross-checks expose a numerical failure that changes the empirical
   HC3 comparison.

3. **An evaluation of extrapolation-related joint release.**
   We connect the fit screen, valid gate decisions, and a registered
   far-field probe on the same fitted models. Explicit conditional
   denominators separate the dangerous models released among screened
   candidates from the fraction of all dangerous fits that survive.

Existing symbolic-regression benchmarks evaluate predictive quality,
expression complexity, and symbolic or structural recovery [5, 6].
Our experimental object is the behavior of attached specification gates
under a fixed candidate family. This focus complements recovery metrics:
it measures the evidence supplied by a diagnostic record when the true
form is unavailable to the analyst. Validation on expressions chosen by
a discovery algorithm is a further step, developed as a research
priority in Section 6.

Sections 3 and 4 define the protocol and report the results. Section 5
draws practical lessons about selecting diagnostics and interpreting
release. Section 6 identifies the experiments needed for stronger
workflow recommendations. Provenance and correction details are
collected in Section 7 and Appendix G.

## 2 Related Work

**Specification testing and robust inference.** Ramsey's RESET [1]
tests powers of fitted values; regressor-power variants examine other
directions of misspecification. Heteroskedasticity-consistent covariance
estimation [2] and HC3 leverage adjustment [3, 4] address inference under
nonconstant error variance. We compare these established tools on
log-transformed, physics-derived designs. The numerical audit in
Section 5.2 and Appendix G concerns implementation of the HC3-Wald
calculation; the statistical correction and its numerical realization
must be distinguished.

**Model comparison as a diagnostic.** Information criteria [12] and
cross-validation [13] can compare a candidate with a richer rival.
Our S5 comparison freezes a quadratic rival and reports both native
and H0-reference-calibrated policies. This establishes a practical
alternative to the single-direction curvature statistic. Its comparison
with A is interpreted using the measured error budgets in Section 4.7.

**Resampling and transformed coordinates.** Wild-bootstrap methods
[10, 11] provide tools for inference under heteroskedasticity.
Residual permutation relies on exchangeability-related conditions.
The family-viability experiment examines how these choices behave after
linearizing transformations: replacing permutation with a wild null
substantially reduces the observed rational-coordinate distortion.
The reported calibration evidence is empirical for the tested designs.

**Symbolic-regression evaluation.** SRBench [5] evaluates predictive
accuracy, expression complexity, and symbolic recovery. SRSD [6]
revises physically motivated sampling configurations and evaluates
structural discrepancy using normalized edit distance. We reuse SRSD
configurations to assess diagnostic decisions for prespecified power
laws. The measurements concern detection, false alarms, and
extrapolation-related release rather than the recovery rate of a
search algorithm.

**Physics-informed discovery.** AI Feynman [7] uses physical units
to constrain search; measurement-driven approaches recover governing
structure from observations [8]; PySR [9] supplies a symbolic-search
framework. Dimensional constraints and specification diagnostics
provide different information. A dimensionally valid non-monomial
law still misspecifies a power-law candidate. Our benchmark examines
that form mismatch after fixing the candidate family.

**Position of the benchmark.** The study contributes a common
experimental record linking selected diagnostic gates to mismatch
directions and complete-workflow outcomes. Its comparative scope
covers permutation curvature, fitted-value and regressor-power RESET,
HC3 inference, and one nested-quadratic rival. Other residual tests,
rival families, and post-selection procedures can be evaluated within
the same accounting framework. The experiments establish findings
for the tested gates and sampling designs; the wider specification
testing literature remains the context for interpreting them.

## 3 Benchmark Design

### 3.1 Design principles

**P1 - Explicit decision policies.** Each comparison identifies
the statistic, threshold, eligible population, and handling of declined
evaluations. Native p-value rules use alpha=.05. The permutation
curvature p-value uses B=999 residual permutations; RESET uses an
F reference. S5 compares native selection thresholds with
H0-reference-calibrated thresholds. These policies are labeled
separately in the results.

**P2 - A common reference target for headline power.** Table 3
calibrates each method to the empirical 5% quantile of its finite H0
p-values within the same sample-size stratum. Section 4.2 states
the inclusive quantile convention and the resulting scope of the
comparison. Native-size results, the mixed policies of S5, and joint
panels retain their measured error rates. Independent evaluation of
the final joint policy is a validation priority in Section 6.

**P3 - Structure-defined labels.** H0 denotes a strict monomial
whose mean is exactly linear in log-log coordinates; H1 denotes
a non-monomial ground-truth structure. The symbolic registry supplies
these labels rather than diagnostic outcomes. S3 family labels and
S4 truth families likewise derive from their generators. Operational
eligibility additionally requires valid data and the recorded
consistency checks.

**P4 - Explicit failure registration.** Invalid evaluations retain
a typed reason: zero-residual degeneracy, rank deficiency, insufficient
degrees of freedom, leverage or covariance degeneracy, statistic
mismatch, generation failure, or unavailable observational data.
Each table states its evaluated denominator. A declined diagnostic
cannot produce a rejection or contribute to joint release. Zero-noise
cells exercise the degeneracy path.

**P5 - Source sampling configurations.** S2/S3 reuse the source
registry's per-variable sampling objects. Tasks incompatible with the
positive-input/output requirement are marked not applicable. The
observational suite uses World Bank Open Data [14], with its indicators
and cached snapshots recorded in the artifact.

### 3.2 Gates under test

All gates operate in log-log space on positive data: let ly = log y and
X = [1, lx₁, …, lxₖ] be the design of the candidate power-law model.

**A. Permutation curvature (form_check).** Adds a single mean-centered
quadratic column q = meanⱼ((lxⱼ − lx̄ⱼ)²) and takes the added coefficient
c₂ as statistic; size evidence comes from B = 999 residual permutations
(vectorized through the pseudo-inverse row), with a statistic-mismatch
sentinel comparing wrapper and native values on every call. Validity
note: permuting *estimated* residuals assumes they are exchangeable,
which OLS residuals are not exactly - the permutation null here is an
approximation whose quality is an empirical question, not a finite-sample
exactness theorem. The benchmark measures this construction's empirical
size on the specified H0 pools, including deviations from .05 (Tables 1-2).
Log-transformed multiplicative noise and favorable designs may matter, and §6 lists the stress
tests (leverage-adjusted residuals, wild and block nulls, cluster
dependence, fixed-design variants) that would bound where the
approximation breaks.

**B. RESET fitted, classical F.** Adds ŷ² and ŷ³; classical F test of
joint significance.

**C. RESET fitted, HC3-Wald.** Same augmented design; HC3 sandwich Wald
statistic on the added coefficients, referred to F(q, n − k).

**D. RESET exog, HC3-Wald.** Built regressors are powers of each
log-regressor; HC3-Wald as in C.

**E. CV-quad (S5).** The "obvious alternative": rival = full quadratic
expansion of the log features (linear terms, per-feature squares, pairwise
cross-terms; minimum-norm least squares absorbs rank deficiency).
Statistic D = RMSE_cv(power) − RMSE_cv(quad) over 5 deterministic folds;
escalate if D > δ. Naive δ = 0; calibrated δ\* chosen on the H0 pool to
FPR = 5%.

**F. AICc-quad (S5).** Same rival; escalate if AICc(power) − AICc(quad) >
m, k = column rank. Naive m = 0; calibrated m\* as in E.

Configurations were frozen before all main runs - powers (2, 3); B = 999;
α = 0.05; the rival family and criteria of E/F; the grids below - and are
identical across suites. The same harness invokes all gates on every
dataset, so comparisons share data, seeds, and failure handling exactly.

### 3.3 Task suites S1-S3

**S1 - Synthetic control (19 scenarios × 500 seeds).** Three monomial H0
scenarios (univariate; two-variable; two-variable with ρ ≈ 0.9 correlated
design). Eleven structural H1 scenarios: quadratic curvature (c = 0.10);
cubic (standardized odd term); multiplicative interaction (log-space
cross term); threshold (exponent 0.4 → 1.0 at the domain midpoint);
exponential law; logarithmic law; rational (saturating) law; omitted
variable; dilution (one of two variables curved); opposing curvature
(+0.12, −0.12 on an exchangeable design, so the mean-square column's
population coefficient is exactly zero); collinear-quadratic (ρ ≈ 0.95
design with curvature). Five noise-layer H1 scenarios have strict-monomial
mean structure but non-Gaussian noise - Student-t(3); SD growing with the
log prediction (0.3× to 1.7×); 2% point outliers (+8σ); 6% group outliers
(+5σ); AR(1) ρ = 0.6 - rejections there are false alarms by construction
and are never pooled into structural power. Grid: n ∈ {50, 200, 1000};
relative Gaussian noise σ ∈ {0, 0.5%, 1%, 3%, 5%} in log space; 500
seeds per cell via a seed rule embedding the grid position.

**S2 - SRSD-Feynman tasks (32 originally fitted tasks, 500 seeds).**
The original positive-input screen selected 13 monomial H0 and
19 non-monomial H1 tasks from the 120-equation registry. Each fit
has a 2,000-sample same-distribution test set. The extrapolation
code constructs 12 log-spaced probe points along the first feature:
six below the observed training minimum (down to minimum/3) and
six above the observed maximum (up to 3 times maximum); all other
log-features are held at their training medians.
The cached far-field metric used in this paper is the relative
error at the last, upper-end point only, against the noiseless
closed-form truth. It is not a one-decade probe or the maximum
error over the probe grid. Its coordinate depends on the observed
training maximum and can change with n and seed.
The generated fitted-run count is 240,000; downstream
eligibility excludes i.30.5 after its sentinel failure and excludes
i.26.2 only when a finite extrapolation outcome is required.

**S3 - Family-stratified Feynman equations (30 equations × 30 seeds).**
Structural family classification from symbolic form: exact power law (12
H0); fractional power (1 H0); exponential (3 H1); rational (6 H1);
trigonometric (5 H1); plus not-applicable members with explicit reasons.
Each H1 equation runs over relative noise {0.5%, 1%, 3%, 5%} ×
{homoskedastic, heteroskedastic} × sampling-range scale {official, ×0.5,
×2}; H0 equations run the noise grid at official range.

### 3.4 Suite S4 - the family-viability battery

The battery evaluates four families through their linearizing
coordinates: power (log x, log y), exponential (x, log y),
logarithmic (log x, y), and rational (1/x, 1/y).
Each member applies A's quadratic construction in its own
coordinates and returns a viable set V={F: p_F ≥ alpha}.
Columns are standardized before squaring. The formal run uses
a wild-bootstrap null with Rademacher signs and B=999.
Section 4.6 compares this choice with residual permutation.

The set-level outcomes are false exclusion of a true family,
inclusion of a wrong family, and escalation when no family
survives. Escalation is correct when the truth is outside all
four families. The surviving set supplies candidates for further
comparison rather than a unique selected model.

The grid has n=200, relative noise from 0.5% to 5%,
500 seeds per cell, and 188 cells (94,000 datasets).
It includes the 19 S1 scenarios with generator-defined
truth-family labels and 14 SRSD H1 equations under
homoskedastic and heteroskedastic noise.

### 3.5 Suite S5 - the gate comparison

S5 compares permutation curvature with the full-quadratic rivals
E (cross-validation) and F (AICc), as specified in gate_compare.py.
Its pools comprise the 19 S1 scenarios and the 32 originally
fitted S2 tasks, with n in {50, 200}, rel=.03, and
200 seeds per cell. A, E, and F evaluate identical datasets.

A uses its native permutation rule. E/F are evaluated at both
native selection thresholds and per-pool, per-n H0-reference
thresholds. Tables 11-12 report measured size and detection
under these policies. Because the resulting budgets differ,
the comparison is interpreted using the reported H0 rates.

Reported outcomes include structural-H1 detection, noise-layer
false alarms, task-level detail, and computational cost.
The S7 reanalysis uses the same records to compare A alone
with the A/E/F panel and reports the accompanying change
in joint H0 rejection.

### 3.6 Suite S6 - real data (World Bank open macro data)

Synthetic and equation-derived suites have ground truth by construction;
the deployment question is whether the synthetic failure-mode map predicts
anything on observational data, where no truth exists. S6 is deliberately
small, prespecified, and divergence-oriented: six macroeconomic relations
(target, features) × four years (2019, 2021, 2022, 2023), cross-country
cross-sections of 184-199 sovereign economies (World Bank aggregates
excluded), from eight indicators pulled from the World Bank API with raw
responses cached for reproducibility - life expectancy ~ GDP per capita
(PPP); infant mortality ~ GDP per capita; GDP ~ population; CO₂ per
capita ~ GDP per capita; internet penetration ~ GDP per capita; and life
expectancy ~ GDP per capita + health expenditure per capita (two-feature).
Gates A-D run at 5 permutation seeds per cell. These seeds are
computational replicates of one p-value estimator on one fixed dataset -
not five independent samples - so they are used only to report Monte
Carlo stability (the empirical p has resolution 1/1000 at B = 999; across
the 24 S6 cells the median p spread is .009 and the max .094), never as a
sample-size or stability claim; a prespecified probe adds three agricultural
and energy relations over three years, including a built-in
control-contrast pair - cereal yield ~ fertilizer (one feature) versus
cereal yield ~ fertilizer + energy use per capita (two features) - chosen
to test, on live data, whether the dilution mechanism of S1 (a second
correlated curved feature halving the mean-square column's signal)
operates outside the laboratory. No calibration or power claim is made on
S6; its contract is external validity of *mechanisms*, not of rates. An
engineering companion bench (Appendix E) extends the same contract to
physical-engineering data: nine open sources, 91 datasets, 150
prespecified relations, run as a single-realization smoke pass.

### 3.7 Suite S7 - joint-silence reanalysis of existing results

S7 evaluates panel non-rejection before imposing the fit screen.
The A-D arm joins S2 records to SVD C/D outputs one-to-one on
(task, n, rel, seed), requires positive noise and a valid statistic
sentinel, and contains 18 H1 tasks and 13 H0 tasks. At each
n in {50, 200, 1000}, a task contributes four noise levels times
500 seeds. Native alpha=.05 and per-gate H0-reference 5% thresholds
define separate policies. Joint silence requires finite p-values for
every gate and no rejection.

The A/E/F arm uses the S5 SRSD pool at rel=.03: 19 H1 tasks
and 13 H0 tasks, with 200 seeds each at n=50 and 200.
A retains its native rule; E/F use their individual H0-reference
thresholds. A alone and A/E/F are compared on the same records,
with the change in joint H0 rejection reported alongside detection.

The two arms have different task pools, noise levels, seeds, and
gate sets. They estimate their own panel behavior; only the paired
within-arm comparison measures the effect of adding E/F. Both
arms reuse existing fitted records. Reference thresholds and the
reported H0 rates use the same pools, so the evidence describes
reference-pool calibration. audit_s7.py records filters and exact counts.

Section 4.10 additionally requires finite fit-quality and far-field
outputs, yielding 17 H1 tasks. This defines the joint-release pool;
S7's gate-only analysis uses all 18 eligible H1 tasks.

## 4 Results

Tables 1-2 report size; Tables 3-4 detection; Table 5 the family
map; Table 6 noise-layer behavior; Table 7 extrapolation; Tables
8-10 the battery; Tables 11-12 the model comparison; Tables
13-14 observational examples; Table 15 joint release; and Table 16
joint silence. Each table defines its event, denominator, task
pool, and failure policy. Computational repeats are not independent
equations. S7 reuses S2/S5 records rather than adding new fits.
Sections 4.2, 4.10, and 4.11 apply the eligibility rules in
Section 3 and the accompanying record-level audit scripts.

### 4.1 Size: false-positive rates on true monomials

**Table 1 - S1 H0 false-positive rates at p < .05. Three scenarios, four positive noise levels, and 500 seeds give 6,000 valid evaluations per n. Entries are rates with fit-level Wilson 95% intervals.**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .0583 [.0527-.0646] | .0457 [.0407-.0512] | .0513 [.0460-.0572] |
| B | .0557 [.0501-.0618] | .0518 [.0465-.0577] | .0472 [.0421-.0528] |
| C (SVD) | .0602 [.0544-.0665] | .0568 [.0513-.0630] | .0478 [.0427-.0535] |
| D (SVD) | .0675 [.0614-.0741] | .0525 [.0471-.0584] | .0502 [.0449-.0560] |

The evaluated-H0 denominator excludes 1,500 zero-noise records
per n, each registered as degenerate_zero_residual. Such a record
has no valid test decision and is not a successful non-rejection.

**Table 2 - S2 H0 false-positive rates at p < .05. Thirteen equations, four positive noise levels, and 500 seeds give 26,000 valid evaluations per n; C/D use SVD.**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .0489 [.0463-.0516] | .0535 [.0508-.0563] | .0517 [.0491-.0545] |
| B | .0499 [.0473-.0526] | .0497 [.0471-.0524] | .0527 [.0500-.0555] |
| C (SVD) | .0565 [.0538-.0594] | .0531 [.0505-.0559] | .0549 [.0522-.0577] |
| D (SVD) | .0604 [.0576-.0633] | .0564 [.0536-.0593] | .0583 [.0555-.0612] |

Both tables are generated from existing positive-noise H0 records
by audit_size.py. Stored rejection flags agree with strict p < .05;
ties at exactly .05 are not rejections. This differs from the
inclusive empirical-quantile convention used in Table 3.

On synthetic Gaussian data, Table 1 shows small-sample
size distortion: at n=50 all four fit-level Wilson intervals lie
above .05, while the pattern is milder at larger n. Calibration is
therefore empirical and imperfect, not exact. On the physics-derived designs (Table 2), A and B are
near .05, while C and D have a modest upward size distortion:
all six displayed C/D fit-level Wilson intervals exclude .05.
The SVD calculation underlying these values was cross-checked
against a separate implementation (Section 5.2).

Generalization over physical equations requires an equation-level
uncertainty analysis, distinct from a fit-level interval on this
fixed task pool. At n=200, resampling the 13 H0 equations gives
C=.053 with a 95% interval [.050, .056] and D=.056 with
[.055, .058]; the corresponding leave-one-equation-out ranges
are [.052, .054] and [.056, .057]. The C lower bound rounds
to .050 but is slightly above .05 before rounding. Appendix C
uses this same SVD result and n. These intervals describe the
observed task pool, not all possible physical laws.

A separate fitted-HC3 implementation check found 8/120 rejections
(.067, Wilson [.034, .126]) on two tasks. This small-sample
check addresses numerical agreement; Table 2 and the equation-level
analysis supply the size estimates for the S2 pool. Appendix G
records the original-path comparison and its interpretation.

### 4.2 Power at a common reference-pool calibration target

Table 3 uses the same 18 sentinel-valid H1 tasks for every method;
feynman-i.30.5 is excluded. Each n has 36,000 records: four
positive noise levels times 500 seeds per task. Thresholds are
the empirical 5% quantiles of each method's finite H0-reference
p-values at the same n. Table 3 rejects when a valid p-value is
less than or equal to that quantile; native decisions instead use
strict p < .05. Non-finite p-values and typed declines are never
rejections or certificates. Quantile discretization and boundary
ties mean that this is a calibration target, not a claim of
exactly equal rejection counts or independent calibration.

**Table 3 - S2 rejection fraction on the common 18-H1-task pool at the H0-reference 5% target; C/D use SVD. Declines remain in the denominator and are not counted as rejections.**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .386 | .486 | .570 |
| B | .187 | .254 | .383 |
| C (SVD) | .128 | .221 | .383 |
| D (SVD) | **.466** | **.657** | **.718** |

For B and C, respectively, 2, 56, and 591 p-values are non-finite
out of 36,000 at n=50, 200, and 1000. These records cannot
contribute a rejection or a silent certificate. The displayed
fractions are unconditional detection yields on the specified
records, not power conditional on successful evaluation.

D has the largest mean rejection fraction among A-D in this
pool. The paired equation-level bootstrap compares D minus A
on these same 18 tasks (10,000 resamples): +.080
[-.083, +.256] at n=50, +.170 [+.042, +.329] at n=200,
and +.148 [+.022, +.302] at n=1000. The first interval crosses
zero; the latter two favor D in this task pool. D exceeds A on
10/18, 13/18, and 8/18 tasks, with 0, 1, and 5 ties,
respectively. A larger pooled mean is therefore not uniform
task-wise dominance.

Table 3, its paired effects, and Figures D.2-D.3 use the same
eligibility rule. S5 evaluates the nested-rival gates on a
separate pool; a paired D-versus-F comparison is a validation
priority. The constructed alternatives below illustrate the
method-specific differences behind the pooled ranking.

**Table 4 - S1 blind-spot cartography (raw α, n = 200, rel = 1%).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| cubic (standardized odd term) | .40 | 1.00 | 1.00 | 1.00 |
| interaction (lx₁·lx₂ cross term) | .21 | 1.00 | 1.00 | .10 |
| opposing curvature (+c, −c, exchangeable) | **.02** | .98 | .98 | 1.00 |
| dilution (one curved variable of two) | 1.00 | .88 | .87 | 1.00 |

Three mechanisms, each verified by construction: (i) A's statistic is the
mean of squared centered columns - under opposing curvature on an
exchangeable design that mean has exactly zero population signal; (ii) a
quadratic-projection statistic loses power on odd-symmetric deviations;
(iii) B's fitted-value powers are a local Taylor surrogate for any smooth
alternative - saturating on interaction where A (.21) and D (.10)
struggle - and, as S3 shows, the same surrogate property near-blinds B on
rational laws.

### 4.3 Family-stratified results and sensitivity

**Table 5 - S3 calibrated power by family (C/D on the verified SVD path).**

| Family | A | B | C (SVD) | D (SVD) |
|---|---|---|---|---|
| exponential | .09 | .06 | .07 | .06 |
| rational | **.65** | .16 | .13 | **.67** |
| trigonometric | .66 | .69 | .58 | **.86** |

B - near-perfect on synthetic quadratics - is nearly blind to rational
laws (.16 calibrated): the fitted-value Taylor surrogate absorbs
precisely the low-order curvature by which a saturating law departs
from a power law. Under the SVD path D is the strongest gate on both
rational (.67) and trigonometric (.86) families; the previous draft's
near-zero C/D entries on trigonometric were the pinv artifact
(thresholds derived from artifact p-mass), and the S3-level
correction agrees in direction with the S2-level one of Table 3.
Exponentials defeat every gate (≤.09 - the linearization a log
transform induces on them is nearly orthogonal to all four
augmentation families); A and D share the rational lead; no method
dominates all three families. Every method's strength sits on
another method's blind spot; Table 5 has no row-major winner.
Sensitivity: scaling every sampling domain ×2 moves exactly one
family-method cell by more than 15 points (rational × A: .66 → .84, the
expected direction); ×0.5 moves none; heteroskedastic replications agree
within .02 on every cell. **These sensitivity readings rest on 30 seeds
per cell and cannot be separated from sampling noise at that size**: a
single cell at p ≈ .66 carries a 95% interval of roughly ±.17, and the
standard error of a difference between two independent 30-seed cells is
≈ .11, so a 15-point trigger cannot distinguish movement from
variation. Family × method intervals are now attached
(`paper/review3_analysis.py` §3.8b). On the within-cell binomial
convention the trigonometric D (SVD) cell is .860 [.848, .871] and
rational D is .669 [.655, .683]; on the honest unit - the equation -
the same cells read .860 [.592, .998] and .669 [.409, .864]. Rational C
moves from .132 [.122, .142] to [.053, .277] the same way. The
equation-level intervals overlap across families and are wide enough to
contain orderings the point estimates contradict, so the family ranking
above is directional, not an ordered effect list; what survives the
wider intervals is the sign pattern - D ≥ C on rational and
trigonometric, every gate low on exponential - not the point ordering.

### 4.4 Noise-layer robustness

**Table 6 - S1 noise-layer H1 (strict-monomial mean structure) rejection rates at α = 0.05 - false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds; C/D on the verified SVD path).**

| Scenario | A | B | C (SVD) | D (SVD) |
|---|---|---|---|---|
| Student-t(3) noise | .034 [.021-.054] | .036 [.023-.056] | .030 [.018-.049] | .030 [.018-.049] |
| variance growing with prediction | .048 [.032-.070] | .058 [.041-.082] | .050 [.034-.073] | .050 [.034-.073] |
| AR(1) errors, ρ = 0.6 | .038 [.024-.059] | .040 [.026-.061] | .042 [.028-.063] | .042 [.028-.063] |
| 2% point outliers (+8σ) | .040 [.026-.061] | .066 [.047-.091] | .036 [.023-.056] | .036 [.023-.056] |
| 6% group outliers (+5σ) | .038 [.024-.059] | .060 [.042-.084] | .054 [.037-.077] | .054 [.037-.077] |

*Wilson 95% intervals use 500 seeds per cell. C/D are computed
from the cached SVD results; paper/review3_analysis.py §3.8c
documents the cell-level check.*

The recorded rates span .030-.066, with Wilson intervals shown
(500 seeds per cell). These selected S1 perturbations do not isolate
the causes of the S2 size distortions. Neither noise-related effects
nor interactions between noise, design, and inference are ruled out.

### 4.5 Do alarms predict extrapolation failure?

The S2 extrapolation pool contains 180,000 records:
13 H0 and 17 H1 tasks, positive noise, and all three sample sizes.
It excludes sentinel-invalid or missing-outcome tasks. Danger here
means the registered upper-end point error exceeds 10%, without
an R² screen; 36,206 records meet it (base rate .2011), including
37 records from true-form H0 tasks.

**Table 7 - Alarm precision and recall for >10% upper-end extrapolation error on the S2 pool, native thresholds; C/D use SVD. Declines produce no alarm and remain in the recall denominator.**

| Method | precision | recall | declines |
|---|---|---|---|
| A | .4226 | .6242 | 0 |
| B | .6379 | .4901 | 649 |
| C (SVD) | .6189 | .4639 | 649 |
| D (SVD) | .4879 | .8844 | 0 |

A and B have 10,296 and 5,442 exclusive true-positive alarms,
respectively, with 12,304 shared. These counts show different
alarm sets, not a claim of independence or chance-level overlap.
D detects more of this pool's dangerous outcomes, but the
precision/recall table does not compare methods at an equal
false-positive budget. Figure D.4 uses the same eligible pool.
An alarm and a non-rejection can both be informative under a
specified sampling distribution; neither supplies a bound on
error throughout an extrapolation domain.

### 4.6 The battery: from rejection to routing

**Table 8 - Probe: permutation vs wild null by family (n = 200, 200 seeds; entries rel = 1% / 3%).**

| Truth | perm null | wild null |
|---|---|---|
| power (size of power member) | .070 / .075 | .075 / .060 |
| rational: false escalation (V = ∅) | **.34 / .31** | .045 / .090 |
| rational: containment | .66 / .70 | **.96 / .91** |
| exponential: containment | .94 / .95 | .92 / .95 |
| logarithmic: containment | .87 / .93 | .89 / .95 |
| threshold (no battery family): escalation - *correct outcome* | 1.00 / 1.00 | 1.00 / 1.00 |

With the permutation null, the battery's own rational member rejects its
own true family a third of the time - nothing about rational laws is at
fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The
wild bootstrap repairs exactly that and nothing else.

**Set-level errors and multiplicity.** The battery's decision is the
viable *set*, so the errors that matter are set-level, and each member
runs at α = .05 without correction - a union decision whose two error
types move oppositely under multiplicity. Two pools, two behaviors.
On H0-side cells (truth is a battery member) the reported containment
is one minus the exclusion rate of the true family's own member
(.0496-.066 across the four truth families; the viable set is never
empty *on this pool*), and per-member α does not inflate that rate
because exclusion depends on one member only. Inflation appears on
the other side: some wrong member survives alongside the truth in
23.3% of power-truth cells (near-neighbour overlap:
exponential/logarithmic members retain power data) and 92% of
logarithmic-truth cells (rational retains log data) - reported as
near-neighbour overlap rather than suppressed in the locally recorded
tables. On the structural-H1 pool (72,000 cells whose truth is no
battery member by design) the viable set is empty - V = ∅, correct
escalation - at the rates of Table 10's first column: universally for
quadratic/threshold/omitted-variable (1.00), rarely for the
interaction/opposing/collinear structures (.001-.05) whose mismatch
directions the member families do not span. Routing there means
*which wrong family is modal* (Table 10, second column); the battery
is a candidate filtering stage whose surviving set the discovery
layer must still rank, not a final selection. Downstream selection
error after the battery is out of scope here and unmeasured.

**Table 9 - Battery on synthetic truths (500 seeds/cell, wild null).**

| Truth family | Containment (truth ∈ V) | False escalation (V = ∅) |
|---|---|---|
| power - size, 8 scenario groups incl. all noise layers | - | .0496 [.046-.053] pooled; per-scenario .043-.054 |
| exponential | .95 | .05 |
| logarithmic | .94 | .005 |
| rational | .93 | .066 pooled (≤ .082 worst rel) |

**Table 10 - Battery on out-of-family structures (escalation to V = ∅ is the *correct* outcome; column 1 is the V = ∅ rate).**

| Scenario | Correct escalation | Modal false V |
|---|---|---|
| quadratic / threshold / omitted variable | 1.00 / 1.00 / 1.00 | - |
| dilution | .92 | rational |
| cubic (odd term) | .24 | power (.68) |
| interaction | .05 | power (.95) |
| opposing curvature | **.00** | power+rational (.84) |
| collinear-quadratic | .03 | rational (.61) |

Table 10 is Table 4 re-expressed as routing errors: the power member
inherits A's documented blind spots, so those structures stay in V under
the power label instead of escalating. On the SRSD side the routing is
reported descriptively (median |V| 1-2 for 11 of 14 equations; i.30.5 →
{exponential} in 78%; ii.38.14 escalates on 67%), because surface-form
family labels dissolve under examination - the battery correctly rejects
Planck-type "exponentials" that no linearizable family owns.

### 4.7 The gate comparison: why not just compare against a nested rival?

**Table 11 - S5 size on H0 pools (200 seeds/cell; naive and calibrated thresholds for E/F).**

| H0 pool | A | E naive | E cal. | F naive | F cal. |
|---|---|---|---|---|---|
| synth n=50 (600 fits) | .080 | .212 | .050 | .132 | .050 |
| synth n=200 (600) | .040 | .183 | .050 | .135 | .050 |
| SRSD n=50 (2,600) | .051 | .167 | .050 | .086 | .050 |
| SRSD n=200 (2,600) | .038 | .152 | .050 | .123 | .050 |

**Table 12 - S5 structural-H1 detection (A at its native threshold; E/F at H0-pool-calibrated thresholds; mean over tasks).**

| | synth n=50 | synth n=200 | SRSD n=50 | SRSD n=200 |
|---|---|---|---|---|
| A (permutation curvature) | .774 | .782 | .348 | .493 |
| E (CV-quad) | .929 | .934 | .462 | .628 |
| F (AICc-quad) | **.930** | **.940** | **.525** | **.669** |

*OR(A, E) at each gate's own 5% threshold - a different, larger joint
error budget (.065-.102 across the four benches), reported separately:
.937/.941/.547/.675. It is not comparable to the rows above and is
excluded from this table.*

The error budgets are not exactly matched in Table 12: A retains
its native threshold, whereas E/F are tuned on the H0 reference
pool. In particular, A's synth n=50 size is .080, and its SRSD
n=200 size is below .05. Interval overlap with a target is not
proof of an identical budget. Table 12 reports measured detection
under these policies, not a uniformly matched-FPR ranking.

The quadratic rivals include per-feature squares and cross-terms,
covering directions absent from A's single mean-square column.
On the recorded S5 pool, F has the largest mean detection in
each displayed column. This is a scoped finding for one rival
family, not a universal default recommendation.

Naive E/F thresholds have elevated measured H0 rejection rates
(Table 11). Reference-pool calibration reduces them to the
reported targets; A requires no separate reference pool to
construct its resampling p-value, but that does not make A
size-valid on every pool. Its .080 small-sample result is one
counterexample.

The earlier SRSD leave-one-H0-equation-out analysis reports
held-out sizes .053/.050 for E at n=50/200 and .051/.052 for F,
with per-equation spread [.010, .095]. This provides a limited
transfer check across the 13 available H0 tasks, not calibration
on all new equations or a replacement for an independent
deployment test set.

For OR(A,E), using the individual thresholds increases the joint
H0 rejection budget to .065-.102. At synth n=50, A alone
already exceeds .05, so tightening E cannot achieve a 5% union
budget while preserving A's rule. The reported joint-budget
recalibration at n=200 gives .940/.664 detection on the
synthetic/SRSD benches, compared with F's .940/.669. Those
specific comparisons do not establish that F is the frontier
at every possible error budget. Combining gates and selecting
thresholds should be evaluated as a complete decision policy.

### 4.8 Real data: does the map predict live behavior?

**Table 13 - S6 World Bank relations: verdicts (median permutation p over 5 seeds for A; deterministic p for B-D; α = .05).**

| Relation (target ~ features) | Years all four gates agree | Divergent cells (which gate, which way) |
|---|---|---|
| infant mortality ~ GDP pc | 4/4 reject | - |
| CO₂ pc ~ GDP pc | 4/4 reject | - |
| internet penetration ~ GDP pc | 4/4 reject (p < 10⁻¹⁶) | - |
| GDP ~ population | 4/4 **quiet** (A .71-.91; B .19-.33) | - |
| life expectancy ~ GDP pc | 0/4 | A alone rejects 2019 (.035) & 2022 (.008); B alone rejects 2022 (.0002); C/D never |
| life expectancy ~ GDP pc + health exp pc | 2/4 quiet | B rejects 2022 (.0002) while A quiet (.108) |

Where signal is strong the gates agree - three relations rejected
unanimously every year, one candidate relation unrejected by all four gates -
and A's five p-values per cell are Monte Carlo replicates of one
estimator whose spread is small (median .009, max .094 at resolution
1/1000; no verdict flips beyond boundary cells). Divergence concentrates
where the synthetic map puts weak curvature: the life-expectancy cells
flip year to year (signal at the detection boundary), and the
two-feature cell shows B rejecting where A stays quiet. We label this
pattern *consistent with* the dilution mechanism (§4.2, Table 4) rather
than a reproduction of it: on observational data the p-value change
after adding a second feature can also reflect feature correlation,
omitted-variable changes, sample composition, measurement error,
cross-country heterogeneity, or spatial dependence - no ground truth
exists to separate these, so the claim is signature-consistency, not
mechanism confirmation.

**Table 14 - S6 dilution probe: the same target with and without a second correlated curved feature (3 years × 5 seeds; detection = median gate verdict).**

| Relation | A | B | D |
|---|---|---|---|
| cereal yield ~ fertilizer (1 feature) | rejects 3/3 years (.011-.031) | rejects 3/3 (~.003) | rejects 3/3 |
| cereal yield ~ fertilizer + energy pc (2 features) | **quiet 3/3 years** (.067-.149) | rejects 2/3 (.025-.092) | rejects 3/3 (≤.018) |

This matches the S1 DILUTION scenario's signature - adding one
correlated curved feature silences A on a relation it was rejecting,
while the RESET variants continue to alarm - a pattern *consistent
with* the dilution mechanism rather than a demonstrated reproduction of
it, since feature correlation, sample composition, and omitted-variable
changes move the same p-values in the same direction. No ground truth
exists on S6 and none is claimed; the contract of §3.6 is signature
consistency in both directions: unanimous cells where the map predicts
signal, and map-consistent divergence where it predicts weakness.

A companion engineering bench (Appendix E) pushes the same probe across
nine physical-engineering sources (147 run cells): unanimous detection of
a documented rating-curve misspecification, mechanism-echo divergences
matching the synthetic map (small-sample power, interaction carried by a
second feature, all three HC3 failure modes live), and mutual misses on
locally approximate laws - all recorded, none adjudicated.

### 4.9 Runtime

The historical S1 runs recorded per-dataset means of 12.0-25.3 ms
for A and 1.0-1.8 ms for B/C/D; the family battery at n=200
recorded 24 ms, and S5 E/F 1.5-3.5 ms. These figures describe
the recorded environments and implementations. No inference-path
runtime equivalence or production-throughput guarantee is claimed.
Figure D.5 is retained as a historical timing panel; the SVD
revision requires a separately labeled timing comparison.

### 4.10 What the full pipeline actually releases

The joint-release analysis links fit quality, native A/B decisions,
SVD C/D decisions, and extrapolation errors one-to-one on
(task, n, rel, seed). It requires positive noise, a valid
statistic-consistency sentinel, and finite R² and far-field outcomes.
The 17 eligible H1 tasks contribute 34,000 records per n.
The i.30.5 statistic-sentinel records and the i.26.2 records without
finite far-field outcomes are excluded under these rules.
audit_results.py records the join and all filters.

**Decision events.** Release requires R² ≥ 0.99, four valid
gate evaluations, and no rejection at the native alpha=.05 rules.
A dangerous wrong model has relative error >10% at the
registered upper-end probe (Section 3.3). A declined diagnostic
withholds release. The error threshold defines the consequence
evaluated here; the probe measures one input rather than an
operating-domain maximum.

**Table 15 - Joint release on the 17-task S2 H1 pool. Each n has 34,000 records. Count symbols S, Sd, Ld, and Da distinguish the event populations from gate labels A-D. Eligibility excludes i.30.5 (statistic-sentinel failure) and i.26.2 (no finite far-field outcome); Table 16 uses the 18-task H1 pool without the finite-far requirement.**

| n | Screened S | Dangerous in S: Sd | Released dangerous Ld | All dangerous Da | Ld / S | Ld / Sd | Ld / Da |
|---|---|---|---|---|---|---|---|
| 50 | 22,060 | 2,208 | 350 | 11,922 | 1.59% | 15.85% | 2.94% |
| 200 | 21,932 | 2,150 | 106 | 12,202 | 0.48% | 4.93% | 0.87% |
| 1000 | 22,000 | 2,045 | 25 | 12,045 | 0.11% | 1.22% | 0.21% |

S counts screened wrong models; Sd counts dangerous wrong
models within S; Ld counts those passing both the screen and
all four valid gates; Da counts every dangerous model in the
eligible H1 pool. The three ratios answer distinct operational
questions: frequency among screened wrong forms, failure to
intercept screened-dangerous forms, and survival among all
dangerous fits. Pooled across n, S=65,992, Sd=6,403,
Ld=481, and Da=36,169. The respective rates are 0.73%,
7.5%, and 1.33%.

**The fit screen changes the release population.** At n=50,
966 dangerous records have four valid non-rejecting gates;
616 fail the fit screen, leaving 350 complete-workflow releases.
At n=200 and 1000, the gate-only and complete-workflow dangerous
release numerators coincide at 106 and 25. Across n, the panel
declines on 2, 56, and 591 eligible H1 records. This accounting
separates rejection, diagnostic withholding, fit screening, and
complete release.

**Detection benefit and residual failure.** Among screened-dangerous
wrong models, 481/6,403 (7.51%) are released. Among screened wrong
models below the danger threshold, 27,025/59,589 (45.35%) are
released. The panel therefore excludes most screened-dangerous
records from release while retaining a substantial residual channel.
These pooled rates summarize the sampled task composition, and
that composition is concentrated. The 481 releases come from two
of the seventeen tasks: bonus.15 (242) and bonus.5 (239). Six
of the remaining fifteen release no dangerous record at any n while
contributing 30,351 dangerous records between them; among these,
all 557 dangerous records of ii.38.14 are intercepted.
Structural mismatch and failure at the chosen probe are different
labels, and their relation varies by equation.

**Sensitivity to the danger threshold.** At n=50, thresholds
of 2%, 10%, and 20% label 34.91%, 10.01%, and 2.76%
of screened H1 models dangerous. The release fractions within
those subsets are 41.30%, 15.85%, and 6.40%. At n=1000
and the 20% threshold, Sd=0 and Ld/Sd is undefined.
The complete 2/5/10/20/50% sweep is in evidence/results.json.

A pooled task-level sweep makes the composition visible. At the
10% threshold, all 481 releases come from bonus.15 (242) and
bonus.5 (239). At 12.5%, the count falls to 176: bonus.15
contributes 5 and bonus.5 contributes 171. At 15%, all 111
releases come from bonus.5; at 20%, all 39 do. The released set
therefore consists of the two tasks throughout, but which task
dominates changes sharply between 10% and 12.5%. These counts
are reproduced by threshold_decomposition.py and recorded in
evidence/threshold_decomposition.json.

**Sample-size pattern.** Ld/Da decreases from 2.94% to 0.87%
to 0.21%, and Ld/Sd decreases from 15.85% to 4.93% to
1.22%. A fixed-label control assigns each matched (task, rel,
seed) triple its n=1000 danger label. Across 34,000 complete
triples, Ld/Sd is then 16.77%, 1.87%, and 1.22%, while
Ld/Da is 2.93%, 0.31%, and 0.21%. The decreasing pattern
survives this label control. Its interpretation still combines
fit-screen behavior, diagnostic sensitivity, and changing probe
coordinates; training samples are not constructed as nested sets.

**Task composition of the released set.** The released records
change identity with n. At n=50 the 350 releases are 239 from
bonus.5 and 111 from bonus.15. At n=200 and n=1000 every release
(106 and 25) comes from bonus.15 alone. Excluding that task,
Ld/Sd is 11.58%, 0%, and 0%. The fixed-label control holds the
danger label but not the task composition: bonus.15 has a median
relative far-field error of 9.2%, with the 75th and 90th
percentiles at 9.5% and 9.8% against the 10% threshold. Its
danger label therefore moves with the threshold while its fitted
form and task identity do not change. The decreasing pattern
combines a panel effect with the position of one task's error
distribution relative to a fixed cutoff.

**Task-level uncertainty.** The record-level ratios treat simulation
replicates as independent observations. A task bootstrap instead
resamples the 17 eligible equations with replacement (B=10,000),
carrying each task's S/Sd/Ld/Da and recomputing the pooled ratios.
Pooled over n, Ld/Sd is 7.5% with a 95% interval of [0%, 49.9%];
12% of resamples contain no released dangerous record, and 15 of the
17 tasks contribute Ld=0. Leave-one-task-out values span
[3.9%, 20.1%]: omitting bonus.5 raises Ld/Sd to 20.1%, while omitting
bonus.15 lowers it to 3.9%, because the remaining task's own release
fraction then dominates. At n=200 and 1000 the interval reaches
[0%, 75.7%] and [0%, 56.8%], since only bonus.15 contributes releases.
This interval measures between-equation heterogeneity, not binomial
precision, and the record-level 7.5% describes the sampled task pool
rather than a task-uniform rate. Counts and intervals are reproduced by
audit_task_intervals.py and recorded in evidence/task_intervals.json.

**Information in a clean diagnostic record.** The corresponding
180,000-record H0/H1 pool has three states: true H0
(78,000 records; 65,934 releases), dangerous wrong H1
(36,169; 481 releases), and non-dangerous wrong H1
(65,831; 27,025 releases). The empirical release probability
is .013299 for dangerous H1 and .845308 for true H0.
Their likelihood ratio is .01573 in the dangerous-versus-true
direction, or 63.56 in the reverse direction.

Across all three states, 481 of the 93,440 released models
are dangerous wrong forms (0.515%). This proportion depends
on the pool's class composition. In a different, explicitly
hypothetical population containing only true H0 and dangerous
H1 with equal priors, Bayes' rule gives a 1.55% dangerous
posterior after release. These are separate populations and
their conditional risks must remain separate. Applying the
decision rule to a new population requires its own validation.

The results support a complete-workflow evaluation: diagnostics
carry information and remove dangerous candidates, while a
clean record leaves measurable residual uncertainty.

### 4.11 Joint silence without the fit screen

Panel non-rejection is measured separately from complete release
using the S7 arms defined in Section 3.7. audit_s7.py reproduces
these counts from the existing records.

**Table 16 - A-D joint silence on S2 records. H1 has 18 eligible tasks (36,000 records per n); H0 has 13 tasks (26,000 per n). Calibrated means a 5% H0-reference target per gate.**

| n | H1 native alpha=.05 | H1 per-gate calibrated | H0 per-gate calibrated |
|---|---|---|---|
| 50 | 13,616 / 36,000 (37.8%) | 13,888 / 36,000 (38.6%) | 22,821 / 26,000 (87.8%) |
| 200 | 10,717 / 36,000 (29.8%) | 10,877 / 36,000 (30.2%) | 23,102 / 26,000 (88.9%) |
| 1000 | 8,355 / 36,000 (23.2%) | 8,478 / 36,000 (23.6%) | 23,108 / 26,000 (88.9%) |

A silent record has finite p-values for all four gates and no
rejection. The H0 rates are measured on the reference pool
used to set the thresholds. The resulting joint H0 rejection
exceeds 5%, illustrating the distinction between individual
calibration targets and the error budget of a panel.

**Paired A versus A/E/F.** On the S5 SRSD pool at n=50,
adding E/F reduces H1 silence from 2,479/3,800 (65.2%) to
1,602/3,800 (42.2%). H0 silence decreases from 2,468/2,600
(94.9%) to 2,352/2,600 (90.5%). At n=200, H1 silence
decreases from 1,926 to 1,160 out of 3,800, and H0 silence
from 2,500 to 2,377 out of 2,600. The paired comparison
demonstrates increased detection together with increased false
rejection. It motivates the equal-joint-budget validation in
Section 6.

This three-gate arm has 19 H1 tasks at one noise level;
the A-D arm has 18 H1 tasks at four noise levels. The
S5 cache has no recorded failures for i.30.5 but lacks
S2's separate replay-sentinel fields. Its within-arm counts
are retained as exploratory evidence. Comparison of complete
A-D and A-F panels requires common records.

**Class-conditional interpretation.** Within the S5 A/E/F
arm at n=50, silence favors true over wrong forms by the
likelihood ratio (2352/2600)/(1602/3800)=2.146.
Under an explicitly hypothetical 50/50 class prior, the
posterior probability of a wrong form is
(1602/3800)/[(1602/3800)+(2352/2600)]=31.8%.
The likelihood ratio measures evidence under this policy;
the posterior also uses the assumed prior.

Table 15 additionally imposes R² ≥ 0.99 and requires a
finite far-field outcome. At n=50, 22,060/34,000 eligible
H1 records pass that screen (64.9%). Gate-only silence
and screened dangerous release therefore answer different
questions with their own eligible populations.

## 5 Discussion

### 5.1 Match diagnostics to the mismatch direction

The controlled alternatives explain why a single method ranking is
insufficient. A's mean-square column can cancel opposing curvature and
has weak signal for odd-symmetric departures; interactions distinguish
fitted-value augmentation from regressor-power augmentation
(Table 4). Resampling calibrates a chosen statistic but does not
expand the directions it examines. The task-level differences
therefore matter alongside the average detection rates.

D is the strongest mean performer among A-D on the common S2
pool under the stated reference target. Its paired advantage over A
is supported at n=200 and 1000. S5 separately shows that a richer
quadratic rival has higher recorded mean detection than A in each
displayed pool. Together these findings motivate testing both
augmentation geometry and model-comparison alternatives. A paired
D-versus-F ranking requires the common-policy experiment in Section 6.

### 5.2 Separate statistical assumptions from numerical stability

Native size is an empirical requirement: Tables 1-2 show departures
from .05 even when the symbolic mean is a true monomial.
The family-viability experiment supplies a second example. In rational
coordinates, changing the resampling null substantially improves
observed containment (Table 8). Transforming coordinates therefore
requires reevaluating the noise and resampling assumptions.

Numerical implementation supplies a distinct failure mode.
A stratified cross-check of exog-HC3 on 176 fits spanning all
13 H0 equations compared the original pinv(X′X) path,
an SVD path, and statsmodels on identical data. The SVD and
statsmodels calculations agreed closely (maximum absolute
p-value difference .0023), whereas the original path produced
large rejection rates on near-collinear designs. Forming X′X
amplifies conditioning problems; the stable path changes the
empirical method comparison. The audit supports using the SVD
results reported here. The complete measurements and earlier
interpretation are preserved in Appendix G.

### 5.3 Measure useful rejection and residual release together

The diagnostics provide useful information. On the S2 evaluable
pool, each gate's alarm precision for the specified far-field
error exceeds the .2011 base rate (Table 7). After fit screening,
481/6,403 dangerous wrong models survive the four-gate policy.
The same policy also releases many wrong forms with a small error
at the chosen probe. These observations distinguish a structural
label from a consequence measured at a particular input.

A useful workflow assessment reports both sides: dangerous models
removed, dangerous models released, correct models retained, and
declined evaluations. Table 15 connects these events on the same
records. Table 16 supplies the complementary gate-only view.
More gates can reduce silence while increasing false rejection;
a common joint budget is needed to decide whether the added
checks improve the trade-off.

Conditional rates retain their populations. The 7.5% screened-dangerous
release rate measures failure to intercept dangerous candidates
that already meet the fit screen, but it is a pooled figure driven by
two tasks; excluding the threshold-borderline task bonus.15, the
release fraction is 0% at n=200 and 1000. A task-level bootstrap over
the 17 equations gives a 95% interval of [0%, 49.9%] for the pooled
7.5%, so the pooled figure is descriptive rather than precision-
estimated. The 0.73% rate measures their frequency among all screened
wrong forms. Deployment risk
additionally depends on task prevalence and on the population of
correct candidates, as the three-state accounting in Section 4.10
demonstrates.

### 5.4 Use family viability and observations as diagnostic follow-up

The battery returns a viable family set after coordinate-specific
testing. On the recorded synthetic family truths, containment is
approximately .93-.95; on out-of-family alternatives, Table 10
identifies structures that escalate and structures that retain a
wrong family. This supplies a filtering stage for further model
comparison, with overlap explicitly recorded.

The World Bank and engineering suites extend the software and
diagnostic record to observed data. Their agreement and divergence
patterns can guide follow-up measurements or alternative models.
Without known form labels or controlled interventions, those
patterns supply compatibility evidence. The strongest empirical
performance conclusions remain those supported by the
ground-truth suites.

## 6 Scope and validation priorities

### 6.1 Scope of the current evidence

The main comparison concerns a prespecified power-law family on
positive data. S4 extends the diagnostic construction to four
linearizable families at n=200; S5 examines one quadratic rival.
The results describe these candidate families, the recorded sampling
configurations, and their threshold policies. Residual permutation
provides an empirical approximation under the evaluated conditions.
More strongly dependent errors, additive-error models, or different
rivals require their own calibration evidence.

S2 contains 18 eligible H1 tasks for detection and joint silence,
and 17 with finite extrapolation outcomes for joint release.
The excluded i.30.5 statistic-sentinel records and the missing
i.26.2 far-field outcome remain in the audit trail. Repeated
simulations estimate within-task behavior; task bootstrap intervals
describe variation within the available equation pool. S3's 30
seeds per cell and S6's Monte Carlo repeats have correspondingly
limited precision and scope.

The far-field outcome is a single probe at three times the observed
training maximum of the first feature, with other log-features held
at training medians. Consequently the evaluation coordinate can
change across sample sizes and seeds. Threshold sensitivity and
fixed-danger-label controls strengthen the descriptive result, but
they do not turn this probe into an operating-domain error bound.

### 6.2 Two priority validation experiments

**Independent evaluation of a common joint policy.** The next
comparison should compute the candidate gates on the same fitted
records, choose thresholds on a separate H0 calibration split,
freeze the policy, and evaluate on independent H0/H1 records.
Strong single gates such as D and F should be compared with
prespecified panels at the same target joint false-rejection
budget. The final test should report actual false rejection,
dangerous release, retained-candidate fraction, and declines,
with paired task-level uncertainty. Holding out equations is
needed when the intended claim concerns new equations.
The existing S5 leave-one-equation checks provide limited
transfer evidence to build on.

**Evaluation on a fixed extrapolation domain.** For each
ground-truth equation, a physically valid external evaluation grid
should be fixed before drawing the training samples. The same
grid and prespecified error summary should be used across
sample sizes, alongside the current endpoint probe. This would
test the stability of the release pattern over an explicit domain.
A paired training-sample design would further support analysis
of the sample-size effect. Neither validation experiment is
reported as completed in this manuscript.

### 6.3 Extension to discovered candidates and engineering decisions

A discovery-workflow claim requires evaluation of the expressions
actually selected by a search algorithm. Candidate generation and
selection should be frozen before independent diagnostic and
prediction evaluation, and the diagnostic null should correspond to
the candidate being assessed. The current main benchmark measures
family diagnostics; algorithm-selected expression risk is a
distinct next experiment.

The observational suites demonstrate data compatibility and
diagnostic behavior. Engineering recommendations additionally
require reference outcomes over the intended operating domain,
relevant noise models, and a loss-based tolerance. These steps
would connect the benchmark's measured failure channels to
a specific application.

## 7 Reproducibility statement

The local benchmark source tree is D:/form_check. It contains
generators, recorded configurations, raw and cached results,
failure labels, and earlier analysis scripts. Appendix A reports
the original fitted-run totals. S7 reuses S2/S5 records and adds
no fitted datasets to those totals.

The manuscript bundle contains manuscript.md, build_pdf.py,
audit_size.py, audit_results.py, and audit_s7.py. The audit scripts
record source-file hashes, eligibility rules, thresholds, and
counts in evidence/size.json, evidence/results.json, and
evidence/s7.json. They reproduce the reported aggregation and
Figures D.1-D.4 from the existing caches. Figure D.5 retains
the historical timing measurements, with its inference-path
scope stated. The build script renders the manuscript and figures.

Within the original source tree, runs/cd_svd/ contains the C/D
SVD outputs; runs/srsd_full/ contains the S2 records and
quality_recompute.parquet. S5 gate definitions are in
gate_compare.py; joint_silence.py and joint_silence_toolbox.py
identify the S7 arms. wb_run.py and wb_blindspot.py document
the observational macro-data workflow; real_data/ contains
engineering data manifests and source-specific license records.

The current validation includes record-level joins, input hashes,
independent aggregation checks, selected implementation comparisons,
and PDF consistency checks. Reproduction of the original fitting
experiments is separate from re-aggregation of their stored outputs.

The artifact is publicly available at
https://github.com/leokynlang/form_check. The data snapshot
accompanying this manuscript is release tag data-2026-10-18.
data_release/manifest.json lists every member path, size, and
SHA-256, and SHA256SUMS.txt covers the archive. The repository README
gives the clone, install, data-restoration, and aggregation commands,
including audit_size.py, audit_results.py, audit_s7.py,
reproduce_r2.py, verify_composition.py, threshold_decomposition.py,
deep_analysis.py, and audit_task_intervals.py.

Third-party data are not uniformly redistributed:
data_release/THIRD_PARTY_DATA.md and dataset_sources.json record
per-source decisions, and copyrighted sources such as WebBook SRD 69
are source-download-only. Appendix G preserves the numerical and
interpretive correction history.

## 8 Conclusion

This benchmark quantifies what selected specification gates detect,
what false alarms they produce, and which dangerous wrong models
survive the complete decision workflow. Controlled alternatives
identify distinct blind spots. On the common 18-task physics pool,
SVD exog-HC3 has the highest mean detection fraction among A-D
at the stated H0-reference target, with positive paired D-A
intervals at n=200 and 1000. Native-size and numerical checks
show why calibration and stable inference must accompany
that comparison.

The joint-release analysis measures the consequence of combining
fit screening with valid diagnostic non-rejection. Of 6,403 screened
dangerous wrong models, 481 pass all four gates: a 7.5% residual
release rate under the specified endpoint criterion. The
corresponding frequency among all screened wrong models is 0.73%.
The released records are concentrated in two tasks; excluding the
threshold-borderline task bonus.15, no dangerous model is released
at n=200 or 1000, so the pooled 7.5% should not be read as a
task-uniform failure rate. The gates remove many dangerous
candidates, yet the remaining failure channel is observable and
depends on the tested structure and decision policy.

The practical contribution is a framework for assessing diagnostic
trade-offs on the same models, with failure registration and explicit
conditional populations. Independent evaluation of common-budget
panels and a fixed extrapolation domain are the next steps toward
workflow recommendations. The benchmark provides concrete
failure cases, comparative evidence, and a reproducible basis for
those tests.

## References

[1] Ramsey, J. B. (1969). Tests for specification errors in classical linear least-squares regression analysis. *Journal of the Royal Statistical Society, Series B*, 31(2), 539-551.

[2] White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4), 817-838.

[3] MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305-325.

[4] Long, J. S., & Ervin, L. H. (2000). Using heteroskedasticity consistent standard errors in the linear regression model. *The American Statistician*, 54(3), 217-224.

[5] La Cava, W., Orzechowski, P., Burlacu, B., de França, F. O., Virgolin, M., Jin, Y., Kommenda, M., & Moore, J. H. (2021). Contemporary symbolic regression methods and their relative performance. *NeurIPS Datasets and Benchmarks Track*.

[6] Matsubara, Y., Chiba, N., Igarashi, R., & Ushiku, Y. (2024). Rethinking Symbolic Regression Datasets and Benchmarks for Scientific Discovery. *Journal of Data-centric Machine Learning Research*. https://data.mlr.press/assets/pdf/v01-3.pdf. Preprint first submitted as arXiv:2206.10540 in 2022.

[7] Udrescu, S.-M., & Tegmark, M. (2020). AI Feynman: A physics-inspired method for symbolic regression. *Science Advances*, 6(16), eaay2631.

[8] Bakarji, J., Champion, K., Kutz, J. N., & Brunton, S. L. (2022). Discovering governing equations from partial measurements with deep delay autoencoders. *Nature Computational Science*.

[9] Cranmer, M. (2023). Interpretable machine learning for science with PySR and SymbolicRegression.jl. *SciPy Conference*.

[10] Liu, R. Y. (1988). Bootstrap procedures under some non-i.i.d. models. *Annals of Statistics*, 16(4), 1696-1708.

[11] Davidson, R., & Flachaire, E. (2008). The wild bootstrap, tamed at last. *Journal of Econometrics*, 146(1), 162-169. https://doi.org/10.1016/j.jeconom.2008.08.003.

[12] Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach* (2nd ed.). Springer.

[13] Arlot, S., & Celisse, A. (2010). A survey of cross-validation procedures for model selection. *Statistics Surveys*, 4, 40-79.

[14] World Bank (2026). World Bank Open Data API, indicators NY.GDP.MKTP.PP.CD, NY.GDP.PCAP.PP.CD, SP.POP.TOTL, SP.DYN.LE00.IN, SP.DYN.IMRT.IN, EN.GHG.CO2.MT.CE.AR5, EN.GHG.ALL.MT.CE.AR5, SH.XPD.CHEX.PC.CD, IT.NET.USER.ZS, AG.YLD.CREL.KG, AG.CON.FERT.ZS, EG.USE.PCAP.KG.OE. https://api.worldbank.org (snapshot cached in the local benchmark artifact; see Section 7).

[15] UCI Machine Learning Repository. archive.ics.uci.edu (datasets: auto_mpg, automobile, concrete, slump, energy_efficiency, yacht, ccpp, air_quality, gas_turbine, tetouan, steel). Accessed 2026-09-06; raw snapshots and SHA-256 hashes recorded in the local artifact; redistribution subject to each dataset's stated license.
[16] NIST Statistical Reference Datasets (StRD). itl.nist.gov/div898/strd. Accessed 2026-09-06; public domain.
[17] NIST Chemistry WebBook, SRD 69. webbook.nist.gov. Accessed 2026-09-06. Copyrighted Standard Reference Database 69; raw pages and extracted tables are not redistributed in the public artifact.
[18] U.S. Geological Survey. National Water Information System (NWIS instantaneous values, 12 gauges, WY2023) and FDSN event service (global M ≥ 5.5, 2000-2024). earthquake.usgs.gov / waterdata.usgs.gov. Accessed 2026-09-06; public domain (U.S. Government work).
[19] Open-Meteo historical weather API (ERA5 reanalysis), 10 wind sites, 2023. open-meteo.com, CC BY 4.0 (ERA5: Copernicus licence). Accessed 2026-09-06; reanalysis, not in-situ measurement.
[20] OpenML. openml.org (datasets: ailerons, cpu_act, elusage). Accessed 2026-09-06; per-dataset licenses vary and are recorded in the local manifest.
[21] U.S. Department of Energy, fueleconomy.gov vehicle fuel-economy data, model years 1984-2025. Accessed 2026-09-06; public domain (U.S. Government work).
[22] RDatasets (Vincent Arel-Bundock), datasets: pressure, cars, rock, trees, attenu, housing, electricity, stormer, whiteside, cement, petrol, women, airquality, puromycin. vincentarelbundock.github.io/Rdatasets. Accessed 2026-09-06; licenses inherit from each source package (mostly GPL/BSD) and are recorded in the manifest.

## Appendix A - Suite composition

| Suite | Cells | Fits per cell | Total fits | Extra evaluations |
|---|---|---|---|---|
| S1 synthetic | 285 | 500 | 142,500 | - |
| S2 SRSD-Feynman | 480 | 500 | 240,000 | 240,000 test sets; 232,500 extrapolation probes |
| S3 family-stratified | 392 | 30 | 11,760 | - |
| S4 battery | 188 | 500 | 94,000 | - |
| S5 gate comparison (v3) | 102 | 200 | 20,400 | - |
| S6 real data + probe (v3) | 33 | 5 | 165 | - |
| **Total** | **1,480** | - | **508,825** | **472,500** |

Outside the seeded totals: the Appendix E engineering bench contributes
150 single-realization smoke cells (147 run), one fit per gate per cell.

## Appendix B - Failure registration taxonomy

| Reason | Meaning | Exercised by |
|---|---|---|
| `degenerate_zero_residual` | Zero-noise cell on exact ground truth; test undefined | S1/S2 σ = 0 cells |
| `rank_deficient` | Design matrix numerically rank-deficient | collinear designs |
| `insufficient_df` | n − k < 1 after augmentation | n = 50 multi-feature cells |
| `hc3_leverage_degenerate` / `hc3_singular_covariance` | HC3 undefined | degenerate built columns |
| `statistic_mismatch` | Wrapper statistic ≠ native statistic (drift sentinel) | continuous audit |
| `checked_false_*` | Native test declines to run (non-positive data, n < 12, non-finite) | not-applicable tasks |
| `runner_error` / generation failure | S5 SRSD-side generation exceptions | S5 |
| indicator archived / no data | S6 indicator not published for a year | S6 (`EN.ATM.CO2E.KT` → AR5 replacement, logged) |

All registered failures remain in the audit trail; denominators are
documented per table. The battery's 94,000 and S5's 20,400 datasets
contain zero data-layer failures.

## Appendix C - Statistical uncertainty and eligible populations

Tables 1 and 2 contain the fit-level Wilson 95% intervals
recomputed by audit_size.py, with 6,000 and 26,000 valid
positive-noise H0 evaluations per n, respectively. The S1
zero-noise declines are excluded from size estimation.
These intervals describe simulation precision on the fixed
scenario/equation pools, not cross-equation generalization.

At S2 n=200, the equation-level analysis in audit_results.py
resamples the 13 H0 equations 10,000 times, with 2,000
records per equation. For C, the mean is .053115 and the
95% task-bootstrap interval is [.050308, .056385].
For D, the mean is .056385 and the interval is
[.054615, .058385]. Leave-one-equation-out ranges are
[.052000, .053542] for C and [.055625, .056667] for D.
Section 4.1 uses the same results with rounded presentation.
A lower bound displayed as .050 is not proof of exact
coverage of .05.

The paired D-A intervals in Section 4.2 instead resample
18 H1 equations, preserving each equation's paired difference.
They must not be substituted for fit-level binomial intervals,
or for results on a 19-task pool. Pointwise Figure D.3 intervals
and paired differences address different quantities.

The joint-release ratios have a separate task-level analysis.
audit_task_intervals.py resamples the 17 eligible H1 equations
(B=10,000) and recomputes pooled Ld/S, Ld/Sd, and Ld/Da from each
task's counts. With releases confined to two equations, the pooled
Ld/Sd of 7.5% has a 95% interval of [0%, 49.9%], and 12% of
resamples contain no release. This interval describes between-
equation heterogeneity; it is not a binomial interval on 34,000
records and does not include threshold-estimation uncertainty.

Seeded computational repeats within one equation do not create
additional independent laws. Macro and record-weighted averages
coincide only where the task design is balanced. The old
pinv-path uncertainty results are historical numerical-path
diagnostics, not current method-performance intervals.

## Appendix E - The real-engineering data bench (companion to S6)

S6 probes mechanism transfer on macroeconomic aggregates. The deployment
domain these gates are actually built for is physical-engineering
measurement - fluid, civil, hydraulic, energy, meteorological, seismic,
chemical, and aerospace data, where power laws are the daily candidate
models. This appendix documents a companion bench of that kind, run under
the same contract as S6: real data has no ground truth, so the bench
registers compatibility and divergence, not error rates; its cells are
single realizations and stand outside Appendix A's seeded totals.

**Composition and discipline.** Nine open sources with source-specific
licenses; 91 datasets; 150 candidate log-log relations registered in the local
manifest with per-relation physical notes; 761,515 raw rows (Table E1).
Raw bytes are cached (~115 MB) with SHA-256 hashes; cleaning is replayable
and imputation-free - rows with missing or non-positive values are dropped
and logged, and the few unit conversions are explicit in the manifest
(e.g., Whiteside °C → K). Simulator datasets that violate the real-data
discipline (a marine-engine and two wave-energy numerical models) are
excluded and listed, as are WebBook fluid pages below the 10-row
rejection-detection floor; one gauge whose instantaneous series was never
published is logged as an empty skip. ERA5 cells are marked as reanalysis,
not direct measurement.

**Table E1 - Composition by source.**

| Source | Datasets | Relations | Rows | Examples |
|---|---|---|---|---|
| USGS NWIS rating curves | 12 | 12 | 442,269 | discharge ~ gage height, 12 gauges |
| UCI | 11 | 38 | 145,920 | CCPP, Tetouan power, steel, yacht, concrete |
| EPA fueleconomy | 1 | 4 | 48,624 | CO₂ ~ MPG; combined MPG ~ displacement |
| Open-Meteo / ERA5 | 10 | 20 | 87,600 | 100 m ~ 10 m wind speed (shear); gust factor |
| OpenML | 3 | 8 | 21,997 | aileron position; CPU utilization |
| USGS FDSN | 1 | 3 | 12,334 | felt reports, stations, intensity ~ magnitude |
| RDatasets | 14 | 26 | 1,349 | cars, trees, puromycin, stormer, rock |
| NIST StRD | 13 | 13 | 968 | certified nonlinear calibration curves |
| NIST WebBook | 26 | 26 | 454 | phase-transition enthalpies ~ temperature |
| **Total** | **91** | **150** | **761,515** | - |

**Protocol.** Each cell applies the §3.2 gates unchanged at α = 0.05 - A
with its B = 999 permutation p, B/C/D as registered - on the
positivity-filtered relation, subsampled to n ≤ 800 (seed 42) to bound
permutation cost. The run is a smoke pass: one realization per cell,
compatibility plus first verdicts, not a rate estimate. Of 150 cells, 147
run; three are registered `insufficient_n` by design (Bennett5 and
Roszman1 have no positive rows after the log-space positivity filter;
Danwood has n = 6) - the bench's honest-failure probes working as
intended.

**Table E2 - Smoke verdicts on the 147 run cells.**

| Gate | Rejections | Registered failures |
|---|---|---|
| A - permutation curvature | 97/147 (.660) | 10 declines, all `checked_false_small_n` (n ≤ 11) |
| B - RESET fitted F | 111/147 (.755) | 0 |
| C - RESET fitted HC3 | 108/147 (.735) | 0 |
| D - RESET exog HC3 | 94/147 (.639) | 3: `insufficient_df` (n = 13, four features), `hc3_leverage_degenerate` (n = 546), `rank_deficient` |

A and B agree on 123/147 cells (92 both reject, 31 both quiet; 5 A-only,
19 B-only). The failure columns are themselves informative: D's three
declines are exactly the three HC3 failure modes of Table 2 and Appendix
B, now on observational data; A's ten are uniform declines-to-run on
n ≤ 11 cells; B and C never fail. Nothing fails silently - every declined
evaluation is registered with a typed reason.

**Designed regimes.** The bench's design document assigns expected regimes
to 66 of the 150 cells before the smoke pass (Table E3). The assignment is
a reading aid, not a scored test: the H0-side laws are approximate by
design, and at n up to 800 an approximate law is a detectable null. The
12 rating-curve cells - whose registered form Q ~ h omits the documented
stage offset h₀ of the true rating law Q = C(h − h₀)ⁿ, a mismatch noted
in the manifest - are rejected unanimously by A and B alike (12/12):
detection of a documented misspecification, not a false alarm. The five
H0-side cells A retains are those with no designed departure mechanism -
the braking-distance v² law, the tree allometry, the road-load relation -
plus the two shear sites with the smallest fitted curvature.

**Table E3 - Designed-regime cells, rejections at α = 0.05.**

| Designed regime | Cells | A | B | C | D |
|---|---|---|---|---|---|
| Approximate power laws (H0-side) | 30 | 25 | 27 | 23 | 23 |
| Known non-power-law physics (H1-side) | 36 | 20 | 24 | 24 | 14 |

(A's H1-side count includes its nine WebBook declines at n ≤ 11 as
non-rejections.) On the H1 side no gate sweeps the board: the misses
concentrate in WebBook enthalpy cells whose measured temperature window
keeps the local power-law approximation good - e.g., toluene Δ_vapH ~ T:
all four gates quiet at n = 26 - a mutual-miss cell recorded as measured.

**Table E4 - Named cells: mechanism echoes and divergences; p-values at α = 0.05 in A / B / C / D order, - = declined with registered reason.**

| Cell | A / B / C / D | Reading |
|---|---|---|
| Wind shear v₁₀₀ ~ v₁₀, 10 ERA5 sites (n = 800 ea.) | 8/10 / 10/10 / 7/10 / 7/10 | B flags every site (p ≤ 2×10⁻⁴); A retains exactly the two smallest-curvature sites (.23, .41) - a graded response, not a threshold artifact |
| Gust factor ~ v₁₀, same sites (n = 800 ea.) | 10/10 ×4 | unanimous rejection at every site; no truth to adjudicate |
| puromycin rate ~ conc (Michaelis-Menten), n = 23 | **.011** / .069 / .169 / .169 | A detects the saturation curvature, B misses at α = .05 - Table 3's small-sample power gap, live |
| stormer Time ~ Viscosity + Wt (calibration law), n = 23 | .001 / 7×10⁻⁴ / 9×10⁻³ / - | single-feature cells quiet (A .94, B .99); the second feature carries the interaction the law contains |
| EPA CO₂ ~ MPG (stoichiometric inverse, exponent −1), n = 800 | .030 / 1.5×10⁻⁴ / .727 / .727 | A and B reject on a near-exact physical law bent by fleet drive-mix curvature |
| EPA combined MPG ~ displacement (road load), n = 800 | .386 / .596 / .599 / .599 | designed power law passes |
| CCPP PE ~ V, n = 800 | .642 / 8×10⁻⁴ / ≈0 / .023 | proxy curvature: V correlates .81 with the strongly curved AT, so B/C's fitted-power columns inherit AT's curvature while A's V-only direction set cannot - the mirror of dilution |
| automobile horsepower ~ engine size, n = 203 | .878 / 1.2×10⁻⁴ / .094 / .094 | divergence under heteroskedastic body-style scatter |
| rock perm ~ area + peri, n = 48 | .835 / 1×10⁻⁵ / 5.6×10⁻⁶ / .596 | A/D quiet, B/C reject; no truth to adjudicate |
| housing price ~ 4 features, n = 546 | .723 / .683 / ≈0 / - | D's leverage-degeneracy failure mode, live at n = 546, while A/B/C return verdicts |

**Scope.** What the bench supports is what S6 supports: compatibility at
scale (147/150 cells run across nine sources, zero silent failures),
failure registration on observational data, and mechanism-echo cells
whose directions match the synthetic map. Its per-gate rejection totals
(97/111/108/94 of 147) are compatibility statistics on single
realizations without ground truth - they rank nothing, and reading them
as performance would repeat exactly the leaderboard confusion the
benchmark exists to correct - A's observed non-rejections (the
shear sites, the road-load law, puromycin at n = 23), A's interaction and
dilution behavior (stormer, CCPP), D's structural fragility (all three
failure modes live). What it cannot support: rates, calibration, or
rankings - single realizations, no ground truth, and time-ordered
NWIS/Open-Meteo series whose permutation-independence assumption is an
approximation marked in the manifest. All four p-values per cell are
recorded locally (`real_data/smoke_results.csv`) with the manifest
(`real_data/clean/manifest.json`, SHA-256 per source) and replayable
fetch/prepare/smoke scripts.

## Appendix G - Audit and correction history

### G.1 Statistical, provenance, and presentation corrections

This ledger supersedes the mixed v2-v5 change notes in the supplied
working draft. The original PDF and Markdown are preserved
unchanged. Historical scripts and reports remain evidence of how
earlier summaries were produced, not automatic authorities for
the current numbers.

1. **Task pools and failures.** Record-level joins identify 7,500
   S2 statistic-sentinel failures for feynman-i.30.5 over the
   full grid. Revised power and joint silence use 18 H1 tasks;
   joint release additionally requires a finite far-field outcome,
   leaving 17 H1 tasks because i.26.2 lacks that outcome.
   Original generated-fit totals are not relabeled as eligible
   observations.

2. **Joint-release accounting.** Table 15 now defines all four
   counts and three conditional rates on the same eligible records.
   Its pooled counts are A=65,992, B=6,403, C=481, D=36,169 (shown as S, Sd, Ld, Da in R2).
   Gate-only silence is not substituted for complete-pipeline
   release. Comparisons between different conditional fractions
   no longer support a claim of danger-concentrated failure.

3. **S7 provenance and gate sets.** S7 is a reanalysis of existing
   S2/S5 records, not new sampling. The S5 arm uses A/E/F,
   not A-F. A paired A versus A/E/F comparison replaces the
   invalid cross-arm interpretation. Likelihood ratios and
   equal-prior posteriors have separate formulas. The S5 cache's
   weaker historical sentinel provenance is disclosed.

4. **SVD power and calibration.** Table 3 and the paired
   bootstrap use one 18-task pool and explicit decline handling.
   D-A effects are recomputed from the SVD data, not obtained by
   reversing old signs. The n=200 H0 task-level intervals are
   consistent between Section 4.1 and Appendix C. The manuscript
   no longer describes all C/D Wilson intervals as covering .05.
   Table 1 is also corrected: the 1,500 zero-noise declines per n
   are excluded, leaving 6,000 evaluated H0 records, rather than
   diluting the rejection rate with undefined tests.

5. **Claims, references, and artifacts.** The introduction and
   related work distinguish existing structural-recovery evaluation
   from this paper's diagnostic-gate experiment. References [6]
   and [11] have corrected authorship. Sample-size effects are
   not attributed solely to screening. Local script presence
   is not described as a public archive or a completed independent
   replication. The original scripts are not modified.

6. **What was checked in this revision.** The accompanying
   audit_size.py, audit_results.py, and audit_s7.py re-aggregate existing caches,
   record input hashes and task filters, and generate machine-readable
   evidence. build_pdf.py renders the corrected manuscript and
   figures. No new model fitting, public release, or engineering
   safety validation was performed. Remaining research limitations
   are in Section 6.

7. **Research presentation revision (R2).** R2 reorganizes the manuscript around detection, false-alarm cost,
and extrapolation-related joint release. The abstract and introduction
state the three experimental contributions; the discussion presents
method-selection lessons; Section 6 specifies two priority validation
experiments and the conditional extension to discovery-generated
candidates. Historical accounting corrections are retained in this
appendix. All numerical data rows in Tables 1-16 and the companion
engineering tables are preserved from R1. This revision adds no model
fits or completed validation experiments. Appendix H reports an
exploratory PySR counterpart outside the main tables.

8. **Review-response revision (R2, 2026-10-18).** The artifact is
public at https://github.com/leokynlang/form_check with data release
tag data-2026-10-18. Section 7 now points to it; reference [17]'s
licensing description is corrected; Section 4.10 adds a pooled
task-level threshold sweep; the abstract, introduction, discussion,
and conclusion state the task concentration of the 481 releases; and
Appendix H adds the exploratory PySR counterpart. No aggregate
(S, Sd, Ld, Da), table row, or figure in Tables 1-16 changed.

### G.2 Numerical implementation audit

The implementation audit compares the original pinv(X′X)
calculation, an SVD calculation using pinv(X), and statsmodels
HC3-Wald with the same F reference. The exog-HC3 check spans
all 13 H0 equations, stratified by equation, sample size, and
noise level (176 fits on identical input data). Maximum absolute
p-value difference between the SVD and statsmodels paths is
.0023; maximum relative statistic difference is .0038;
covariance Frobenius relative difference is at most 3×10⁻⁴.
These measurements establish close agreement on the tested fits.

The original path rejects at rate 1.00 on seven tasks whose
measured cond(X′X) is at least 4×10¹⁷, with measured
cond(X) from 6.5×10⁸ to 3.2×10¹⁰. On these small audit
samples, both stable paths have task rejection rates from 0
to .10. The two better-conditioned tasks ii.27.16 and ii.27.18,
with cond(X) near 10, show no corresponding failure.
A fitted-HC3 check on two tasks separately gives 8/120
rejections (.067, Wilson [.034, .126]).

Forming X′X squares the spectral condition number before
the pseudoinverse cutoff is applied. On the observed near-collinear
designs, discarded small singular directions alter the hat and
sandwich calculations. The audit attributes the large original-path
distortion to this numerical realization; it does not establish
a failure of the HC3 correction. The SVD path supplies the C/D
method comparisons reported in the main tables.

The earlier S2 H0 summaries for the original path were
C=.417-.422 and D=.563-.565 across sample sizes.
These historical values are withdrawn as evidence about HC3's
statistical performance. The implementation audit was a separate,
one-off comparison rather than a rolling software check on every
fit. Close implementation agreement and empirical size are
different validation questions.

## Appendix H - Exploratory PySR counterpart

This appendix reports an exploratory check that applies the four gates
to expressions emitted by a symbolic-regression search. It is separate
from the main benchmark: candidate generation is not frozen, the sample
is small, and the counts are compatibility evidence rather than rate
estimates. The data are in evidence/deep_analysis.json and
runs/pysr_gate/wide_svd/ in the data bundle.

**Protocol.** PySR searched the log-space training design with a
restricted arithmetic operator set (`+`, `-`) matching the power-law
hypothesis class. Each emitted expression was refitted as a prespecified
candidate and passed through gates A-D under the same native thresholds.
The run produced 360 candidates across the evaluated tasks.

**Results.** All 360 candidates: 87 jointly silent (24.2%). In the
R² ≥ 0.99 layer (183 candidates), 70 are jointly silent (38.3%), close
in magnitude to the 38.6% S2 four-gate figure at n=50 although the
candidate pools differ. In the far-field error >10% layer (142
candidates), only 3 are jointly silent (2.1%); the R² ≥ 0.99 and
far-field >10% intersection contains 2 candidates, both rejected by
gate D (p ≈ 5×10⁻²⁹) while A, B, and C are silent. Per-gate rejection
over all candidates: A 182 (50.6%), B 103 (28.6%), C 82 (22.8%),
D 189 (52.5%); A is the sole rejecting gate on 45 candidates and D
on 41.

**Scope.** These counts describe one search configuration and 360
candidates, of which only two meet the joint danger criterion. They do
not establish a pass rate for discovery algorithms, do not rank the
gates, and are not merged into Tables 1-16. Their purpose is to show
that the same diagnostic record can be attached to search-selected
expressions and that the main benchmark's failure channel has a
counterpart in this pilot.

