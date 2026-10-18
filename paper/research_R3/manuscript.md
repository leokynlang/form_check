# Silent Acceptance: A Benchmark of Specification Gates and Extrapolation Failures in Physical Model Assessment

**Revision status - R3, 2026-10-18. Working manuscript.**
This revision adds independently calibrated release policies (P0-1 v2)
and a fixed-domain, structural-group evaluation (P0-2). Historical
S1-S7 results retain their original populations and policies.
The proposed expansion of physical coverage is future work and supplies
no additional experimental outcomes. Appendices G and I record the
correction history and evidence boundary.

## Abstract

Specification gates challenge candidate equations, but a non-rejecting
record provides protection only against departures the selected policy
can detect. We evaluate power-law candidates using residual-permutation
curvature, RESET variants, quadratic-rival comparisons, and a fixed
radial-basis-function rival. We connect fit screening, independently
calibrated diagnostic policies, and error over external regions fixed
before training data are sampled.

The historical benchmark identifies directional blind spots and
numerical sensitivity; its single-probe joint-release result is retained
as a separate descriptive analysis. New validation first uses fresh data
for 20 independent reference/calibration pairs on the known equations.
An H0-only robust affine score removes artificial empirical-CDF tail
saturation without eliminating intrinsic statistic ties. The main
fixed-domain evaluation freezes structural groups, positive valid
regions, grids, weights, and outcomes before sampling. It separates
historically evaluated equations, source-derived physical tasks,
designed mechanisms, and reserve groups withheld from study outcomes;
calibration and
confirmation nulls test design transfer within the same monomial family.

The complete fixed-domain study contains 246,240 datasets across
three sample sizes and four noise levels. At 3% log-noise, the
quadratic AICc gate releases 41.59%, 21.88%, and 9.36% of
screened-dangerous wrong fits at n=50, 200, and 1000, versus
64.41%, 32.99%, and 15.65% for the A/E/F panel at the same
5% calibration budget. These primary dangerous cases come from
designed mechanisms. The fit screen intercepts all observed dangerous
physical-source fits, leaving no screened-dangerous case to estimate
the gates' conditional release rate in this six-source-task design.

Equal calibration budgets do not imply equal test false-rejection
rates. Group and calibration-round uncertainty remains conditional on
the registered library and realized test bank. The evidence supports
measuring residual release and correct-model losses together; it does
not establish universal gate rankings, broad physical-law coverage, or safety
of expressions selected by a discovery algorithm.

## 1 Introduction

A candidate equation is useful only to the extent that the available
evidence supports its intended use. Prediction error describes agreement
on sampled inputs; specification diagnostics ask whether selected
departures from an assumed form remain detectable. A workflow may fit
a candidate, require a high training score, run diagnostics, and release
models that trigger no alarm. The practical question is what errors
that complete rule removes and what errors it still releases.

The distinction is concrete in a controlled friction example. A power
law fitted to Prandtl-Karman data reaches R² = 0.99, passes
cross-validation and an information-criterion comparison, yet incurs
about -20% error one decade beyond training. The four specification
tests studied here detect that synthetic departure at power ≥ 0.99.
A fit-only failure can therefore be detectable, while other mismatch
directions survive the same checks.

We assess prespecified power-law candidates in log coordinates.
Permutation curvature compresses departure from log-linearity into
one quadratic direction; RESET examines powers of fitted values or
regressors; rival models add quadratic features or a frozen RBF
residual correction. Each construction has a sensitivity geometry.
Its usefulness also depends on false alarms, numerical stability,
the fit screen, and the policy by which several statistics are combined.

The historical S1-S7 suites map these effects using controlled
alternatives, physics-derived equations, family viability, rival
comparisons, and observational compatibility probes. They also expose
the limits of a cached endpoint analysis: the external coordinate
depends on the training sample, threshold crossings concentrate in
two equations, and individual reference thresholds do not control a
joint panel's independent-test cost. We retain those results as
historical evidence rather than treating them as the final policy
comparison.

Two prospective validation layers address these limitations.
P0-1 v2 uses wholly fresh samples and 20 independent reference/
calibration pairs, preserving the known equation identities and
legacy endpoint. It separates calibration sensitivity from task
composition and fixes an artificial tail ceiling in the initial score
mapping. P0-2 then freezes task grouping, valid external regions,
evaluation grids, and a region-level primary outcome before training
sampling. Historical tasks and algebraic near-twins remain in the
development pool. Newly evaluated physical-source tasks and designed
positive mechanisms are reported separately; variants and parameter
windows remain inside their structural groups.

Three findings organize the study. First, sensitivity follows the
departure direction, so a single diagnostic ranking does not describe
all alternatives. Second, joint policy calibration changes the value
of adding a gate: a panel's threshold must account for the complete
policy, and equal calibration targets still allow different
independent-test H0 costs. Third, gate silence and extrapolation
consequences must be linked on the same fits. Fixed-region release
and whole-workflow correct-model retention make that link explicit.

We make three contributions:

1. **A benchmark with explicit decision populations.** Form labels,
   fit screening, alarms, declines, release, and error outcomes are
   registered separately. Historical single-probe analyses and new
   fixed-region analyses keep their distinct denominators.

2. **A paired comparison of independently calibrated policies.**
   Single gates and prespecified panels share data and an H0-only
   joint calibration procedure. Repeated reference/calibration draws
   quantify policy sensitivity without counting repeated decisions
   on the shared test bank as independent observations.

3. **A fixed-domain structural-group assessment.** Regions, grids,
   weights, and truth are frozen before training draws. Near-twins
   do not cross evaluation roles. Structural-group and calibration-
   round resampling distinguish equation coverage from simulation
   volume, while physical-source and designed-mechanism evidence
   remain identifiable.

Symbolic-regression benchmarks evaluate predictive quality, expression
complexity, and symbolic or structural recovery [5, 6]. Our main
experimental object is the behavior of attached specification policies
under a fixed candidate family. Validation of expressions actually
selected by a discovery algorithm remains a separate requirement.
The proposed expansion of physical coverage in Section 6 has not
yet generated experimental outcomes.

Sections 3 and 4 define the protocols and results. Section 5 discusses
diagnostic and workflow implications; Section 6 states the remaining
coverage and selection limits. Section 7 and the appendices provide
provenance, correction history, and the boundary between completed
evidence and registered future work.

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

### 3.1 Design principles and historical policies

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
panels retain their measured error rates. These are historical reference-pool results. Sections 3.8-3.9
define the separate, independently calibrated validation policies.

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

Historical A-F configurations were fixed for S1-S7 - powers (2, 3); B = 999;
α = 0.05; the rival family and criteria of E/F; the grids below - and are
identical across their applicable suites. Within each comparison the
harness invokes its specified gates on the same datasets. The new
validation evaluates A-G together; it does not merge datasets across
historical suites or equate their different threshold policies.

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

### 3.8 P0-1 v2 - independently calibrated joint release

P0-1 v2 separates reference fitting, threshold calibration, and testing.
It uses the 32 historically evaluated equations (13 H0, 19 H1),
n ∈ {50, 200}, and Gaussian log-noise SD 0.03. Twenty independent
reference/calibration pairs contain 16,640 reference and 33,280
calibration datasets; one separate shared test bank contains 8,192
datasets (128 test replicates per equation and sample size). The two
sample sizes have disjoint test datasets (4,096 each); every fitted
policy is evaluated on the 4,096 datasets generated at its own n, so
the reported denominators are 128 records per equation per sample
size. All 58,112 datasets are newly generated. The historical
equations are known tasks; new seeds alone do not establish structural
generalization.

The first calibration pilot mapped statistics through a finite empirical
H0 CDF. Distinct values beyond the reference maximum then shared an
artificial ceiling, which could make a strict-budget panel unable to
alarm. The frozen v2 protocol replaces this map. For A-D,
u = -log10(max(p, float64.tiny)); for E/F/G, u is the raw rival statistic.
A valid H0 reference pool supplies each method's median center and
IQR/1.3489795 scale, falling back to normalized MAD and then population
SD. An unusable scale or invalid member declines the policy. Valid
standardized scores are not clipped to the reference range.
A's finite 999-permutation resolution and floating-point p-value ties
remain intrinsic limitations and are counted separately.

Five policies are prespecified: D; F; the A/B/C/D panel (AD);
the A/E/F panel (AEF); and G. A panel uses the maximum of its member
standardized scores. An independent H0 calibration pool fixes one
threshold for the complete policy at withholding budgets
1%, 5%, and 10%. Withholding includes alarms and declines; strict
score > threshold defines an alarm. If declines alone exceed the
integer calibration budget, the policy is infeasible rather than
silently treated as compliant. Release requires a training log-space
R² ≥ 0.99 screen, valid policy evaluation, and no alarm. Calibration
uses the screened H0 population. The same target budget is not a
guarantee of equal or bounded independent-test H0 withholding.

G is a fixed competitor: a 32-landmark RBF residual model, five-fold
cross-validation, training-fold-only scaling, squared bandwidth equal
to feature count, and ridge penalty 0.01 per observation. Its statistic
is the power-law model's CV error improvement under this competitor.
The kernel, complexity, and penalty are not tuned to H1 outcomes.
G is a finite fixed-rival diagnostic; we make no omnibus-consistency
claim.

The v2 outcome retains the sample-dependent legacy endpoint so that
the calibration change is distinguishable from a domain change.
The primary estimator is the equal mean of round-specific pooled-record
ratios. We report calibration-round SD and empirical percentile ranges
as sensitivity summaries, conditional paired task bootstrap intervals,
round-only intervals, and intervals crossing rounds with H0/H1-
stratified task weights (2,000 draws). All policies share resampling
weights; undefined ratios remain undefined. These intervals condition
on the known H0 calibration equation identities, realized calibration
rounds, and the test bank. They neither redraw the calibration equation
population nor include an inner record-level test bootstrap.

### 3.9 P0-2 - fixed external regions and structural-group holdout

P0-2 freezes a registry of 110 tasks before noisy training sampling:
32 historical tasks, six newly evaluated source equations with explicit
positive-domain restrictions, 48 designed positive-mechanism tasks,
and 24 new H0 designs. All historical source tasks and algebraic
near-twins belong to development. H1 structural groups are assigned
whole to development, confirmation, stress, or reserve. Global
prefactors, renamed variables, nonlinear-core variants, and parameter
windows cannot turn a known structure into an independent holdout.
The source inventory is visible during protocol design, so
confirmation means newly evaluated grouped tasks, not unseen formula
metadata drawn from an unknown population.

The confirmation pool has 38 H1 tasks in eight structural groups and
12 H0 tasks in four design groups. Six of the H1 tasks derive from
three physical-source groups: Gaussian profiles, a positive logarithmic
source formula, and hyperbolic response. The remaining 32 tasks are
controlled positive mechanisms, identified as such; they are not
presented as additional physical laws. Each designed mechanism has
a positive multiplicative amplitude input and four prespecified
training-window variants. Its confirmation mechanisms include erf,
atan, elliptic and Bessel functions, damped resonance, a nested
exponential, softplus, and sech; the full formulas and bounds are
in the frozen registry. Amplitude variation changes log-target
variance and therefore R² screening; the resulting challenge pool
is not a natural-frequency sample of physical mechanisms.
Eight tasks in two reserve
groups have no study samples, candidate fits, or gate outcomes.
Deterministic software/domain checks do not count as experimental
evaluation of those reserve groups.

Every H0 is a strict power law. Calibration uses 12 newly defined
log-uniform H0 tasks; confirmation uses 12 distinct exponent designs
sampled from Beta(2,2) in log coordinates. Groups correspond to the
one-to-four-feature sampling designs. This is held-out design transfer
within the monomial null family, not holdout of different null
functional families. The fixed R² screen can itself discard correct
low-signal models; conditional diagnostic cost and whole-workflow
H0 retention are both reported.

For each task, the registry states train and outer bounds, coordinate
scale (log, linear, or integer), registered domain restrictions, and an
analytic validity argument. Positive inputs and positive finite truth
are required throughout the complete registered outer region.
Constants and dimensional source expressions are retained.
Variables declared integer by the sampling schema use integer support;
that schema is not a complete audit of every textual physical type.
Stable expm1 evaluation
replaces cancellation-prone exp(t)-1 without changing the mathematical
source. For historical I.41.16, stable truth differs by at most 3.42%
from the cancellation-prone upstream callback and agrees with a
65-digit calculation to 1.68×10⁻¹⁵ relative error in preflight.
This task remains development, and the old experimental files are
unchanged. No task, domain, window, or point is adjusted
after observing fit quality, gates, or extrapolation error.

A deterministic external grid and its noiseless truth are saved before
training draws. The recipe specifies eight coordinates per axial
direction and 256 joint Sobol candidates plus outer corners, yielding
40-320 distinct external points after geometric filtering and merging.
Axial external points receive half the mass; joint
Sobol external points and outer corners receive half. Coordinates
inside the registered training box are excluded by a fixed geometric
rule. Duplicate coordinates merge with their weights summed.
The same grid is used for every n, noise level, and training replicate.
The primary endpoint is weighted mean pointwise relative error,
sum w_j |exp(pred_log_j - log truth_j)-1|; dangerous H1 means this
quantity >10%. Registered secondary endpoints are weighted p95 and
maximum relative error; mean-error thresholds 2%, 5%, and 20% are
sensitivity analyses. No denominator floor or finite-point deletion
is used. Error targets the positive noiseless function (the log-noise
location), not the lognormal arithmetic mean. These are finite-grid
summaries, not continuum error bounds. A legacy endpoint is evaluated
only when it lies inside registered valid support, and remains
compatibility-only.

The full frozen grid is n ∈ {50, 200, 1000} crossed with log-noise SD
{0.005, 0.01, 0.03, 0.05}. Stage 1 covers n=50 and 200 at SD 0.03;
stage 2 covers the other ten strata without rerunning stage 1.
There are 100 independent test datasets per task and stratum.
Ten independent reference/calibration pairs per stratum contain
32 reference and 64 calibration draws per calibration-null task.
P0-1 v2's score transformation, policies, and budgets are reused;
reference scales and thresholds are refitted only to these new H0
datasets. Formula registry, grouping, grids, seeds, and the complete
noise/sample-size grid remain frozen across stages.

Five named evaluation pools separate confirmation, physical-source
confirmation, designed-control confirmation, development, and stress.
The physical and designed subpools share confirmation H0 controls and
may share structural groups; they are not independent replications.
The stress threshold mechanism includes departures absent from some
training windows and does not determine the main policy ranking.

Primary points again average round-specific pooled-record ratios.
The primary pool is confirmation-all at the registered 5% calibration
budget. Resampling crosses the ten calibration rounds with
H0-design/H1-structure groups (2,000 paired draws), keeping all
variants in a group together and group weights shared across methods
and selected rounds. Conditional-group and round-only intervals are
also retained. An equal average of defined group-specific dangerous-
release ratios is secondary; its omitted groups are counted.
Zero risk denominators and undefined bootstrap draws are recorded,
never converted to zero risk. Calibration percentile ranges describe
ten-round sensitivity, not confidence coverage. These intervals remain
conditional on the fixed calibration-design pool and realized test
records; simulation volume cannot replace additional physical groups.

## 4 Results

Sections 4.1-4.11 retain historical S1-S7 results with their original
policies. Sections 4.12-4.13 supply the new validation evidence.
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
separate pool; the paired D-versus-F validation is reported separately
in Sections 4.12-4.13. The constructed alternatives below illustrate the
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

### 4.10 Historical joint release at a sample-dependent endpoint

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

### 4.11 Historical joint silence without the fit screen

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
Sections 4.12-4.13.

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

### 4.12 Fresh-data joint-policy validation and calibration uncertainty

All 600 P0-1 v2 policy configurations are feasible, giving 2,457,600
decisions on a shared independent test bank. AEF at n=50 and the
1% target has nonzero H1 alarms in all 20 rounds. This verifies removal
of the artificial reference-CDF ceiling. It does not remove statistic
ties from A's finite permutation resolution or floating-point p-values.

**Table 17 - P0-1 v2 at the registered 5% calibration budget. H0
withholding is alarm or decline conditional on R² screening.
Dangerous release is conditional on screened H1 with legacy endpoint
error >10%. Points average 20 calibration-round ratios. Brackets are
crossed task/round 95% percentile intervals.**

| n | Policy | Test H0 withheld | Dangerous H1 released | Crossed interval |
|---|---|---|---|---|
| 50 | D | 4.12% | 34.72% | [0.00%, 56.42%] |
| 50 | F | 4.68% | 7.98% | [0.00%, 48.15%] |
| 50 | AD | 4.40% | 21.35% | [0.00%, 48.15%] |
| 50 | AEF | 4.67% | 7.98% | [0.00%, 48.15%] |
| 50 | G | 5.46% | 7.98% | [0.00%, 44.60%] |
| 200 | D | 5.06% | 5.56% | [0.00%, 87.50%] |
| 200 | F | 4.79% | 5.56% | [0.00%, 100.00%] |
| 200 | AD | 4.69% | 4.41% | [0.00%, 77.81%] |
| 200 | AEF | 4.22% | 5.42% | [0.00%, 95.62%] |
| 200 | G | 5.36% | 4.83% | [0.00%, 86.25%] |

The independent-test costs are measured outcomes, despite common
calibration targets. The known-task evidence is too narrow for a
general ranking. Screened-dangerous denominators come from six
tasks at n=50 and four at n=200; AEF's dangerous releases still all
come from bonus.15. At n=50 F has 7.98% dangerous release in every
round. Its round-only interval is correspondingly degenerate, whereas
its crossed interval extends to 48.15%. Stable recalibration therefore
does not establish adequate equation coverage. At n=200, 19/2,000
risk bootstrap draws have no denominator and are recorded as undefined.
All v1 files and v2 outputs remain separately preserved.

### 4.13 Fixed-region validation across the complete registered grid

Stage 2 completes the ten strata omitted from stage 1, adding
205,200 independent datasets in 3,300 execution cells. The full
frozen n/noise grid therefore contains 246,240 datasets: 46,080
reference, 92,160 calibration, and 108,000 test. The primary
confirmation pool comprises 60,000 distinct test cases (14,400 H0,
45,600 H1); development and stress account for the other 48,000.
No stage-1 dataset, formula, domain, weight, endpoint, or policy
definition is replaced. The eight reserve tasks remain without
study sampling, fitting, or gate outcomes.

All 1,800 calibrated policy configurations are feasible. Their
16,200,000 decisions are stored in separate stage artifacts; repeated
decisions on one shared test bank are not additional independent
datasets. Independent stage audits reproduce fixed-grid error from saved
coefficients, H0 reference scales and joint thresholds, all decisions,
group event counts, and primary/equal-group/paired uncertainty.
A separate complete-grid audit confirms nonoverlapping data streams,
exact stratum/table integration, and preservation of the prior artifacts.
These checks use recorded A-G statistics; they do not independently
rerun every statistic from the original training observations.

**Table 18 - Complete P0-2 confirmation grid at the same 5% joint
calibration budget. Values are mean dangerous H1 release conditional
on R² ≥ 0.99 and fixed-region weighted mean error >10%. Ten
calibration rounds share the same independent test bank in each
stratum. G is the frozen RBF rival. The last column counts H1
structure groups contributing to this primary danger denominator.**

| n | Log-noise SD | D | F | AD | AEF | G | H1 groups |
|---|---|---|---|---|---|---|---|
| 50 | 0.5% | 10.04% | 7.47% | 12.54% | 17.54% | 11.73% | 6 |
| 50 | 1% | 19.91% | 15.96% | 23.14% | 32.62% | 22.68% | 6 |
| 50 | 3% | 50.35% | 41.59% | 55.75% | 64.41% | 53.67% | 7 |
| 50 | 5% | 72.39% | 53.99% | 76.45% | 79.92% | 74.09% | 7 |
| 200 | 0.5% | 0.55% | 0.39% | 0.89% | 2.55% | 1.30% | 6 |
| 200 | 1% | 8.40% | 6.33% | 10.89% | 14.25% | 10.22% | 6 |
| 200 | 3% | 26.51% | 21.88% | 31.51% | 32.99% | 31.12% | 7 |
| 200 | 5% | 43.79% | 35.30% | 49.66% | 51.93% | 48.26% | 7 |
| 1000 | 0.5% | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% | 6 |
| 1000 | 1% | 0.09% | 0.00% | 0.35% | 0.79% | 0.56% | 6 |
| 1000 | 3% | 11.45% | 9.36% | 13.59% | 15.65% | 13.30% | 6 |
| 1000 | 5% | 22.02% | 17.75% | 26.18% | 25.79% | 25.81% | 6 |

At the 3% log-noise level, F's mean dangerous release decreases
from 41.59% at n=50 to 21.88% at n=200 and 9.36% at n=1000.
AEF gives 64.41%, 32.99%, and 15.65%, respectively. At 5% noise,
F still releases 17.75% of screened-dangerous H1 at n=1000,
compared with 25.79% for AEF. All five policies have zero observed
primary dangerous release in the n=1000, 0.5%-noise stratum.
Across the complete grid, F has the lowest point in eleven strata,
and all methods tie at zero in the remaining stratum. This describes
this registered challenge library and budget; it is not a universal
ranking or a test-cost-matched dominance claim. Figure D.6 displays
the full-grid estimates with uncertainty.

The paired AEF-minus-F dangerous-release differences at 3% noise
are +22.82 [17.93, 29.64], +11.11 [8.67, 13.17], and
+6.29 [2.47, 9.00] percentage points for n=50, 200, and 1000.
Corresponding D-minus-F differences are +8.76 [5.23, 12.15],
+4.63 [3.35, 5.54], and +2.08 [0.75, 3.24] points.
These are crossed group/calibration-round 95% percentile intervals;
all 2,000 draws are defined for these comparisons. At n=1000
and 3% noise, F's individual interval is [3.54%, 13.76%]
and AEF's is [6.27%, 22.23%]. Intervals are pointwise rather
than simultaneous guarantees across strata. They remain conditional
on the fixed calibration-design pool and realized test records.
The limited conclusion is that adding panel members does not
guarantee lower dangerous release at a common calibration budget
on these controlled mechanisms.

Each noise level shows decreasing release points as n increases,
but different n values use independent training samples and can have
different screened-dangerous composition. The point pattern is not
a nested-sample causal estimate of sample-size benefit. Six or seven
H1 groups contribute primary danger in each stratum; all primary
screened-dangerous cases come from designed mechanisms. Their
distinct parameter windows remain grouped for uncertainty.

**Table 19 - Correct-model cost and physical-source denominator
coverage for the same grid. H0 columns use the confirmation null
designs. F/AEF withholding is conditional on the fit screen; F
retention uses all H0 test cases. Physical Sd counts screened-dangerous
physical-source H1 cases: the primary mean-error criterion and the
registered p95 secondary criterion, both >10%. These are event
counts, not physical risk rates.**

| n | Log-noise SD | H0 screen | F withheld | AEF withheld | F retained | Physical Sd | p95 Sd |
|---|---|---|---|---|---|---|---|
| 50 | 0.5% | 100.00% | 4.47% | 4.51% | 95.53% | 0 | 100 |
| 50 | 1% | 97.67% | 3.46% | 5.21% | 94.28% | 0 | 100 |
| 50 | 3% | 70.67% | 5.24% | 4.80% | 66.97% | 0 | 100 |
| 50 | 5% | 49.92% | 3.52% | 3.67% | 48.16% | 0 | 102 |
| 200 | 0.5% | 100.00% | 5.03% | 4.46% | 94.97% | 0 | 100 |
| 200 | 1% | 99.33% | 4.77% | 4.98% | 94.59% | 0 | 100 |
| 200 | 3% | 69.58% | 3.95% | 5.58% | 66.83% | 0 | 100 |
| 200 | 5% | 48.08% | 3.80% | 6.00% | 46.26% | 0 | 100 |
| 1000 | 0.5% | 100.00% | 7.13% | 5.42% | 92.87% | 0 | 100 |
| 1000 | 1% | 99.92% | 4.71% | 5.05% | 95.21% | 0 | 100 |
| 1000 | 3% | 68.67% | 4.25% | 3.25% | 65.75% | 0 | 100 |
| 1000 | 5% | 49.50% | 4.71% | 7.76% | 47.17% | 0 | 100 |

Correct-model loss is substantial at larger noise, even before a gate
acts. At 3% noise the H0 screen passes 70.67%, 69.58%, and 68.67%
for n=50, 200, and 1000; F's whole-workflow retention is 66.97%,
66.83%, and 65.75%. At 5% noise the screen passes only
48.08%-49.92% across n. Conditional H0 withholding also varies
despite identical calibration targets: F ranges from 3.46% to
7.13% and AEF from 3.25% to 7.76%. In the low-noise n=1000
stratum, both release no observed primary dangerous H1 but F's
test H0 cost is 7.13% versus AEF's 5.42%. The measured costs
are integral to the comparison. The 5% calibration budget does
not imply 95% overall retention or an independent-test size
guarantee. Figure D.7 separates this loss from physical coverage.

The physical-source pool's primary conditional rate remains
undefined in all twelve strata because its screened-dangerous
denominator is zero. The complete workflow nevertheless releases
0 of 3,600 observed primary dangerous physical-source fits across
the grid, an overall secondary result: all three Gaussian tasks'
1,200 fits apiece are dangerous and all fail the R² screen.
The other three source tasks' 3,600 fits all clear the screen
and none meets the primary danger criterion. Thus the fit filter
intercepts the observed dangerous cases; this pool supplies no
evidence on how the gates handle high-fitting but externally
dangerous physical candidates. It is neither a post-screen zero-risk
estimate nor a physical-population safety statement.

The registered physical p95 secondary criterion has 100
screened-dangerous cases in eleven strata, all from I.44.4.
At n=50 and 5% noise it has 102, including two II.35.21 cases
from a second group. At 3% noise, F's p95-dangerous release is
97.0%, 19.2%, and 0% across n=50, 200, and 1000; AEF's is
95.4%, 47.3%, and 0%. These are gate/workflow consequences
for a secondary endpoint supported by one source group at that
noise level. They do not promote that endpoint to the primary
criterion or establish broad physical performance. At n=1000
the p95 intervals degenerate to zero on the fixed test bank, and
only 1,393/2,000 crossed resamples have a defined denominator.
Undefined draws remain undefined. Such single-group intervals
do not acquire physical-population precision from the number of
simulation replicates.

All data-generation, candidate-fit, and fixed-grid outcomes are valid.
B and C each decline as rank_deficient on the same 513
development-test datasets; A, D, E, F, and G have no recorded
failures. These are 513 affected datasets, not 1,026 independent
failures. The declines occur in II.11.28 and III.12.43 and do
not affect the confirmation pool. They remain in the registered
development decision accounting, including AD withholding; none is
deleted or regenerated. Full operational, threshold, endpoint, and
danger-cutoff tables are retained in the stage and complete-grid
artifacts.


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
augmentation geometry and model-comparison alternatives. The new paired common-budget comparison in Section 4.13 evaluates
D and F on the same records, while preserving its task-library and
test-cost qualifications.

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

### 5.3 Measure useful rejection, correct-model loss, and residual release

The historical endpoint analysis establishes an observable residual
channel: 481 of 6,403 screened-dangerous H1 fits pass its native
four-gate policy. It also shows that those releases concentrate in
two equations. P0-1 v2 makes independent joint calibration explicit,
and P0-2 fixes the external coordinates and structural-group units.
These layers answer related questions with different pools and
outcomes; their percentages cannot be pooled into a single risk rate.

Adding members changes a panel's maximum-score null distribution and
its jointly calibrated threshold. Consequently, extra members can
reduce sensitivity to some departures even though a panel at unchanged
member thresholds would alarm more often. The empirical comparison
must therefore include the target budget, realized independent-test
H0 withholding, and the H1 consequence on the same fits. A lower
dangerous-release point at one budget does not establish dominance
at identical test cost or on an arbitrary new population.

The R² screen deserves the same scrutiny as the diagnostics. The
calibration budget concerns diagnostic withholding among screened
H0 candidates; it does not guarantee 95% retention of all correct
models. Whole-workflow H0 retention, fit-screen passage, alarms,
declines, and dangerous release are separate outputs. An undefined
screened-dangerous denominator is also a result about identifiability:
it neither estimates the post-screen gate risk nor establishes
interception by those gates. Overall dangerous release can still be
measured against the pre-screen dangerous population, as in the
physical-source results in Section 4.13.

Region-level error connects the decision rule to intended use more
directly than one moving endpoint. The registered weighted mean,
p95, and maximum still summarize a finite grid, use a specified
relative-error tolerance, and depend on the registered domain.
No diagnostic record supplies a distribution-free extrapolation bound.
Interpretation must retain which physical or designed groups actually
contribute dangerous screened candidates.

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

### 6.1 Scope of the completed evidence

The main experiments assess a fixed power-law candidate family on
positive data. They do not evaluate the expression distribution of a
search algorithm. S4 extends the coordinate construction to four
linearizable families; the new policies add a quadratic rival and one
frozen RBF rival. Strongly dependent errors, additive noise, other
rivals, and application-specific tolerances need their own evidence.
Residual permutation remains an empirical approximation under the
evaluated conditions.

The historical S2 detection and gate-only pools contain 18 eligible
H1 equations, and its endpoint-release pool contains 17. That outcome
uses a sample-dependent first-coordinate probe and remains labeled
as historical. P0-1 v2 retains known equation identities while changing
sampling and calibration. P0-2 fixes domains and keeps structural
groups together, but its primary confirmation pool deliberately
contains both source-derived physical tasks and designed mechanisms.
The number of formulas or simulation replicates is not the number
of independent physical mechanisms. The additional amplitude input
in designed controls also changes the population that clears R²;
their screened-dangerous frequency is not a physical prevalence estimate.

Source sampling metadata also require a physical interpretation audit.
A positive finite mathematical region does not by itself establish
all physical admissibility conditions or unit/type semantics. The
source description of I.39.22 treats n as an integer, whereas its
benchmark sampler is continuous. This is a historical development
task, and its frozen P0-2 records are retained under the registered
sampling interpretation. It is not one of the six physical-source
confirmation tasks. The future source census flags such contract
conflicts before any new physical sampling; the current study does
not claim that every legacy domain has completed that broader audit.

The physical-source pool has no screened-dangerous primary H1 record
in any of the twelve registered strata. Its corresponding conditional
release rate is undefined throughout. The complete registered workflow
nevertheless releases 0 of 3,600 observed dangerous physical-source fits
across the grid: all are stopped by the R² screen. This overall
secondary result describes successful filtering of these particular
cases, with no gate interception required. The missing evidence concerns
high-fitting physical candidates that are dangerous over the external
region; it does not make every physical workflow outcome unmeasurable
or establish physical-population safety. The p95 secondary endpoint contributes
primarily one source group, with two additional records from a second
group in one stratum. It remains secondary and cannot establish a
broad physical policy ranking. The stronger multi-group primary
comparison concerns the designed mechanisms.

All H0 tasks share the monomial family; their held-out grouping tests
sampling-design transfer. Test H0 withholding is measured rather
than guaranteed by the reference target. Intervals condition on the
calibration-design identities, empirical calibration-round bank, and
realized test data; they omit inner test-record resampling and
uncertainty from drawing a different physical library. Ten rounds
give an empirical sensitivity description, not a precise estimate
of the full distribution of possible recalibrations. Finite-grid
mean or maximum error is not a continuum guarantee, and training
samples at different n are independent rather than nested.

### 6.2 Outcome-independent plan to expand physical coverage

A separate physical-coverage design, p03-physical-coverage-source-census-v1,
was developed after P0-2 stage-1 results were available, while the
registered full grid was being completed. It is not a retrospective
preregistration of the current results. Its admission procedure reads the complete pinned
120-equation source catalog, source sampling and constraints, structural
exposure identities, and validity certificates. It forbids use of P0-1/
P0-2 R² values, gate statistics, error or danger labels, screen-passing
rates, or method rankings to admit a formula or choose its domain.
The source count is a census, not a selected set of successful tasks.
The scientific dimension limit is eight inputs; eight is also the
maximum in this catalog, so that limit excludes no current source.

One common domain rule conditions the original sampling support on
positive inputs and takes its central 5th-95th percentile box in the
original sampling coordinate. Original integer support and sign branches
are preserved before conditioning. The complete outer box must satisfy
source constraints and strictly positive finite signed truth. Constants,
coordinates, supplied unit metadata, and the source expression remain
intact; absolute values, sign flips, output shifts, hand-picked positive
lobes, and equation-specific domain rescue are forbidden. Missing unit
or source-contract evidence is recorded explicitly, not inferred from
variable names. Failure to establish a whole-box certificate is an
unresolved obligation, not proof that the equation has no valid physical
region anywhere.

The census records 86 mathematically admissible equations under the
registered box rule (35 H0, 51 H1),
31 with unresolved validity or source obligations, two with nonpositive
target throughout the box, and one without original positive-input
support. Mathematical admission is separate from executable evaluation:
eight eligible equations have more than four inputs and require adapter
validation; two Gaussian laws need stable log-truth treatment rather
than exclusion for floating-point underflow. These counts do not imply
that physical/unit provenance and numerical checks are complete.
Physical source-contract review remains independently pending for
all 120 rows. A uniform textual-type scan finds seven documented
integer-to-real relaxations and flags I.39.22 as
source_type_review_pending because no such relaxation agreement
resolves its integer description and continuous sampler.

The primary grouping conservatively merges variable renamings,
global and separable monomial prefactors, fixed output powers or
reciprocals, source variants, and registered mechanism families.
Exposure propagates through connected components. It yields 13
components across the catalog; the eligible H1 equations occupy only
four already exposed groups under this conservative rule. The
grouping is coarser than P0-2's structural grouping. Counts from the
two protocols describe their registered resolutions, not an objective
number of physical mechanisms, and are not directly comparable.
There are therefore zero
undisputed prospective eligible H1 groups in this catalog under the
primary grouping. A finer 14-component description leaves bonus.16
as one disputed combination candidate, retained in conservative
development. Neither counting more source equations nor treating
that disputed case as many mechanisms supplies independent physical
confirmation.

The proposed future factorial crosses centered training widths
25%, 50%, and 75% of each fixed source-coordinate outer span with
all four log-noise SDs {0.005, 0.01, 0.03, 0.05} and n ∈ {200, 1000},
using 100 test replicates per complete cell. Every admitted group and
variant uses the same design; difficult screens and numerical declines
remain in the operational denominators. The sample-size choice is a
structural feasibility rule: at the maximum eight inputs, a full
quadratic has 45 columns and five-fold training at n=200 has 160 rows.
It is not a result-based power selection. R² ≥ 0.99 is only a
downstream release screen, never a source-admission or stopping rule.

This is a metadata-level planning registration with status
pre_execution_registration_incomplete, not an executable experimental
freeze. Exact task roles, source/physical certificates, stable log-truth
and higher-dimensional adapters, width-specific fixed grids/weights,
dimension-matched monomial calibration controls, source snapshots,
and new seed namespaces must be completed before sampling. A new
external physical-source catalog has not yet been chosen; its complete
census and group-exposure review are prerequisites for multi-mechanism
prospective confirmation. The present registration has generated
zero new physical study samples, candidate fits, or diagnostic outcomes.
It contributes a reviewable inclusion/exclusion procedure, not a
completed extension of the physical empirical evidence.

### 6.3 Extension to discovered candidates and engineering decisions

The next algorithmic test must evaluate expressions actually selected
by a discovery workflow. Candidate generation and selection must be
fixed before independent diagnostic and external prediction evaluation.
The null attached to the diagnostic must correspond to the selected
candidate; a power-family diagnostic applied after search is not a
general candidate-specific validity certificate. The exploratory PySR
counterpart in Appendix H does not replace that experiment.

The observational suites show software applicability and diagnostic
compatibility patterns, while their true equations and external loss
are unknown. Engineering recommendations require reference outcomes
over an intended operating domain, relevant error dependence, and an
application loss tolerance. Those requirements go beyond detecting
one form of statistical misspecification.

## 7 Reproducibility statement

The local source tree D:/form_check contains generators, recorded
configurations, results, failure labels, and earlier analyses.
The R2 manuscript, PDF, figures, and audit artifacts are preserved.
The R3 manuscript integrates new evidence without replacing the
historical experimental files. Appendix A reports original S1-S6
fitted-run totals; S7 reuses existing records. New P0-1/P0-2 dataset
totals are additional and are not mixed into those historical counts.

Historical aggregation evidence is retained in the R2 bundle, including
evidence/size.json, results.json, s7.json, task_intervals.json, and the
record-level provenance audits. Figures D.1-D.5 describe those
historical runs; Figure D.5's timing cannot be used as a timing
comparison for the new RBF policy or updated inference paths.

P0-1 v1 and v2 reside in runs/p01_pilot_v1 and runs/p01_v2.
The v2 frozen protocol, manifest, source snapshot, raw_policy.parquet,
reference scales, policy thresholds, per-task counts, round summaries,
and audit.json make fresh-sample calibration and uncertainty
reconstructible. Source code is under experiments/release_validation;
validation/audit_release_policy_v2.py independently recomputes policy
metrics and uncertainty without importing their experimental
implementation.

P0-2 resides in runs/p02_fixed_domain. protocol_frozen.json,
task_registry.json, group_isolation.json, grid_manifest.json, and
grids/ record the pre-sample registry and fixed truth endpoints;
manifest and source snapshots bind the run to its implementation.
Stage 1 and stage 2 retain their separate raw records, checkpoints,
decisions, thresholds, and summaries. The complete-grid aggregation
combines the two nonoverlapping stage grids and does not replace
either stage. The independent fixed-domain auditor reconstructs grid
error from saved coefficients/truth and policy/count/uncertainty
outputs from recorded raw diagnostic statistics; it does not
independently recompute every A-G statistic from original training
arrays. Registered secondary summaries retain their own endpoints
and are not substituted for the primary outcome.

The proposed physical expansion is documented in
docs/PHYSICAL_COVERAGE_PROTOCOL.md and
experiments/physical_coverage_v1/protocol.json. The complete
120-row ledger, structural components, rules snapshot, and review
hashes are retained under runs/p03_physical_protocol. Its manifest
states incomplete pre-execution registration and zero generated
study outcomes. An independent anti-selection/source-fidelity review
passes for this snapshot; it checks the full census, uniform domain
rule, and exposure metadata, without certifying all interval bounds
or physical/unit/type contracts. Hashing this planning snapshot does not certify
that the future study is runnable or that its pending physical
coverage has been achieved.

Preservation audits compare full file-tree hashes for existing runs.
Seed and recorded dataset-hash checks separate the new studies from
the preserved P0-1 runs and from one another; source-file integrity
is distinct from reconstructing every historical dataset lacking a
saved data hash. Reproduction of fitted experiments is also distinct
from reaggregation of saved raw statistics.

The original R2 artifact identifies repository
https://github.com/leokynlang/form_check and historical data tag
data-2026-10-18. This revision does not claim that its new local runs
have been published in that historical release. The local manifests
are the provenance record for the new evidence. Third-party
redistribution decisions remain governed by data_release/
THIRD_PARTY_DATA.md and dataset_sources.json; copyrighted NIST
WebBook SRD 69 source data remain source-download-only. No new
external references or external empirical claims are added here.

The R3 bundle distributed with this manuscript contains the summary
evidence for the new studies: revision_facts.json, per-task count and
policy-threshold tables, round-level calibration summaries, paired
and interval tables, and the independent audit outputs. The complete
P0-1 v2 and P0-2 runs are larger and remain in runs/p01_v2 and
runs/p02_fixed_domain. In P0-1 v2, decisions.parquet and
raw_policy.parquet alone total about 41 MB, and the run also stores
protocol_frozen.json, manifest.json, source_snapshot.json, the per-cell
files, and audit.json. P0-2 keeps two nonoverlapping stage trees, a
1.8 MB grid_manifest.json, and freeze_independent_audit.json. These
record-level files, including the 58,112 P0-1 v2 datasets and the
246,240 P0-2 datasets, are available as a frozen copy on request to
the author. Shipping only the summary tables is a size decision, not
a claim that the new runs are already public.

## 8 Conclusion

Specification policies can remove dangerous candidates while leaving
an interpretable residual release channel. The historical benchmark
maps directional blind spots, records false alarms, and shows why
stable numerical inference matters. Its native-panel single-probe
result remains informative but concentrated in two equations.

New validation separates H0 reference fitting, independent threshold
calibration, and testing. The v2 affine score removes artificial
CDF tail saturation, and repeated recalibration reveals why a stable
decision rate can coexist with wide task-composition uncertainty.
The fixed-domain study then evaluates the same fitted models over
pre-sample external grids and holds related formulas inside structural
groups.

The complete registered n/noise grid shows substantial residual
release on the controlled mechanisms even at n=1000 and larger
noise. F has the lowest primary release point in eleven strata,
with an all-policy zero tie in the remaining stratum, while actual
independent-test H0 costs differ. A larger panel therefore does
not guarantee lower dangerous release at a common calibration
budget in this library.

The physical-source result has a different boundary. Its observed
dangerous fits are all stopped by R² screening, so overall dangerous
release is zero in those cases; no screened-dangerous case remains
to identify the gates' conditional physical performance. Correct-model
losses are substantial at the larger noise levels, largely through
the fit screen. The new source-census protocol
records this coverage problem without selecting formulas or
narrowing domains from observed performance, but its mathematical
admission ledger and pending source contracts are not completed
multi-mechanism physical confirmation.

The practical contribution is a reproducible way to compare
complete decision policies using actual H0 costs, correct-model
retention, dangerous-release denominators, and group/calibration
uncertainty together. Broader physical coverage and independent
validation after expression selection remain necessary before
translating these controlled comparisons into general scientific
discovery or engineering recommendations.

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

Sections G.1-G.2 retain the historical correction record through R2.
Statements about work performed in those revisions do not describe
the additional fitted experiments in R3; the R3 additions are listed
in Appendix I and Sections 3.8-3.9 and 4.12-4.13.

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


## Appendix I - New validation and proposed-work boundary

P0-1 v2 is a known-equation fresh-sample experiment; P0-2 is a
prospectively fixed-domain and grouped-task experiment. The latter's
physical-source and designed-mechanism labels remain separate in
the artifact and in interpretation. A newly executed source formula
does not become an independent mechanism merely because its
coefficient, variable names, or positive restriction differs.

Stage 1's two strata and stage 2's ten strata are nonoverlapping
parts of the same frozen 12-cell grid. Their joint raw record count
is 246,240, with 108,000 independent test datasets. Reusing each
test bank across calibration rounds changes decisions, not the
independent sample size. The recorded B/C rank declines remain in
development accounting; historical R2 tables and all prior run
artifacts are preserved. The complete-grid summaries are derived
outputs and do not replace either stage's raw results.

The proposed expansion of physical coverage is a separate future
study. Its written design supplies no empirical evidence that more
physical structures pass a fit screen, become dangerous externally,
or change the relative gate ranking. Any later software or eligibility
change must preserve this record and state its timing relative to
outcome access.
