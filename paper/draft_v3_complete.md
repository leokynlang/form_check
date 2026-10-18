# Silent Acceptance: Goodness-of-Fit Alone Cannot Certify Functional Form — Blind Spots of Fit-Based Acceptance and Specification Gates in Physical Model Discovery

**Complete draft v3 — working manuscript, not submitted.** Supersedes
`draft_v2_complete.md`. Changes from v2: two suites are added — **S5, the
gate-comparison suite** (the benchmark's most obvious attack, "why not just
compare against a nested rival?", answered with data; `runs/gate_cmp`;
20,400 datasets) and **S6, a real-data suite** on World Bank open macro
data (`runs/wb_real`; 165 datasets, external validity + a live dilution
signature) — and the paper is **repositioned** accordingly: no detector is
claimed as the gate; the deliverable is a measured failure-mode map plus a
deployment ordering. All v2 numbers are retained as audited; the change
notes, including both revision rounds, are Appendix G.

## Abstract

Symbolic regression and physical model discovery pipelines certify
candidate models with goodness-of-fit metrics — R², cross-validation
error, information criteria, coefficient standard errors. This paper
measures what those metrics can and cannot certify, in both
directions, and maps where the classical specification gates share
blind spots. The claims are deliberately bounded: fit quality cannot
*by itself* certify a functional form — within the single R² ≥ 0.99
stratum, true power laws extrapolate at 0.1–0.2% error while wrong
forms reach 100% >10% error on our synthetic bench — a design
property of that bench, whose deviations are built to be detectable,
not an estimate — and 9.7% of the wrong forms screened on the pool
synthesized from real Feynman equations are dangerous; and the
evaluated specification tests (RESET and
variants, nested-rival comparison) remain jointly blind to some
form mismatches — wrong forms passing an R² ≥ 0.99 screen are
jointly missed by all four gates at 0.67% overall, and among
*dangerous* wrong forms (extrapolation error >10%) the joint
release rate is 6.9%. With that screen removed the same four gates
are jointly silent on **37.8%** of wrong forms drawn from real
equations at n = 50 — roughly two orders of magnitude above the
synthetic bench's 0.15–0.22% on the same unscreened basis — and
adding selection-side gates raises rather than lowers the calibrated
figure (38.6% → 42.2%, against 90.5% silence on true forms; a silent
certificate therefore carries only a 2.2× likelihood ratio). We do
not claim, and our data contradicts,
stronger readings: fit screening combined with the gates carries
real information (a full-green certificate is 65× less likely for
a dangerous wrong form than for a true one: .013 against .845); the failure is
non-certification and non-detection of specific mismatch
directions, not statistical independence.

On the synthetic bench, within the single R² ≥ 0.99 stratum, true
power laws extrapolate at ≈0.1–0.2% error while wrong forms sit at
100% >10% error: the same certificate, opposite function truth. The
paper's headline quantity is the *joint* event on the same models —
a dangerous wrong form (extrapolation error >10%) that passes the
R² ≥ 0.99 screen AND is jointly missed by all four gates: 0.67% of
screened wrong forms (481/71,992), but 6.9% of dangerous ones; and
the release rate of dangerous forms falls from 7.8% (n = 50) to
0.2% (n = 1000) counting all dangerous fits, and from 13.2% to 1.2%
counting screened dangerous models — a fall driven by the fit screen
tightening at scale, not by any single gate (Section 4.10b). A
full-green certificate carries real information (a dangerous wrong
form is 65× less likely to hold one than a true form: .013 against
.845); what it cannot do is *exclude* the failure mode.

We introduce the falsification-side benchmark that quantifies this.
Six suites total 1,480 grid cells and 508,825 fitted datasets, with
472,500 further held-out evaluations and a quality re-evaluation
over 297,000 of these fits: (S1) 19 controlled synthetic
scenarios; (S2) 32 tasks from the SRSD-Feynman registry under official
sampling configurations; (S3) 30 Feynman equations stratified by
structural family; (S4) a 188-cell family-viability battery extending
the same machinery from rejection to routing; (S5) a gate comparison
placing the benchmark's best test against the "obvious alternative"
it had never faced — a nested quadratic rival adjudicated by 5-fold
cross-validation or AICc; and (S6) six macroeconomic relations from
World Bank open data across four years and ~190 economies per cell.
A seventh suite (S7) re-samples S2's 18 non-monomial equations on an
independent seed grid — 108,000 further fits — to measure the panel's
joint-silence rate with the fit screen removed; its counts are
reported separately and are not folded into the 508,825 above.
All methods run under one protocol: a common significance level; power
compared only at matched empirical false-positive rate; preregistered
labels; explicit registration of every degenerate or rank-deficient
evaluation.

Five findings stand out. **First**, the phenomenon is real and silent:
a power law fitted to Prandtl–Karman friction data certifies itself at
R² = 0.99 through every metric a pipeline computes, while extrapolating
at −20% one decade out. **Second**, the benchmark caught — and its revision retracted — its own
most striking claim: our harness's HC3-Wald path rejects at rate 1.00 on
5 of 13 physics-derived tasks whose built designs are near-collinear
(pooled 42–56% false alarms, flat in n), but an independent SVD-based
implementation (statsmodels) is nominal-compatible (.067, Wilson
[.034,.126]) on the same data and agrees with ours to three digits
elsewhere: the failure is a silent
numerical pathology of pinv-squared-conditioning inference paths, not of
the robustness correction — a deployment hazard discovery code should
test for. **Third**, no test is
uniformly best; rankings invert under matched false-positive rates and
again across task families, and each method's strength sits on another's
blind spot. **Fourth**, the battery extends rejection to routing after
its null is rebuilt per family — null behavior is coordinate-dependent
(the residual permutation is empirically calibrated in log-log
coordinates on every pool we measured; 30–34% false-escalation in
rational coordinates until a wild bootstrap repairs it). **Fifth** — the
benchmark's most uncomfortable results — the nested quadratic rival is
the better default detector wherever a reference pool exists, and the
naive OR union of the permutation gate with the CV rival, previously the
best-recall cell, owes that advantage to a larger joint false-positive
budget: recalibrated to a 5% joint budget, the union matches but does
not beat AICc-quad. On real macroeconomic data the gates agree where
signal is clear, diverge where the map predicts weakness, and show an
observational pattern consistent with the synthetic dilution
signature; S6's six relations carry no ground truth and license
signature-consistency claims only, not rates.
The deployment ordering: default AICc-quad when a reference pool is
available; the permutation gate when a self-calibrating p-value is
required; recall-critical deployments should budget for AICc-quad, not
for an uncalibrated union. We release all generators, per-cell results,
failure registrations, and content-hashed configurations under fixed
seed rules.

## 1 Introduction

A regression model can be simultaneously excellent on every certification
metric a discovery pipeline computes and wrong in its functional form. A
power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes
k-fold cross-validation, wins its information-criterion comparison, and
reports small standard errors — while extrapolating with −20% error one
decade outside the training domain. This is not hypothetical: the failure
is reproduced as a controlled scenario in our laboratory notebook and
appears in this benchmark as scenario S1 (Section 3.3), where all four
specification tests detect it at power ≥ 0.99. Fit quality and form
correctness are different axes, and current discovery practice measures
only the first.

This paper sharpens that observation into falsifiable, bounded claims
and measures both directions. Claim one — fit quality cannot by
itself certify a functional form. On our synthetic bench (Section
4.10), within the single R² ≥ 0.99 stratum, true power laws
extrapolate at 0.1–0.2% error while wrong forms sit at 100% >10%
error: the same certificate, opposite function truth. Claim two —
the evaluated specification gates remain jointly blind to some form
mismatches. The joint event on the same models (Section 4.10):
a dangerous wrong form (extrapolation error >10%) that passes the
R² ≥ 0.99 screen and is missed by all four gates occurs for 0.67%
of screened wrong forms but 6.9% of dangerous ones; the extreme
case, equation ii.11.28, is jointly missed at 83–84% (n = 50/200)
despite fitting at R² ≈ 0.002. Both numbers are necessary-negation
results, one clean counterexample each would establish the logical
point; what 297,000 re-evaluated fits add is the *extent* — where
the blind spots sit (odd-symmetric and interaction mismatch
directions), how often they fire, and how they scale with sample
size. We do not claim the absence of a statistical relationship:
fit screening plus the gates carries real information (a dangerous
wrong form is 65× less likely to earn a full-green certificate than
a true one), and the alignment diagnostic corr(ŷ², residual)
predicts gate behavior at Spearman +0.47 to +0.76. The failure mode
is specific: no combination of these tools *excludes* being wrong
in a direction they do not span — certificates of correctness
cannot be built from fit magnitude alone.

The statistics literature has owned the form axis for over fifty years.
Ramsey's RESET [1] tests whether powers of the fitted values carry
additional signal; built-regressor variants use powers of the regressors;
White [2] and MacKinnon–White [3] provide heteroskedasticity-robust
inference. And the workflow alternative that requires no new statistics at
all — fit a nested rival family and let cross-validation or an
information criterion arbitrate — has been standard practice in applied
regression for as long. Yet none of this is measured by symbolic
regression benchmarking. SRBench [5] measures whether discovery methods
find the right expression across 252 problems; SRSD [6] redesigns
Feynman-derived tasks with realistic sampling ranges. In neither framework
does anything measure whether an attached falsification layer would
reject a wrong form. The asymmetry has practical consequences: a discovery
pipeline with no working falsification layer silently accepts misfit
forms, and its downstream use — extrapolation, engineering design,
scientific claim-making — inherits the error silently.

Measuring the falsification side faithfully requires resolving a protocol
problem that, to our knowledge, no prior work addresses for this test
family: the candidate gates reject on incommensurable bases. A
bootstrap-t threshold ("alarm when |t| > 3"), an asymptotic F p-value
("reject when p < 0.05"), and a cross-validated RMSE difference ("escalate
when the rival wins") measure different things; comparing them directly
measures threshold placement, not detection. Our first pilot compared a
native |t| > 3 alarm against p < 0.05 and produced inverted rankings that
vanished under a common basis. The protocol below (Section 3.1) rests on
three commitments: one decision basis per comparison; power compared only
at a matched empirical false-positive rate estimated on true-monomial
tasks; and explicit registration of every failure.

This paper asks its benchmark the three questions a reviewer would ask,
in the order they would ask them. (i) *Do the classical tests work at all
on physical data distributions?* — Sections 4.1–4.5: calibration, power,
blind zones, extrapolation semantics. (ii) *A rejection terminates
inquiry — then what?* — Section 4.6: a family-viability battery that
routes instead of refusing, once its null is rebuilt per coordinates. (iii)
*Why not just fit a nested rival and compare?* — Section 4.7: because
until it is measured, that question has no answer; measured, it wins, and
the benchmark says so. A real-data suite (Section 4.8) then asks the
external-validity question: does the synthetic failure-mode map predict
anything about live observational data? It does, including one live
reproduction of a synthetic blind spot.

The question this benchmark answers is narrow and falsifiable: **when the
ground truth is known and the candidate form is wrong, which gates detect
the misspecification — at what false-alarm cost, on what reference budget,
and do the synthetic failure modes transfer to real data?** We contribute:
(i) the benchmark — six core suites, 508,825 fitted datasets plus 472,500
held-out evaluations, released generators with seed rules making any cell
independently reproducible; (ii) the protocol — matched-basis comparison
with empirical-FPR calibration and failure registration; (iii) the
failure-mode map — including a silent numerical failure of pinv-based
HC3-Wald inference on near-collinear physics-derived designs
(independently verified; the correction itself is not implicated), the
near-blindness of fitted-F RESET to rational laws, and the measured
dominance and calibration cost of the nested-rival gates; (iv) the
family-viability battery — 94,000 datasets showing rejection extends to
routing once — and only once — the null is rebuilt per coordinates; and
(v) a real-data suite showing the map's mechanisms operating outside the
laboratory; and (vi) the joint-release analysis — Section 4.10,
which connects the fit screen, the four gates, and extrapolation
risk on the same models: what the full pipeline releases, how
dangerous the releases are, and how the release channel scales
with sample size. What we
deliberately do not contribute is a winning detector;
the benchmark's own results (Section 4.7) demote its best test from
default gate to portfolio member, and the paper is organized around that
honesty.

## 2 Related Work

**Specification testing.** RESET [1] augments the regression with powers
of the fitted values and tests the joint significance of the added
coefficients; the exog variant uses powers of the regressors instead.
Heteroskedasticity-consistent covariance estimation [2] and its HC3
leverage-weighted refinement [3, 4] are the standard small-sample
corrections and are the default in several mainstream econometrics and
statistics packages. These results are derived for linear models under
fixed or quasi-fixed designs; none of the classical treatments
contemplates log-transformed physical laws whose regressors span only
a few orders of magnitude — the regime where our harness's own
pinv-based inference path silently failed (Section 5, Q3) and where
the verified SVD implementation of the same tests behaves nominally
(Section 4.1).

**Model selection as a specification gate.** Comparing a null model
against a nested rival by information criterion [12] or cross-validation
[13] is the oldest workflow-level misspecification check in applied
practice. It is so obvious that a reviewer of any new test will ask for
it; it is so workflow-dependent (which rival? which criterion? which
threshold?) that benchmarking it requires exactly the protocol this paper
introduces — a frozen rival family, frozen decision statistics, and
thresholds calibrated on a reference pool. Suite S5 is, to our knowledge,
the first head-to-head we know of between a resampling-null specification test and
nested-rival selection gates under matched false-positive rates.

**Resampling nulls under non-i.i.d. errors.** The wild bootstrap [10, 11]
was developed for inference under heteroskedasticity, and its permutation
cousin requires exchangeability. The battery's mechanism finding
(Section 4.6) is a compact new instance of this classical dichotomy: a
log transform homoskedasticizes multiplicative noise, so residual
permutation is valid in coordinates ending in log y and invalid in
coordinates ending in y or 1/y — the null's validity inherits the noise
geometry of the coordinates the test runs in.

**Symbolic regression benchmarks.** SRBench [5] established the
discovery-side evaluation: 14 SR methods and 7 ML baselines on 252
problems, scored by exact recovery and held-out accuracy. SRSD [6]
reconceived the task distribution itself — realistic per-variable ranges
reviewed against the physics — and is the source of our task
configurations. Both measure the discovery side; a pipeline that returns
a wrong expression and one that flags it score identically. Our
generators deliberately reuse SRSD's official sampling objects so that
the falsification-side pools are the same distributions the discovery
side already standardizes on.

**Physics-informed and dimensionally-aware discovery.** AI Feynman [7]
uses physical units to constrain the hypothesis space;
measurement-driven discovery distills governing equations and structure
directly from data [8]; PySR [9] represents the modern search side.
Dimensional correctness is necessary but not sufficient for form
correctness: every ground-truth law in our suites is dimensionally
self-consistent, yet each non-monomial law is a form misspecification for
a power-law model. The falsification side complements these constraints —
units narrow what to search; specification gates check what was found.

**Honest-failure evaluation.** Existing SR benchmarks report no analogue
of a false-positive rate because they contain no explicit rejection
mechanism to calibrate. We are not aware of a symbolic regression
benchmark that measures attached falsification behavior, per family, with
preregistered labels and failure registration — or one that benchmarks
its own proposed test against the nested-rival alternative its reviewers
would name. This is the gap the present benchmark occupies.

**The wider specification-testing literature.** Econometrics and
statistics hold many more form-specification tools than the gates tested
here: smooth residual tests, generalized likelihood-ratio and conditional
moment tests, rainbow and related one-degree-of-fitting checks,
spline/GAM residual comparison, kernel-based lack-of-fit statistics,
cross-fitted and sample-splitting variants, and conformal-style
selective-risk diagnostics. S5's nested-quadratic rival is the
discovery community's canonical objection, not a survey of this
literature; LOOCV, BIC, held-out likelihood, or non-nested rivals may
rank differently. The gates were chosen because they are what discovery
pipelines actually attach today (RESET variants, robust inference) plus
the workflow alternative a reviewer would name (nested-rival selection);
the harness is gate-agnostic — a new gate is one statistic-plus-decision
function over the same fitted cells — and extending the benchmark to
these tests is future work, not a closed set.

## 3 Benchmark Design

### 3.1 Design principles

**P1 — Common rejection basis.** Within any comparison, all gates decide
on the same basis at the same significance level α = 0.05. The
permutation curvature test's empirical p-value is computed from B = 999
permutations of the residual-permutation null; RESET variants use
asymptotic F reference distributions; S5's selection gates use thresholded
statistics, with thresholds at naive zero and at an H0-pool-calibrated
value — and naive versus calibrated results are reported side by side,
because raw tables are what current practice would see.

**P2 — Power only at matched empirical FPR.** Headline power comparisons
calibrate each method's decision threshold to a 5% empirical
false-positive rate, estimated on pooled true-monomial tasks within the
same suite, sample-size stratum, and noise regime.

**P3 — Preregistered labels.** H0 (strict monomial: the log-log mean
function is exactly linear) versus H1 (any non-monomial structure) is
fixed from symbolic structure of the ground truth before any test runs.
The family stratification of S3, the battery's truth registry, the S5
gate definitions, and the S6 relation list are likewise frozen in code
before their runs; no label or configuration is revised after seeing
results.

**P4 — Explicit failure registration.** Every evaluation that cannot
produce a valid decision is registered with a typed reason — zero-residual
degeneracy, design rank deficiency, insufficient residual degrees of
freedom, leverage degeneracy, singular covariance, statistic mismatch,
generation failure (S5 SRSD side), indicator archival (S6) — and excluded
from rate numerators while remaining in the denominators' audit trail.
The zero-noise grid cells exist to exercise this path.

**P5 — Official data configurations.** S2/S3 tasks reuse the official
per-variable sampling objects of the source registry; tasks whose inputs
or outputs are not strictly positive under the official configuration are
marked not-applicable and excluded — we do not shift or rescale data to
make them log-fittable, since that changes the hypothesis under test. S6
uses a public observational source (World Bank Open Data API [14]) with
the relation list frozen before any gate ran; the United Nations
Foundation itself publishes no dataset API, and the World Bank is the
UN-system institution with the open macro data of record.

### 3.2 Gates under test

All gates operate in log-log space on positive data: let ly = log y and
X = [1, lx₁, …, lxₖ] be the design of the candidate power-law model.

**A. Permutation curvature (form_check).** Adds a single mean-centered
quadratic column q = meanⱼ((lxⱼ − lx̄ⱼ)²) and takes the added coefficient
c₂ as statistic; size evidence comes from B = 999 residual permutations
(vectorized through the pseudo-inverse row), with a statistic-mismatch
sentinel comparing wrapper and native values on every call. Validity
note: permuting *estimated* residuals assumes they are exchangeable,
which OLS residuals are not exactly — the permutation null here is an
approximation whose quality is an empirical question, not a finite-sample
exactness theorem. What the benchmark establishes is calibration of this
construction on our H0 pools (Tables 1–2); log-transformed multiplicative
noise and favorable designs plausibly help, and §6 lists the stress
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

Configurations were frozen before all main runs — powers (2, 3); B = 999;
α = 0.05; the rival family and criteria of E/F; the grids below — and are
identical across suites. The same harness invokes all gates on every
dataset, so comparisons share data, seeds, and failure handling exactly.

### 3.3 Task suites S1–S3

**S1 — Synthetic control (19 scenarios × 500 seeds).** Three monomial H0
scenarios (univariate; two-variable; two-variable with ρ ≈ 0.9 correlated
design). Eleven structural H1 scenarios: quadratic curvature (c = 0.10);
cubic (standardized odd term); multiplicative interaction (log-space
cross term); threshold (exponent 0.4 → 1.0 at the domain midpoint);
exponential law; logarithmic law; rational (saturating) law; omitted
variable; dilution (one of two variables curved); opposing curvature
(+0.12, −0.12 on an exchangeable design, so the mean-square column's
population coefficient is exactly zero); collinear-quadratic (ρ ≈ 0.95
design with curvature). Five noise-layer H1 scenarios have strict-monomial
mean structure but non-Gaussian noise — Student-t(3); SD growing with the
log prediction (0.3× to 1.7×); 2% point outliers (+8σ); 6% group outliers
(+5σ); AR(1) ρ = 0.6 — rejections there are false alarms by construction
and are never pooled into structural power. Grid: n ∈ {50, 200, 1000};
relative Gaussian noise σ ∈ {0, 0.5%, 1%, 3%, 5%} in log space; 500
seeds per cell via a seed rule embedding the grid position.

**S2 — SRSD-Feynman tasks (32 tasks × 500 seeds).** From the
120-equation registry, 32 tasks pass the positivity screen under official
configurations: 13 strict-monomial H0 and 19 non-monomial H1. Each
training fit (n ∈ {50, 200, 1000}) is accompanied by a 2,000-sample
same-distribution test set and a far-field probe: the first regressor
extrapolated ×3 beyond both training domain edges, scored against the
closed-form truth in a numerically stable form.

**S3 — Family-stratified Feynman equations (30 equations × 30 seeds).**
Structural family classification from symbolic form: exact power law (12
H0); fractional power (1 H0); exponential (3 H1); rational (6 H1);
trigonometric (5 H1); plus not-applicable members with explicit reasons.
Each H1 equation runs over relative noise {0.5%, 1%, 3%, 5%} ×
{homoskedastic, heteroskedastic} × sampling-range scale {official, ×0.5,
×2}; H0 equations run the noise grid at official range.

### 3.4 Suite S4 — the family-viability battery

Each test of §3.2 evaluates a single null form; a rejection terminates
inquiry. The battery registers four families with their linearizing
coordinate pairs — power (log x, log y); exponential (x, log y);
logarithmic (log x, y); rational (1/x, 1/y) — and for each family F runs
A's construction in F's coordinates, returning the viability set
V = {F : p_F ≥ α} rather than a verdict. Columns are standardized before
squaring so features at different scales weight equally, and the null is
a wild bootstrap (Rademacher signs; B = 999) rather than a permutation —
for the reason §4.6 documents: **null validity is a coordinate-dependent
property**. The prediction (preregistered before any run) was that
permutation, valid wherever the v-coordinate is log y, would break in
coordinates ending in y or 1/y; the probe measured exactly that, and the
formal run uses the wild null throughout. Three errors are reported
separately and never merged: false exclusion (truth absent from V),
false inclusion (wrong family present), and escalation on structures no
family owns — the correct outcome, never scored as a win. Grid: n = 200;
rel ∈ {0.5–5%}; 500 seeds; 188 cells; 94,000 datasets, all 19 S1
scenarios with truth-family labels frozen in code plus 14 preregistered
SRSD H1 equations × {homoskedastic, heteroskedastic}.

### 3.5 Suite S5 — the gate comparison

S5 exists because of an objection this benchmark could not answer from
its own data: *why not skip the test statistic entirely and let a nested
rival arbitrate?* Gates E and F (§3.2) are that objection, made runnable.
Design, frozen in `gate_compare.py`'s docstring before the run: benches =
S1's 19 scenarios and S2's 32 tasks; grid n ∈ {50, 200}, rel = 3%, 200
seeds per cell (seed bases 300,000 / 400,000 — disjoint from all prior
suites); rival = full quadratic expansion; E and F thresholds naive and
H0-pool-calibrated per bench × n; A run under its standard protocol on
identical data. Metrics: H0 size (naive and calibrated), structural-H1
detection at calibrated FPR 5%, noise-layer fire rates (should stay
quiet), per-scenario and per-equation detail (never macro-averaged away),
and cost in fit units. One asymmetry is inherent and is reported rather
than hidden: A's permutation null is self-calibrating, while E/F at naive
thresholds inherit nested-comparison tie bias; the calibrated variants
price what a reference pool buys.

### 3.6 Suite S6 — real data (World Bank open macro data)

Synthetic and equation-derived suites have ground truth by construction;
the deployment question is whether the synthetic failure-mode map predicts
anything on observational data, where no truth exists. S6 is deliberately
small, preregistered, and divergence-oriented: six macroeconomic relations
(target, features) × four years (2019, 2021, 2022, 2023), cross-country
cross-sections of 184–199 sovereign economies (World Bank aggregates
excluded), from eight indicators pulled from the World Bank API with raw
responses cached for reproducibility — life expectancy ~ GDP per capita
(PPP); infant mortality ~ GDP per capita; GDP ~ population; CO₂ per
capita ~ GDP per capita; internet penetration ~ GDP per capita; and life
expectancy ~ GDP per capita + health expenditure per capita (two-feature).
Gates A–D run at 5 permutation seeds per cell. These seeds are
computational replicates of one p-value estimator on one fixed dataset —
not five independent samples — so they are used only to report Monte
Carlo stability (the empirical p has resolution 1/1000 at B = 999; across
the 24 S6 cells the median p spread is .009 and the max .094), never as a
sample-size or stability claim; a preregistered probe adds three agricultural
and energy relations over three years, including a built-in
control-contrast pair — cereal yield ~ fertilizer (one feature) versus
cereal yield ~ fertilizer + energy use per capita (two features) — chosen
to test, on live data, whether the dilution mechanism of S1 (a second
correlated curved feature halving the mean-square column's signal)
operates outside the laboratory. No calibration or power claim is made on
S6; its contract is external validity of *mechanisms*, not of rates. An
engineering companion bench (Appendix E) extends the same contract to
physical-engineering data: nine open sources, 91 datasets, 150
preregistered relations, run as a single-realization smoke pass.

### 3.7 Suite S7 — the joint-silence rate

Section 4.10 asks what the *pipeline* releases (fit screen, then
gates). S7 asks a different question: with the fit screen removed,
**how often does the entire gate panel stay silent on a wrong form?**
The two must never be pooled — S7's denominator is every evaluated H1
fit, not the screened subset — and they price different deployment
questions. Section 4.10 prices the pipeline a practitioner would
actually run; S7 prices the gates themselves.

Design. The 18 non-monomial H1 equations of S2 that survive the
far-field probe (`feynman-i.30.5` withdrew, v5 App. B). Per task and
n ∈ {50, 200, 1000}, 2,000 datasets — 36,000 per n, on a seed grid
disjoint from S2's 500-seed draw. Gates A/B at native thresholds,
C/D on the SVD path. Two threshold regimes: naive α = .05, and
H0-pool calibrated to the same per-gate 5% budget. Joint silence
requires that **no gate rejects and no gate declines**: a `fail_*`
registration is an honest refusal of the evaluation, not a
certificate, and counting it as silence would flatter the panel.

A second arm adds the model-comparison gates (CV-quad,
AICc-quad — the S5 members) on `gate_cmp`'s SRSD arm (3,800
datasets; a subset of the 20,400-cell S5 bench, not the whole of it)
to test whether selection-side gates seal the leaks the
specification gates leave open. Its H0 reference is that same
panel's silence rate on true power laws. Because the arms rest on
different benches, **their levels are not comparable across arms**;
only within-arm contrasts are used below.

Conditional cost is measured separately: among jointly-silent wrong
forms, the far-field (one-decade) extrapolation error distribution.
A silent certificate whose median far-field error is small is not
thereby safe; the claim is about the tail.

## 4 Results

Table numbering is consecutive: Tables 1–2 size, Tables 3–4 power,
Table 5 the family map, Table 6 noise, Table 7 extrapolation, Tables 8–10
the battery, Tables 11–12 the gate comparison, Tables 13–14 real data,
Table 15 the no-entailment summary, Table 16 the joint-silence rate. All
rates are rejection fractions; Wilson 95% intervals for headline rates
appear in Appendix C. Seed counts per cell: 500 (S1, S2, S4), 30 (S3),
200 (S5), 5 (S6 — permutation-p stability, not rate estimation), 2,000
per task (S7; equations, not datasets, are the independent unit there —
see §6).

### 4.1 Size: false-positive rates on true monomials

**Table 1 — S1 synthetic H0, false-positive rate at α = 0.05 (Gaussian noise pooled over σ ∈ {0.5%…5%}; 4,500 fits per cell).**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A — permutation curvature | .047 [.041–.053] | .037 [.032–.043] | .041 [.036–.047] |
| B — RESET fitted F | .045 [.039–.051] | .041 [.036–.047] | .038 [.033–.044] |
| C — RESET fitted HC3 | .048 [.042–.054] | .045 [.039–.051] | .038 [.033–.044] |
| D — RESET exog HC3 | .054 [.048–.061] | .042 [.036–.048] | .040 [.034–.046] |

**Table 2 — S2 physics-derived monomial tasks (13 SRSD-Feynman equations; synthetic draws from real equations, not measurements), false-positive rate at α = 0.05 (26,000 fits per cell; C/D on the verified SVD path).**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .049 [.046–.052] | .053 [.050–.055] | .052 [.049–.054] |
| B | .050 [.047–.053] | .050 [.047–.052] | .053 [.050–.055] |
| C (SVD) | .057 [.054–.060] | .053 [.051–.056] | .055 [.052–.058] |
| D (SVD) | .060 [.058–.063] | .056 [.054–.059] | .058 [.056–.061] |

*(For reference, the retracted pinv(X′X) path on the same data: C
.417–.422, D .563–.565 at every n — the numerical-implementation
failure case of Section 5 Q3, kept out of method comparison.)*

On synthetic Gaussian data, all four tests are calibrated (Table 1).
On the physics-derived designs the picture splits (Table 2, C/D on
the verified SVD path): all four remain within Wilson contact of
nominal — A .049–.053, B .050–.053, C .053–.057, D .056–.060, flat
in n. The 42–56% breakage this benchmark originally reported for
C/D on these tasks was a silent numerical failure of our own
harness's pinv(X′X)-based inference path on near-collinear built
designs — not a property of the HC3 correction. That failure, its
independent verification, and its retraction are reported as a
numerical-implementation case in Section 5 Q3; the method
comparisons of this paper use the SVD path throughout.

The pooled rates also conceal granularity, which matters
for any claim of generality over "physical laws": under the
original pinv path the failure was equation-specific, not uniform
(5 of 13 tasks at .998–1.00 for C while the other 8 sat at
.048–.060); under the SVD path no equation family shows systematic
deviation. The
Wilson intervals of Table 2 treat all 26,000 fits as independent draws;
because seeds repeat within a fixed equation, the honest uncertainty
unit for generalization over equations is the task. Task-level bootstrap
(resampling the 13 equations; Appendix C) gives the SVD-path C a 95%
interval of [.045, .068], against [.197, .707] under the retracted pinv
path; leave-one-equation-out moves the SVD pooled rate over
[.048, .062], against [.053, .054] for A. One qualification follows.
What makes the pinv-failing subset special is design conditioning
(Section 5, Q3). An independent-implementation check run in revision
(Section 5, Q3) reattributes the failure: a standard SVD-based
HC3-Wald implementation (statsmodels) is nominal-compatible on the
very tasks where our pinv-based path rejects at 1.00 (8/120
rejections = .067, Wilson [.034,.126] over 60 seeds × 2 tasks —
sufficient to reject the old rate of 1.00, not precision evidence of
exact 5% calibration), and the two implementations agree to ≥3
significant digits on well-conditioned tasks.

### 4.2 Power: raw, then calibrated

Raw power at α = 0.05 on S2 H1 (n = 200, pooled over rel ≥ 0.5%):
A .52, B .28, C .25, D .68 (C/D SVD). Calibrating every method to
a common 5% empirical FPR:

**Table 3 — S2 power at matched empirical FPR 5% (C/D on the verified SVD path; thresholds re-derived from SVD p-values on the H0 pool).**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .386 | .486 | .570 |
| B | .187 | .254 | .389 |
| C (SVD) | .129 | .240 | .422 |
| D (SVD) | **.494** | **.675** | **.733** |

The SVD correction reorders the field. Under the retracted pinv
path, C/D "power" was computed against thresholds derived from
artifact-contaminated p-value scales (mass piled at 0 and ~1e-195;
Section 5 Q3), making those rows meaningless — including the
previous draft's claim that A leads at small and medium n. With
verified C/D, **D (exog-HC3) has the highest calibrated pooled
power at every sample size**, and the cross-equation check
confirms it: task-level paired bootstrap (18-task pool,
i.30.5 excluded; 10,000 replicates) gives Δ(D−A) = −.080
[−.258, +.081] at n = 50 (crossing zero), **−.170 [−.329, −.041]
at n = 200 and −.148 [−.304, −.019] at n = 1000 — intervals that
exclude zero in D's favor**; D beats A on 44% of tasks at n = 50
and 72–78% at n = 200/1000. A remains the strongest gate that
requires no reference pool (its permutation null is empirically
size-valid everywhere measured, §4.7); B is dominated everywhere;
C adds nothing over D. The one paired difference outside the
four-gate family remains Δ(F−A) = +.177 [+.054,+.323] at n = 50
and +.176 [+.056,+.320] at n = 200 (F wins 63–68% of tasks; F
cannot be paired with D without integrating D into the S5
harness — outstanding). τ thresholds describe decision boundaries
and are no substitute for these intervals. No method exceeds .73 —
most of S2's H1 mass is hard for every test, which Section 4.3
decomposes by family.

**Table 4 — S1 blind-spot cartography (raw α, n = 200, rel = 1%).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| cubic (standardized odd term) | .40 | 1.00 | 1.00 | 1.00 |
| interaction (lx₁·lx₂ cross term) | .21 | 1.00 | 1.00 | .10 |
| opposing curvature (+c, −c, exchangeable) | **.02** | .98 | .98 | 1.00 |
| dilution (one curved variable of two) | 1.00 | .88 | .87 | 1.00 |

Three mechanisms, each verified by construction: (i) A's statistic is the
mean of squared centered columns — under opposing curvature on an
exchangeable design that mean has exactly zero population signal; (ii) a
quadratic-projection statistic loses power on odd-symmetric deviations;
(iii) B's fitted-value powers are a local Taylor surrogate for any smooth
alternative — saturating on interaction where A (.21) and D (.10)
struggle — and, as S3 shows, the same surrogate property near-blinds B on
rational laws.

### 4.3 Family-stratified results and sensitivity

**Table 5 — S3 calibrated power by family (C/D on the verified SVD path).**

| Family | A | B | C (SVD) | D (SVD) |
|---|---|---|---|---|
| exponential | .09 | .06 | .07 | .06 |
| rational | **.65** | .16 | .13 | **.67** |
| trigonometric | .66 | .69 | .58 | **.86** |

B — near-perfect on synthetic quadratics — is nearly blind to rational
laws (.16 calibrated): the fitted-value Taylor surrogate absorbs
precisely the low-order curvature by which a saturating law departs
from a power law. Under the SVD path D is the strongest gate on both
rational (.67) and trigonometric (.86) families; the previous draft's
near-zero C/D entries on trigonometric were the pinv artifact
(thresholds derived from artifact p-mass), and the S3-level
correction agrees in direction with the S2-level one of Table 3.
Exponentials defeat every gate (≤.09 — the linearization a log
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
rational D is .669 [.655, .683]; on the honest unit — the equation —
the same cells read .860 [.592, .998] and .669 [.409, .864]. Rational C
moves from .132 [.122, .142] to [.053, .277] the same way. The
equation-level intervals overlap across families and are wide enough to
contain orderings the point estimates contradict, so the family ranking
above is directional, not an ordered effect list; what survives the
wider intervals is the sign pattern — D ≥ C on rational and
trigonometric, every gate low on exponential — not the point ordering.

### 4.4 Noise-layer robustness

**Table 6 — S1 noise-layer H1 (strict-monomial mean structure) rejection rates at α = 0.05 — false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds; C/D on the verified SVD path).**

| Scenario | A | B | C (SVD) | D (SVD) |
|---|---|---|---|---|
| Student-t(3) noise | .034 [.021–.054] | .036 [.023–.056] | .030 [.018–.049] | .030 [.018–.049] |
| variance growing with prediction | .048 [.032–.070] | .058 [.041–.082] | .050 [.034–.073] | .050 [.034–.073] |
| AR(1) errors, ρ = 0.6 | .038 [.024–.059] | .040 [.026–.061] | .042 [.028–.063] | .042 [.028–.063] |
| 2% point outliers (+8σ) | .040 [.026–.061] | .066 [.047–.091] | .036 [.023–.056] | .036 [.023–.056] |
| 6% group outliers (+5σ) | .038 [.024–.059] | .060 [.042–.084] | .054 [.037–.077] | .054 [.037–.077] |

*Wilson 95% intervals per cell (500 seeds each). The C/D columns were the
last main-text cells still on the pre-correction path; re-computed from the
released S1 cells they reproduce the earlier values cell-for-cell
(`paper/review3_analysis.py` §3.8c), which closes that audit item.*

Every method, every non-Gaussian regime: .030–.066, with Wilson
intervals in the table (500 seeds per cell). The table's function
is exclusionary: the miscalibrations of Table 2 do not come from noisy
data; the design, not the noise, is what specification tests are fragile
to.

### 4.5 Do alarms predict extrapolation failure?

Base rate of far-field relative error above 10% across S2 fits: .198.

**Table 7 — Alarm precision and recall for far-field error > 10% (S2, all tasks, all n, rel ≥ 0.5%).**

| Method | precision | recall |
|---|---|---|
| A | .39 | .63 |
| B | .56 | .49 |

An A alarm raises the probability of severe extrapolation failure from
.198 to .39 (2.0×); a B alarm to .56 (2.8×). Of the true-positive
catches, 10,560 are A-only and 5,447 B-only against 12,575 shared — the
alarm sets overlap barely more than chance. But alarms are ordinal
evidence, not magnitude certificates: the suite's two worst extrapolators
(i.39.11, i.43.43, both 98% median far-field error) carry the least
in-domain curvature and are missed by every method, while ii.38.14 is
rejected by A essentially always yet extrapolates benignly (7.7%). Alarm
⟹ form evidence; silence ⟹ nothing. No rejection is evidence of
extrapolation safety.

### 4.6 The battery: from rejection to routing

**Table 8 — Probe: permutation vs wild null by family (n = 200, 200 seeds; entries rel = 1% / 3%).**

| Truth | perm null | wild null |
|---|---|---|
| power (size of power member) | .070 / .075 | .075 / .060 |
| rational: false escalation (V = ∅) | **.34 / .31** | .045 / .090 |
| rational: containment | .66 / .70 | **.96 / .91** |
| exponential: containment | .94 / .95 | .92 / .95 |
| logarithmic: containment | .87 / .93 | .89 / .95 |
| threshold (no battery family): escalation — *correct outcome* | 1.00 / 1.00 | 1.00 / 1.00 |

With the permutation null, the battery's own rational member rejects its
own true family a third of the time — nothing about rational laws is at
fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The
wild bootstrap repairs exactly that and nothing else.

**Set-level errors and multiplicity.** The battery's decision is the
viable *set*, so the errors that matter are set-level, and each member
runs at α = .05 without correction — a union decision whose two error
types move oppositely under multiplicity. Two pools, two behaviors.
On H0-side cells (truth is a battery member) the reported containment
is one minus the exclusion rate of the true family's own member
(.0496–.066 across the four truth families; the viable set is never
empty *on this pool*), and per-member α does not inflate that rate
because exclusion depends on one member only. Inflation appears on
the other side: some wrong member survives alongside the truth in
23.3% of power-truth cells (near-neighbour overlap:
exponential/logarithmic members retain power data) and 92% of
logarithmic-truth cells (rational retains log data) — reported as
near-neighbour overlap rather than suppressed in the released
tables. On the structural-H1 pool (72,000 cells whose truth is no
battery member by design) the viable set is empty — V = ∅, correct
escalation — at the rates of Table 10's first column: universally for
quadratic/threshold/omitted-variable (1.00), rarely for the
interaction/opposing/collinear structures (.001–.05) whose mismatch
directions the member families do not span. Routing there means
*which wrong family is modal* (Table 10, second column); the battery
is a candidate filtering stage whose surviving set the discovery
layer must still rank, not a final selection. Downstream selection
error after the battery is out of scope here and unmeasured.

**Table 9 — Battery on synthetic truths (500 seeds/cell, wild null).**

| Truth family | Containment (truth ∈ V) | False escalation (V = ∅) |
|---|---|---|
| power — size, 8 scenario groups incl. all noise layers | — | .0496 [.046–.053] pooled; per-scenario .043–.054 |
| exponential | .95 | .05 |
| logarithmic | .94 | .005 |
| rational | .93 | .066 pooled (≤ .082 worst rel) |

**Table 10 — Battery on out-of-family structures (escalation to V = ∅ is the *correct* outcome; column 1 is the V = ∅ rate).**

| Scenario | Correct escalation | Modal false V |
|---|---|---|
| quadratic / threshold / omitted variable | 1.00 / 1.00 / 1.00 | — |
| dilution | .92 | rational |
| cubic (odd term) | .24 | power (.68) |
| interaction | .05 | power (.95) |
| opposing curvature | **.00** | power+rational (.84) |
| collinear-quadratic | .03 | rational (.61) |

Table 10 is Table 4 re-expressed as routing errors: the power member
inherits A's documented blind spots, so those structures stay in V under
the power label instead of escalating. On the SRSD side the routing is
reported descriptively (median |V| 1–2 for 11 of 14 equations; i.30.5 →
{exponential} in 78%; ii.38.14 escalates on 67%), because surface-form
family labels dissolve under examination — the battery correctly rejects
Planck-type "exponentials" that no linearizable family owns.

### 4.7 The gate comparison: why not just compare against a nested rival?

**Table 11 — S5 size on H0 pools (200 seeds/cell; naive and calibrated thresholds for E/F).**

| H0 pool | A | E naive | E cal. | F naive | F cal. |
|---|---|---|---|---|---|
| synth n=50 (600 fits) | .080 | .212 | .050 | .132 | .050 |
| synth n=200 (600) | .040 | .183 | .050 | .135 | .050 |
| SRSD n=50 (2,600) | .051 | .167 | .050 | .086 | .050 |
| SRSD n=200 (2,600) | .038 | .152 | .050 | .123 | .050 |

**Table 12 — S5 structural-H1 detection (A at its native threshold; E/F at H0-pool-calibrated thresholds; mean over tasks).**

| | synth n=50 | synth n=200 | SRSD n=50 | SRSD n=200 |
|---|---|---|---|---|
| A (permutation curvature) | .774 | .782 | .348 | .493 |
| E (CV-quad) | .929 | .934 | .462 | .628 |
| F (AICc-quad) | **.930** | **.940** | **.525** | **.669** |

*OR(A, E) at each gate's own 5% threshold — a different, larger joint
error budget (.065–.102 across the four benches), reported separately:
.937/.941/.547/.675. It is not comparable to the rows above and is
excluded from this table.*

*"Matched FPR" applies to three of the four columns: A's own measured
size is within Wilson contact of nominal at synth n = 200 (.040), SRSD
n = 50 (.051), and SRSD n = 200 (.038), but not at synth n = 50, where
A fires at .080 — roughly 1.6× the budget available to E and F. That
column is not budget-matched and is reported as measured: F's margin
there (.930 vs .774) is read net of A's head start, not as a matched
comparison. The same standard the OR row is held to is therefore
applied to A: only the three budget-matched columns support a
like-for-like ranking.*

Three findings. **(a) The obvious alternative is the better default
detector, and the margin is exactly A's blind zone.** The quadratic
expansion carries per-feature squares *and cross-terms*, so INTERACT
(A .24/.26) and OPPOSING (A .03/.02) are caught at 1.00 by both rival
gates, and the two SRSD extrapolation-blind equations of §4.5 (i.39.11,
i.43.43) are detected at .99–1.00. The direction set, not the null
construction, was the binding constraint. CUBIC curvature defeats all
three gates equally (~.3): odd curvature is near-orthogonal to every
even-column design — the honest general form of the blind-zone lesson.
**(b) The rival gates do not come calibrated.** Because the rival nests
the null, ties and noise break toward escalation: naive thresholds
false-alarm at .15–.21 (CV) and .09–.14 (AICc), and calibration to 5%
required an H0 reference pool — knowledge a production gate does not
have. A's residual-permutation null is an approximation (§3.2) but was
empirically size-valid on every pool we measured (.038–.080; synth
n = 50 sits at .080 and is reported as measured) without any reference
set — A's remaining structural advantage. **(c) Does the calibrated
threshold transfer?** Calibrating E and F on all but one H0 equation and
evaluating on the held-out equation (13-fold per bench × n) gives
held-out sizes of .053/.050 (E, n = 50/200) and .051/.052 (F), with
per-equation spread [.010, .095]: reference-pool calibration is not an
in-sample artifact — it transfers across equations, though any single
held-out equation can land within a ±.04 band of nominal. **(d) The OR
union's recall advantage is a joint-budget artifact.** OR(A, E) with
each gate at its own 5% runs a joint false-positive budget of .065–.102
across the four benches (Table 11) — its TPR edge over E (+.085/+.047
on SRSD) is purchased with that larger budget, not with information.
Recalibrating the *combination* to a 5% joint budget (A's size fixed,
E's threshold tightened on the H0 pool) is only defined where A's own
size permits it: on synth n = 50, A alone fires at .080 — above the
5% target — so no tightening of E can bring the union to 5%, and the
union *is* A there (TPR .774, reported as uncalibratable, not as a
5%-budget cell). Where the recalibration is defined (n = 200, both
benches) the union reaches .940/.664 against F's .940/.669 — parity
with F, never a gain.
reading: at any common error budget F is the frontier; the union is
defensible only when a workflow explicitly budgets for its larger naive
size. Noise-layer scenarios leave all three gates quiet (.04–.07): the
gates measure form, not noise assumptions, and the §4.4 separation
survives its most obvious stress test. Cost in fit units: A ≈
2 full fits + 1 vectorized pseudo-inverse (999 permutations); F = 2 fits;
E = 12 — all negligible against the candidate flood a discovery pipeline
already fits, which is the actual status quo's cost point: the workflow
runs *none* of them today.

### 4.8 Real data: does the map predict live behavior?

**Table 13 — S6 World Bank relations: verdicts (median permutation p over 5 seeds for A; deterministic p for B–D; α = .05).**

| Relation (target ~ features) | Years all four gates agree | Divergent cells (which gate, which way) |
|---|---|---|
| infant mortality ~ GDP pc | 4/4 reject | — |
| CO₂ pc ~ GDP pc | 4/4 reject | — |
| internet penetration ~ GDP pc | 4/4 reject (p < 10⁻¹⁶) | — |
| GDP ~ population | 4/4 **quiet** (A .71–.91; B .19–.33) | — |
| life expectancy ~ GDP pc | 0/4 | A alone rejects 2019 (.035) & 2022 (.008); B alone rejects 2022 (.0002); C/D never |
| life expectancy ~ GDP pc + health exp pc | 2/4 quiet | B rejects 2022 (.0002) while A quiet (.108) |

Where signal is strong the gates agree — three relations rejected
unanimously every year, one true scaling relation cleared unanimously —
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
cross-country heterogeneity, or spatial dependence — no ground truth
exists to separate these, so the claim is signature-consistency, not
mechanism confirmation.

**Table 14 — S6 dilution probe: the same target with and without a second correlated curved feature (3 years × 5 seeds; detection = median gate verdict).**

| Relation | A | B | D |
|---|---|---|---|
| cereal yield ~ fertilizer (1 feature) | rejects 3/3 years (.011–.031) | rejects 3/3 (~.003) | rejects 3/3 |
| cereal yield ~ fertilizer + energy pc (2 features) | **quiet 3/3 years** (.067–.149) | rejects 2/3 (.025–.092) | rejects 3/3 (≤.018) |

This matches the S1 DILUTION scenario's signature — adding one
correlated curved feature silences A on a relation it was rejecting,
while the RESET variants continue to alarm — a pattern *consistent
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
locally approximate laws — all recorded, none adjudicated.

### 4.9 Runtime

Per-dataset wall time (S1 means over 500-seed cells): A 12.0 ms (n = 50)
to 25.3 ms (n = 1000); B/C/D 1.0–1.8 ms; the full four-family battery at
n = 200 runs in 24 ms; S5 gates E/F 1.5–3.5 ms. All gates are effectively
free at the scale of a discovery pipeline that already fits thousands of
candidates; the bottleneck resource is not compute but the decision
discipline.

### 4.10 What the full pipeline actually releases

Sections 4.1–4.9 measured the gates in isolation. This section
connects the three stages a discovery pipeline actually applies to
the *same* fitted model — the fit screen, the specification gates,
and the extrapolation consequence — and asks the joint question:
**which wrong models does the whole pipeline release, and how
dangerous are they?** The analysis runs on the S1 quality bench
(57,000 fits, 200 seeds) and the S2 pool with recomputed quality
(240,000 fits, 500 seeds), C/D on the SVD path throughout (the
pinv-based C/D of Tables 2–3 is a numerical-implementation failure
case, Section 5 Q3, not a method result). Failed or declined gate
evaluations are counted separately and never treated as acceptance.
Every number is regenerated by released scripts (`svd_main_results.py`,
`no_causal_link.py`) whose anchor checks against previously verified
reports fail the build on drift.

**Definition (released dangerous model).** A fitted model is
*released* when it passes the fit screen (R² ≥ 0.99, the threshold a
pipeline would actually apply), no gate rejects (A/B native, C/D
SVD), and no gate declines (fail_* registrations are honest
rejections of the evaluation, not acceptance). It is *dangerous*
when its far-extrapolation error against the registered ground
truth exceeds 10%.

**Both thresholds are conventions, and the two stages spend their
error budgets on different scales.** The 10% danger bound carries no
downstream error budget or decision loss — it was fixed before the
runs and kept for comparability — so the headline pair (0.67% of
screened wrong forms, 6.9% of dangerous ones) is threshold-specific
and moves with it. The R² ≥ 0.99 screen is likewise uncalibrated
against the gates' 5% error budget: an absolute goodness screen and a
significance-level gate are not commensurable, and part of "the gates
do not fire at large n" is this asymmetry rather than gate behavior.
The sensitivity curve (`paper/review3_analysis.py` §3.3a) runs the
threshold over {2, 5, 10, 20, 50}%. At n = 50 the dangerous share of
screened wrong forms moves 39.6% → 11.0% → 2.5% across the 2 / 10 / 20%
thresholds, and the release rate among screened dangerous models reads
33.4% → 13.2% → 6.4%; at n = 1000 the 20% threshold leaves no dangerous
model at all (B = 0). A headline pair therefore moves by roughly a
factor of two in each direction over a threshold range a reader would
consider reasonable, and by much more at 2%. That is why every rate in
this section is printed with its threshold attached and none is offered
as a constant. A screen calibrated to the gates' own 5% error budget
remains outstanding.

**Table 15 — Joint release on the same models (S2: synthetic draws from the
registered Feynman equations, not measurements; released as
`svd_main_results.md` §C). Every rate names its denominator.**

| n | Screened wrong forms (A) | Dangerous among them (B) | Released dangerous (C) | C / A | C / B | C / D |
|---|---|---|---|---|---|---|
| 50 | 24,060 | 2,655 (11.0%) | 350 | 1.5% | 13.2% | 7.8% |
| 200 | 23,932 | 2,296 (9.6%) | 106 | 0.4% | 4.6% | 0.9% |
| 1000 | 24,000 | 2,046 (8.5%) | 25 | 0.1% | 1.2% | 0.2% |

*A = wrong forms passing the R² ≥ 0.99 screen; B = those of them whose
far-extrapolation error exceeds 10% (dangerous); C = dangerous models that
pass the screen and are released by all four gates. D = all dangerous H1 fits
at that n, screened or not (12,369 / 12,348 / 12,046). C/B is the release rate
among screened dangerous models, C/D among all dangerous models — the two
coincide at n = 200 and n = 1000 (106/106, 25/25) but diverge at n = 50, where
616 further dangerous fits that no gate rejected fail the screen and therefore
never reach the pipeline's output. Earlier drafts printed C/B in this table
and C/D in the running text under the same label, "release rate of
dangerous"; both are shown here.*

Three readings. **(a) The screen concentrates, the gates miss.**
Among wrong forms that pass R² ≥ 0.99, 9.7% are dangerous — the
stratum mixing of the previous subsection — and of those dangerous
survivors, 6.9% (481/6,997 pooled) walk through all four gates
untouched. The joint event is rare *unconditionally* (0.67% of
screened wrong forms) but six times more frequent *conditioned on
being dangerous*: the wrong forms most worth catching are not the
ones most likely to be caught. **(b) The release rate of dangerous forms falls with n — through the
screen, not through the gates.** With all dangerous H1 fits as
denominator the rate is 7.8% → 0.9% → 0.2% (966/12,369, 106/12,348,
25/12,046); with screened dangerous models as denominator it is
13.2% → 4.6% → 1.2%. Both are reported because they answer different
questions (Table 15). At n = 200 and n = 1000 every released
dangerous model has already passed the R² ≥ 0.99 screen (106/106,
25/25): at scale the binding constraint is the fit screen, not any
gate. **This gradient is not, by itself, evidence that sample size
closes the release channel through detection.** Three effects move
the same counts in the same direction — (i) gate power rises with n;
(ii) numerator and denominator both change with n; (iii) the best
L₂ approximation to a wrong form *within the power-law class* is
estimated better as n grows, so its extrapolation error contracts
toward a fixed projection error — an approximation-theoretic effect,
not a testing property. The dangerous share of screened wrong forms
falls 11.0% → 9.6% → 8.5% in the same direction and is consistent
with (iii) alone. Decomposing (i)–(iii) requires holding the danger definition fixed
across n. Two controls (`paper/review3_analysis.py` §3.2). **(1) Paired
label flips.** On the 36,000 (task, rel, seed) triples carrying all
three n, the 10% label flips almost symmetrically between n = 50 and
n = 200 — 966 triples dangerous only at n = 50 against 945 dangerous
only at n = 200 — so at that step the label chases noise around the
threshold rather than contracting; the wider steps are one-sided (994
vs 671 from n = 50 to n = 1000; 417 vs 115 from n = 200 to n = 1000).
**(2) Fixed danger label.** Holding the label at its n = 1000 value and
re-stratifying, the release rate among *screened* dangerous models
reads 16.8% → 1.9% → 1.2%, and among all dangerous fits under the same
fixed label 2.9% → 0.3% → 0.2%. **The gradient survives the control** —
it is if anything steeper at n = 50, because the fixed-label denominator
drops the fits that merely look dangerous at small n (and those are
released least often). So (iii) is real but is not what produces the
fall: the fall is the screen and the gates acting on a fixed set of
dangerous models, not the danger definition migrating. The defensible
reading remains the screen's: a practitioner
applying an R² ≥ 0.99 gate at n ≥ 200 is protected mostly by the
screen, and the gates contribute the residual at small n, where the
screen is weakest. **(c) The certificate
carries information — bounded information.** P(full-green
certificate | dangerous wrong form) = .013 against P(certificate |
true form) = .845: a dangerous wrong form is 65× less likely to
hold the certificate (likelihood ratio .015, exact reading in
`svd_main_results.md` §C.3), so the certificate is real evidence — but
P(dangerous | certificate) is not zero — 2.2% under the pool's own
mixing prior (19 H1 of 32 tasks), 1.5% under a 50/50 prior; the
pool-dependent reading is in `svd_main_results.md` §C.3 — and 481
dangerous models hold one. The correct semantics:
the certificate shifts the odds substantially and licenses no
exclusion.

**What the low-fit cases do and do not show.** Wrong forms below
R² = 0.5 are jointly missed at 57.5% (n = 50) — but a model at
R² = 0.002 would not survive any sane pipeline's fit screen; its
interest is mechanistic (the mismatch direction is orthogonal to
every augmentation direction the toolbox spans, ii.11.28), not
epidemiological. The joint-release accounting above is the
pipeline-relevant quantity; the low-fit silences belong to the
gates' blind-spot map (Section 4.4, Table 4).

**Direction geometry as diagnostic, one reading among others.**
Detection correlates with the alignment between the misspecification
direction and the test's augmentation space: corr(ŷ², residual)
against fitted-F detection at Spearman +0.47 (n = 50) to +0.76
(n = 1000), against R² at −0.06 to −0.15. We offer this as a
*diagnostic*, not an identified mechanism: the alignment statistic
is itself built from the same fitted values RESET augments with, so
its predictive success is partially circular — an independent test
would vary alignment while holding misfit magnitude fixed, and the
full augmentation space (ŷ³, exog powers, cross-terms), not ŷ²
alone, defines the geometry. What the correlation establishes is
narrower: detection is predictable from quantities orthogonal to
fit magnitude, and no fit statistic dominates that prediction.

**The synthetic bench is the control.** There, wrong forms passing
R² ≥ 0.99 are jointly silent on 0–0.1% of samples (1/811 at n = 50,
0 beyond): synthetic deviations were designed to be detectable,
which is why the real-formula release numbers are the finding and
the synthetic ones the foil.

### 4.11 The joint-silence rate: how often does the whole panel stay quiet?

Section 4.10 measured what the pipeline releases *after* the screen.
Strip the screen away and what remains is a different quantity: the
gate panel's own detection record. That record — not the screen's
selectivity — is what a practitioner inherits whenever a fit passes on
grounds other than an R² threshold: a hand-built correlation, a
literature formula being re-certified, a model carried in from an
earlier study. The two quantities are not ordered as larger or
smaller; they answer different questions, and §4.11 contrasts them
only within their own conditioning.

**Table 16 — Joint silence on S2's 18 non-monomial H1 equations (2,000 datasets per task per n; 36,000 per n; gates A–D, C/D on the SVD path). Joint silence = no gate rejects AND no gate declines.**

| n | naive α = .05 | H0-pool calibrated |
|---|---|---|
| 50 | 37.8% | 38.6% |
| 200 | 29.8% | 30.2% |
| 1000 | 23.2% | 23.6% |

Three readings.

**(a) The synthetic foil does not transfer.** On the synthetic bench
(S1, same four gates, **no screen applied**) joint silence runs
0.15–0.22% at n = 50. On equations drawn from real physics the same
four gates are jointly silent on **37.8%** of wrong forms at the same
n — a factor of roughly 170–250 on a like-for-like basis. That
like-for-like contrast is the one that carries a claim; comparison
against Section 4.10's screened rate would not, and is not made here.
The synthetic gallery was built with deviations designed to be
detectable; real formulas were not built at all. On this bench this
is the benchmark's largest genuine-versus-synthetic gap, and it is
offered as such — a gap between two benches, not an estimate of any
real-world error rate.

**(b) Adding selection gates does not close it.** At matched
calibration the four-gate panel is silent on **38.6%** of wrong forms
at n = 50 (Table 16, calibrated column). Adding the model-comparison
gates (CV-quad ∪ AICc-quad on `gate_cmp`'s SRSD arm) raises that to
**42.2%** rather than lowering it, because the added gates carry
their own 5% error budget and the union's silence on true power laws
is **90.5%**. The likelihood ratio a silent certificate can carry is
therefore only 42.2 / (42.2 + 90.5) ≈ **2.2×**; at a 50/50 class
prior the posterior that a silent fit is nevertheless wrong is ≈
**32%**. Silence is weak evidence, not safety.

**(c) The cost of silence is tail-loaded.** Among jointly-silent
wrong forms at n = 50, far-field error is *small at the median* (0.8%)
— which is exactly why the panel is silent — but **7.3%** of silent
certificates carry >10% far-field error (p90 = 8.7%). At n = 200 and
n = 1000 that share falls to 1.0% and 0.3%. The panel's silence is
not uniform in consequence; it concentrates where the misfit is small
in-domain and damaging one decade out.

**Relation to Section 4.10 — the denominators must not be pooled.**
These are not rival estimates of one quantity. Section 4.10's rates
are *conditioned on* the R² ≥ 0.99 screen and answer "what does the
pipeline emit"; S7's are *unconditioned* and answer "what does the
panel alone detect". The n = 50 pair — 13.2% of screened dangerous
models released (§4.10, C/B) versus 37.8% of all wrong forms jointly
silent (Table 16) — differs in large part because the screen admits
only about a tenth of wrong forms (24,060 of 240,000 S2 fits, §4.10)
before the gates ever run. Quoting either as "the release rate"
without its conditioning would repeat the denominator error Section
4.10 documents, one level up.

## 5 Analysis

**Q1. Does the low-df curvature test beat RESET on univariate quadratic
deviations?** No — it ties; both saturate from n = 50. What low degrees
of freedom buy is calibration, not power — which is what makes A's
calibrated numbers trustworthy rather than large. The asymmetries that do
appear are structural (odd terms, interactions), where A's quadratic
projection is mathematically the wrong instrument.

**Q2. Does averaging squared columns miss feature-specific curvature and
interactions?** Yes, by construction — and the benchmark quantifies when
that matters, in three escalating venues: synthetic power cells (Table 4),
routing errors (Table 10), and live observational data (Table 14). The
dilution mechanism in particular is now a three-level result: constructed
(S1), routed (S4), and signature-consistent on observational data (S6). A
synthetic blind spot stops being a laboratory artifact the moment it
anticipates which real relation goes quiet — though on S6 that
correspondence remains pattern-level, not mechanism-proof.

**Q3. Does HC3 fix false alarms under heteroskedasticity — and was the
observed failure HC3's or ours?** Heteroskedasticity was never the
disease (Table 6). The observed failure was **ours — a numerical
implementation failure, now fully retracted as a method result and
reclassified as a numerical-implementation case study**. The method
comparisons of this paper (Tables 2, 3, 5, 15) use the SVD path
throughout; the pinv(X′X) path survives only here, as the case
study. Table 2 (original draft) showed *our* HC3-Wald path breaking
on near-collinear built designs; the original draft attributed this
to the correction itself ("the robust estimator is the failure
mode").
independent-software checks required in revision retract that
attribution, and the second round extends them from fitted-HC3 (2 tasks
× 60 seeds) to **exog-HC3 across all 13 H0 equations with stratified
sampling by equation, sample size, and noise level (176 fits)**, on
identical data frames, comparing three paths: our released pinv(X′X)
path, our own SVD path (hat and sandwich via pinv(X) of the augmented
design), and statsmodels HC3-Wald with the same F reference. Results:
(a) the SVD path and statsmodels agree everywhere — maximum |Δp| =
.0023, maximum relative statistic difference .0038, covariance-matrix
Frobenius relative difference ≤ 3×10⁻⁴ — so the mathematical definition
is implemented identically and the divergence is numerical; (b) the
pinv(X′X) path rejects at rate 1.00 on the 7 tasks whose measured
cond(X′X) ≥ 4×10¹⁷ (measured cond(X) 6.5×10⁸–3.2×10¹⁰), while both SVD
paths stay at 0–.10 on every task including those — compatible with
nominal, not precision-calibrated; (c) the two quiet tasks with cond(X)
≈ 10 (ii.27.16, ii.27.18) show no failure on any path. Mechanism,
stated precisely: the harness calls NumPy 2.5.2 `np.linalg.pinv`,
whose default cutoff zeroes singular values below rcond = 10⁻¹⁵ ×
σ_max(the max(m,n)·ε default belongs to SciPy's pinv, which we do not
use); forming X′X squares the condition number, and per NumPy's own
conditioning documentation a condition number 10ᵏ loses ≈ k digits in
the worst case — at measured cond(X′X) = 10¹⁷–10²⁰ the smallest
singular values fall below the cutoff and are silently discarded, so
the hat diagonal and sandwich are computed from a destroyed spectrum.
The earlier draft's informal precision arithmetic is superseded by
these measured values. There is no per-call rolling cross-check in the
harness — the independent comparison is the released one-off script
described above. Corrected conclusion: a silent numerical failure of
pinv(X′X)-based inference on squared-condition-number designs; the
practical hazard is real (research code commonly pinv-solves; the false
alarms are total and silent), but it is an implementation robustness
lesson, not evidence against HC3. Accordingly, all method comparisons
in this revision (Tables 2, 3, 5, 15) run C/D on the SVD path; the
pinv numbers above are the case study, and any residual difference
between the paths on suites not yet re-run is an open item (§6).

**Q4. Do alarms predict extrapolation error?** As binary evidence,
usefully (2.0–2.8× enrichment, complementary alarm sets); as magnitudes
or guarantees, no. The correct reading is asymmetric: alarm ⟹ form
evidence; silence ⟹ nothing.

**Q5. Can the same machinery answer "not a power law — then what?"** Yes,
as a viability set, after one repair the main benchmark could not have
predicted: nulls are coordinate-dependent, and the wild bootstrap is the
repair. The formal run delivers the contract (containment .93–.95,
pooled size .0496, correct escalation 1.00 on quadratic/threshold/
omitted-variable), with near-neighbour overlap reported rather than
suppressed and A's blind spots wearing routing clothes (Table 10).

**Q6. Why not just compare against a nested rival?** Measured, it wins.
F (AICc-quad) leads every cell of Table 12 and removes A's blind spots,
because the rival's direction set (per-feature squares and cross-terms)
covers the directions A's single mean-square column cannot project onto.
The purchase price is a threshold: naive E/F false-alarm at .09–.21 until
an H0 reference pool calibrates them; the calibration transfers across
equations (leave-one-equation-out held-out sizes .050–.054 on average,
§4.7(c)), while A was empirically size-valid on every pool we measured
without one (§3.2 records the approximation this rests on). The OR
combination, originally reported as the best-recall S5 cell, owes that
edge to a larger joint false-positive budget (.065–.102); recalibrated
to a 5% *joint* budget it matches but never beats F (§4.7(d)). The
portfolio reading follows: a pipeline that can afford a reference pool
should default to F; one that must report a p-value per candidate with
no reference data should default to A; and a recall-critical deployment
should budget for F rather than pay for an uncalibrated union. The
benchmark's best test is thereby demoted from default gate to
self-calibrating portfolio member, by the benchmark's own protocol; we
consider this the paper's most credibility-relevant result, because it is
the one a competing paper would have had an incentive to omit.

**Q7. Does any of this transfer to real data?** Where the synthetic map
predicts strong signal, four gates built for physical laws agree on
macroeconomic relations they were never designed for (Table 13) — cheap
evidence that the map is not a laboratory artifact. Where the map
predicts weakness, observational divergence appears with the predicted
signature (dilution-consistency, Table 14) — consistent-with, not
confirmation-of: too many observational confounds and no truth to
adjudicate. What S6 and Appendix E cannot supply — rates, calibration,
extrapolation truth, mechanism identification — they do not attempt;
real data here is a compatibility probe, not a leaderboard, and their
rejection counts are compatibility statistics, not performance
rankings.

**Q9. Does fit quality certify functional form — or does the full
pipeline catch what the screen misses?** Neither certificate nor
guarantee, and the honest quantity is the joint one (Section 4.10).
Fit quality cannot by itself certify form: within R² ≥ 0.99 the two
benches mix true forms (extrapolation error ≈0.1%) and wrong forms
(100% >10% error on the synthetic bench — a design property of that
bench; 9.7% of screened wrong forms are dangerous on the
Feynman-equation pool) with the
certificate unable to separate them. The full pipeline releases
dangerous wrong forms at 0.67% of screened wrong forms but 6.9% of
dangerous ones — the joint miss concentrates exactly on the models
most worth catching — and the release rate of dangerous forms falls
with sample size (7.8% → 0.2% with all dangerous fits as denominator,
13.2% → 1.2% among screened dangerous models); Section 4.10(b) shows
that fall running through the fit screen rather than through any
gate. The
certificate nonetheless carries real, bounded information: a
dangerous wrong form is 65× less likely to hold a full-green
certificate than a true one. What no configuration of these tools
does is *exclude* the failure mode; certificates of correctness
must be anchored in which alternative directions were tested, not
in fit magnitude.

## 6 Limitations

**Single candidate family, and no symbolic-regression method (main
benchmark).** All gates of §3.2–3.5
evaluate a power law; the battery addresses routing, not selection —
ranking within V is the discovery layer's job and is not benchmarked.
No symbolic-regression method is run anywhere in this benchmark: every
gate is applied to a pre-specified candidate, so what is measured is
gate behavior on a power-law null, not the pass rate of
*method-emitted* wrong forms. The SRBench/SRSD comparison of §1 is
therefore an argument for the missing side, not a claim of parity with
discovery-side benchmarks, and a suite that attaches these gates to at
least one SR method's output is the next step — named here so that no
reader mistakes this paper for having taken it.
**One rival, one criterion each (S5).** The rival is the quadratic
expansion and the criteria are 5-fold CV and AICc — the reviewer's
canonical objection, not a survey of selection gates; LOOCV, held-out
likelihood, or non-nested rivals may rank differently, and the
calibration asymmetry (reference pool vs self-calibrating null) is a
workflow property whose cost we price but do not optimize. **Battery
scope.** n = 200 only; no member for interaction or trigonometric
structures; SRSD side descriptive pending a symbolically exact
linearizability registry. **Real-data scope (S6).** Six relations, one
source, four years, no ground truth: no rate claims are possible, the
relation list is preregistered to prevent cherry-picking but remains six
choices, and the year-to-year verdict flips in weak-curvature cells
signal genuinely boundary-level signal rather than gate defect — these
cells are reported as measured, not smoothed. **Fixed augmentation.**
Built powers fixed at (2, 3). **Log-space scope.** All gates operate in
log-log space; additive-error models are out of scope and marked
not-applicable. **Statistical precision and clustering.** S1/S2/S4 rates rest on 500
seeds; S3 on 30 (directional); S5 on 200; S6's 5 seeds are Monte Carlo
replicates of one estimator, not samples (median p spread .009, max
.094). Seeds repeat within equations, so Table 2's Wilson intervals
understate cross-equation uncertainty; task-level bootstrap intervals
and leave-one-equation-out ranges (Appendix C) are the honest units for
generalization over equations, and the HC3 failure is concentrated in 5
of 13 equations. No multiple-comparison correction is applied across the
1,480 cells; rates are descriptive, and family-level set decisions are
analyzed as sets (§4.6) rather than as four independent tests.

**Joint-silence suite (S7).** The suite's two arms rest on different
benches — A–D on the S2 18-equation pool (36,000 datasets per n) and
the model-comparison extension on `gate_cmp` (3,800) — so their
levels are reported separately and are **not** comparable across
arms; only within-arm contrasts are used. The four-gate arm's H0
reference was not re-measured on the S7 seed grid, so Table 16 is an
H1-side measurement: the 90.5% true-form silence rate quoted in
§4.11(b) belongs to the **model-comparison arm alone** and is not
the four-gate arm's null. Table 16 carries no Wilson intervals — the
2,000 datasets per task are draws within 18 equations, so equations,
not datasets, are the independent unit, and the per-task spread is
wide by design (per-task detail is retained rather than
macro-averaged away). Two exclusions predate the suite:
`feynman-i.30.5` (withdrawn, v5 App. B) and `rel = 0` cells (no
certificates exist there). **Provenance.** S7's cached outputs are
`runs/joint_silence/`. The summary table in that directory's
`REPORT.md` disagrees with the per-task CSV it was computed from
(32.0% vs 37.8% at n = 50 after weighting; 70% vs 90.5% for the H0
baseline) and carries a dated in-file annotation to that effect.
Every number in Table 16 is regenerated from the per-task CSV
(`joint_silence_srsd_per_task.csv`), and the discrepant summary is
not used anywhere in this paper.
**HC3 independent verification.** After two revision rounds: fitted-HC3
(C) is verified on 2 failing tasks × 60 seeds; exog-HC3 (D) is verified
across all 13 H0 equations with stratified sampling (176 fits) — in
both cases our SVD path matches statsmodels to ≤ .0023 in p and ≤ 3×10⁻⁴
in covariance, while the X′X path fails only at measured
cond(X′X) ≥ 4×10¹⁷. This revision goes further than the previous one:
all main-text method comparisons (Tables 2, 3, 5, 15) run C/D on the
SVD path, and the previous draft's path-specific footnotes are replaced
by re-derived numbers. The S1 noise-layer rows of Table 6 — the last
main-text C/D cells still on the pre-correction path — are now re-run
on the SVD path and reproduce the reported values cell-for-cell (§4.4),
with Wilson intervals attached. Still outstanding: the S5 integration
of D, and regeneration of the Appendix D figure panels (below). **Permutation validity.** The
residual-permutation null is an approximation whose finite-sample
behavior under leverage, dependence, or stronger heteroskedasticity is
unmeasured here; leverage-adjusted, wild, block, cluster, and
fixed-design null variants are the required stress set. **OR budget.**
Union gates are reported at both naive and joint-5% calibration; any
deployment comparison must fix the budget first (§4.7(d)). **Threshold
conventions and the n-gradient (§4.10).** The 10% danger bound and the
R² ≥ 0.99 screen are conventions fixed before the runs, not calibrations
against a downstream error budget, and the two stages therefore spend
their error budgets on incommensurable scales. Both headline
joint-release quantities are threshold-specific. The release gradient
across n mixes at least three effects — rising gate power, a changing
numerator and denominator, and the contracting best L₂ approximation to
a wrong form within the power-law class — and this paper does not
separate them: the threshold-sensitivity curve and the n-decomposition
are registered as outstanding (`paper/review3_analysis.py`), and §4.10's
gradient claim is weakened to the screen's binding-constraint reading
for that reason. **Figure provenance (Appendix D).** The Appendix D panels are regenerated on the verified SVD path (`paper/review3_analysis.py`): D.1 plots the H0 p-value CDFs of both pools, D.2 the per-task calibrated power, D.3 calibrated power versus n with task-level intervals, and D.4 the alarm precision/recall and error distribution of §4.5 — D.4's base rate (0.198) and A/B precision-recall reproduce Table 7 exactly. The pre-correction figure set, whose captions had drifted from the published images, is retired; D.5 (timing) is unchanged because timing does not depend on the inference path.

## 7 Reproducibility statement

Every generator is seeded from a documented rule embedding the grid
position, so any *numeric* cell of any table regenerates
bit-identically, and the Appendix D panels are regenerated from the released cells by `paper/review3_analysis.py` (Section 6, "Figure provenance"). The
harness records content hashes of the code under test — including the
exact `form_check` implementation — and of the harness itself; a
statistic-mismatch sentinel verifies on every call that the wrapper's
recomputation agrees with the native implementation. The battery is
implemented at the experiment layer over the unmodified core code path,
so audited core-code hashes remain valid. S5's gate definitions and grid
are frozen in `gate_compare.py`'s preregistration docstring (seed bases
disjoint from all prior suites). S6 pulls from the World Bank API with
every raw response cached under the run directory (`runs/wb_real/raw/`),
making the observational suite reproducible against the cached snapshot
and refreshable against the live API; the relation and probe lists are
frozen in `wb_run.py` / `wb_blindspot.py` before any gate ran. The
revision's re-analyses are scripted and released with their raw inputs:
`paper/revision_analysis.py` regenerates every number in §4.1's
task-level paragraph, §4.6's set-level errors, §4.7's OR-budget and
transfer analyses, and §4.8's Monte-Carlo spreads from the released
per-cell parquet/CSV files; `paper/verify_round2.py` regenerates the
second-round D-verification, Wilson and paired-bootstrap numbers; the
independent HC3 check (Q3) uses statsmodels 0.14 on identical data
frames. The joint-release analysis of Section 4.10 and the SVD-path
Tables 2/3/5 numbers are regenerated by `svd_main_results.py` (from
runs/cd_svd and the released quality benches, anchor-checked against
the verified summaries). The no-entailment analysis of Section 4.10 (and its twelve
anchor checks against the previously verified reports) is regenerated
by `no_causal_link.py` from the released quality benches
(`runs/quality50/`, `runs/srsd_full/quality_recompute.parquet`) and
the SVD-path C/D rerun (`runs/cd_svd/`), writing `no_causal_link.md`
and `figures/no_causal_link.png`; the anchor list includes the
joint-silence rates of the companion audit, so silent drift between
the paper's two evidence chains fails the build. The engineering companion
bench of Appendix E is released the same way (`real_data/`): raw
snapshots cached with SHA-256 hashes and access dates (2026-09-06), a
manifest of preregistered relations with per-relation notes and
per-source license fields, and replayable fetch/prepare/smoke scripts;
every declined evaluation carries a typed reason. The benchmark artifact
(repository archive with all generators, per-cell results, and
content-hashed configurations) is deposited with the publication; a
permanent DOI is assigned at release. All cells' raw results (parquet and CSV), failure registrations
with typed reasons, selection and exclusion tables, analysis scripts, and
figure sources are released. Zero-noise cells deliberately exercise the
failure-registration path; the battery's 94,000 datasets and S5's 20,400
datasets contain zero data-layer failures, which is itself a registered
outcome.

## 8 Conclusion

We asked what the fit-quality metrics that certify models in
discovery pipelines can and cannot deliver, and measured both
directions on two benches (508,825 fitted datasets plus a quality
re-evaluation over 297,000 of these fits). Fit quality cannot by
itself certify a functional form: within the single R² ≥ 0.99
stratum, true forms and wrong forms share the certificate while
their extrapolation behavior differs by two orders of magnitude.
The evaluated specification gates share blind spots: on real
formulas, dangerous wrong forms (extrapolation error >10%) that
pass the R² ≥ 0.99 screen are jointly missed by all four gates at
6.9% — the joint miss concentrates on the models most worth
catching — and the release rate of dangerous forms falls with sample
size (7.8% → 0.2% with all dangerous fits as denominator, 13.2% →
1.2% among screened dangerous models), a fall Section 4.10 attributes
to the fit screen rather than to any gate. The
certificate carries real, bounded information (a dangerous wrong
form is 65× less likely to hold one than a true form); what it
cannot do is exclude the failure mode. These are bounded claims:
we quantify the extent, frequency, and sample-size scaling of the
blind spots — not the absence of a statistical relationship, which
our own likelihood-ratio accounting contradicts.

The falsification-side benchmark built to measure this is itself a
deliverable: seven suites, one protocol that puts every gate on a
common basis, and full registration of every failure — with all
C/D method comparisons on the independently verified SVD
implementation. Under that correction the gate ranking changes
substantively: D (exog-HC3, SVD) has the highest calibrated pooled
power at every sample size on the Feynman-equation pool, with the
cross-equation bootstrap confirming the advantage at n ≥ 200 — the
previous draft's small-sample lead for A was an artifact of
retracted numerics. The benchmark's most striking claim did not
survive its own revision: the 42–56% HC3 false alarms were a silent
numerical failure of a pinv-based inference path on near-collinear
built designs — an SVD-based implementation is calibrated exactly
where it fails — a retraction that is itself the benchmark working
as designed, and the reason this revision replaces rather than
footnotes the affected numbers. Raw power rankings still invert
across families: no gate covers all three families of Table 5 —
exponentials defeat every gate — and the nested quadratic rival adjudicated by AICc remains the
strongest default where a reference pool exists — paying for it
with a naive 9–14% false-alarm rate until a reference pool
calibrates it — while the permutation curvature gate remains the
strongest gate requiring no reference pool, and the OR union's
apparent recall premium is a joint-error budget artifact that
vanishes at a common budget. The battery extends the machinery
from rejection to routing, after a one-line null repair that
carries a general lesson: a null distribution's validity is a
coordinate-dependent property. And on macroeconomic data the map
does what a map is for: it predicts where four gates will agree,
and which gate goes quiet on which relation.

For the pipeline designer the results compose into one rule:
certificates must be anchored in which alternative directions were
tested and found silent — never in fit magnitude, and never in the
absence of alarms from a gate set whose blind directions have not
been mapped. In scientific discovery, the model you release on the
strength of its fit and its clean diagnostics is the one that
breaks your extrapolation.

## References

[1] Ramsey, J. B. (1969). Tests for specification errors in classical linear least-squares regression analysis. *Journal of the Royal Statistical Society, Series B*, 31(2), 539–551.

[2] White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4), 817–838.

[3] MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305–325.

[4] Long, J. S., & Ervin, L. H. (2000). Using heteroskedasticity consistent standard errors in the linear regression model. *The American Statistician*, 54(3), 217–224.

[5] La Cava, W., Orzechowski, P., Burlacu, B., de França, F. O., Virgolin, M., Jin, Y., Kommenda, M., & Moore, J. H. (2021). Contemporary symbolic regression methods and their relative performance. *NeurIPS Datasets and Benchmarks Track*.

[6] Mundhenk, T. N., Landajuela, M., Glatt, R., Santiago, C. P., faissol, D. M., et al. (2022). Rethinking symbolic regression datasets and benchmarks for scientific discovery. *arXiv:2206.10540*.

[7] Udrescu, S.-M., & Tegmark, M. (2020). AI Feynman: A physics-inspired method for symbolic regression. *Science Advances*, 6(16), eaay2631.

[8] Bakarji, J., Champion, K., Kutz, J. N., & Brunton, S. L. (2022). Discovering governing equations from partial measurements with deep delay autoencoders. *Nature Computational Science*.

[9] Cranmer, M. (2023). Interpretable machine learning for science with PySR and SymbolicRegression.jl. *SciPy Conference*.

[10] Liu, R. Y. (1988). Bootstrap procedures under some non-i.i.d. models. *Annals of Statistics*, 16(4), 1696–1708.

[11] Davidson, A. C., & Flachaire, E. (2008). The wild bootstrap, tamed at last. *Journal of Econometrics*, 146(1), 162–169.

[12] Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach* (2nd ed.). Springer.

[13] Arlot, S., & Celisse, A. (2010). A survey of cross-validation procedures for model selection. *Statistics Surveys*, 4, 40–79.

[14] World Bank (2026). World Bank Open Data API, indicators NY.GDP.MKTP.PP.CD, NY.GDP.PCAP.PP.CD, SP.POP.TOTL, SP.DYN.LE00.IN, SP.DYN.IMRT.IN, EN.GHG.CO2.MT.CE.AR5, EN.GHG.ALL.MT.CE.AR5, SH.XPD.CHEX.PC.CD, IT.NET.USER.ZS, AG.YLD.CREL.KG, AG.CON.FERT.ZS, EG.USE.PCAP.KG.OE. https://api.worldbank.org (cached snapshot released with the benchmark).

[15] UCI Machine Learning Repository. archive.ics.uci.edu (datasets: auto_mpg, automobile, concrete, slump, energy_efficiency, yacht, ccpp, air_quality, gas_turbine, tetouan, steel). Accessed 2026-09-06; raw snapshots and SHA-256 hashes released with the benchmark; redistribution per each dataset's stated license.
[16] NIST Statistical Reference Datasets (StRD). itl.nist.gov/div898/strd. Accessed 2026-09-06; public domain.
[17] NIST Chemistry WebBook, SRD 69. webbook.nist.gov. Accessed 2026-09-06; public domain (U.S. Government work).
[18] U.S. Geological Survey. National Water Information System (NWIS instantaneous values, 12 gauges, WY2023) and FDSN event service (global M ≥ 5.5, 2000–2024). earthquake.usgs.gov / waterdata.usgs.gov. Accessed 2026-09-06; public domain (U.S. Government work).
[19] Open-Meteo historical weather API (ERA5 reanalysis), 10 wind sites, 2023. open-meteo.com, CC BY 4.0 (ERA5: Copernicus licence). Accessed 2026-09-06; reanalysis, not in-situ measurement.
[20] OpenML. openml.org (datasets: ailerons, cpu_act, elusage). Accessed 2026-09-06; per-dataset licenses vary and are recorded in the released manifest.
[21] U.S. Department of Energy, fueleconomy.gov vehicle fuel-economy data, model years 1984–2025. Accessed 2026-09-06; public domain (U.S. Government work).
[22] RDatasets (Vincent Arel-Bundock), datasets: pressure, cars, rock, trees, attenu, housing, electricity, stormer, whiteside, cement, petrol, women, airquality, puromycin. vincentarelbundock.github.io/Rdatasets. Accessed 2026-09-06; licenses inherit from each source package (mostly GPL/BSD) and are recorded in the manifest.

## Appendix A — Suite composition

| Suite | Cells | Fits per cell | Total fits | Extra evaluations |
|---|---|---|---|---|
| S1 synthetic | 285 | 500 | 142,500 | — |
| S2 SRSD-Feynman | 480 | 500 | 240,000 | 240,000 test sets; 232,500 extrapolation probes |
| S3 family-stratified | 392 | 30 | 11,760 | — |
| S4 battery | 188 | 500 | 94,000 | — |
| S5 gate comparison (v3) | 102 | 200 | 20,400 | — |
| S6 real data + probe (v3) | 33 | 5 | 165 | — |
| **Total** | **1,480** | — | **508,825** | **472,500** |

Outside the seeded totals: the Appendix E engineering bench contributes
150 single-realization smoke cells (147 run), one fit per gate per cell.

## Appendix B — Failure registration taxonomy

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

## Appendix C — Headline Wilson 95% intervals

S2 H0, n = 200 (26,000 fits per cell; C/D on the verified SVD path):
A .053 [.050–.055]; B .050 [.047–.052]; C .053 [.051–.056]; D .056
[.054–.059]. S1 H0, n = 200 (4,500 fits per cell): A .037 [.032–.043];
B .041 [.036–.047]; C .045 [.039–.051]; D .042 [.036–.048]. Battery
power-family size: .0496 [.046–.053]. S5 H0 pools (600–2,600 fits per
cell): intervals within ±.013 of every Table 11 entry. Intervals for
all other cells are computable from released per-cell counts.

**Task-level uncertainty (S2, n = 200).** Seeds repeat within an
equation, so fit-count intervals understate cross-equation
uncertainty; the equation is the honest unit for generalization over
physical laws. Resampling the 13 H0 equations against the SVD-path
p-values (10,000 bootstrap replicates; `paper/review3_analysis.py`
§3.8a): C .053 [.050, .056] and D .056 [.055, .058], against A .0535
[.051–.056] and B .0497 [.048–.052] as reported above;
leave-one-equation-out moves C over [.052, .054] and D over
[.056, .057]. All four gates' cross-equation intervals are of the same
order on the corrected path — the wide interval belongs to the pinv
path (below), not to the test.

*The same two statistics under the retracted pinv path read C .417
[.197–.707] and D .563 [.276–.782], with leave-one-out ranges
[.368, .448] and [.527, .606]. This is the appendix's widest swing and
the clearest single demonstration that a task-level interval can be a
property of the numerical path rather than of the test: the
path-conditional width for C (.51) exceeds the SVD path's
cross-equation width (.02) by more than an order of magnitude. Under
that path C's per-equation rates split 5 equations at .998–1.00
(i.34.27, i.38.12, ii.27.18, iii.12.43, iii.15.14) against 8 at
.048–.060, and no single equation contributes more than 18.5% of its
false alarms — the failure was widespread across the near-collinear
subset, not one bad task. All method comparisons in this paper run on
the SVD path (Section 5, Q3).*

Macro and micro averages coincide here only because the design is
balanced (2,000 fits per equation); on unbalanced pools they should be
reported separately.

## Appendix E — The real-engineering data bench (companion to S6)

S6 probes mechanism transfer on macroeconomic aggregates. The deployment
domain these gates are actually built for is physical-engineering
measurement — fluid, civil, hydraulic, energy, meteorological, seismic,
chemical, and aerospace data, where power laws are the daily candidate
models. This appendix documents a companion bench of that kind, run under
the same contract as S6: real data has no ground truth, so the bench
registers compatibility and divergence, not error rates; its cells are
single realizations and stand outside Appendix A's seeded totals.

**Composition and discipline.** Nine open, license-free sources; 91
datasets; 150 candidate log-log relations preregistered in the released
manifest with per-relation physical notes; 761,515 raw rows (Table E1).
Raw bytes are cached (~115 MB) with SHA-256 hashes; cleaning is replayable
and imputation-free — rows with missing or non-positive values are dropped
and logged, and the few unit conversions are explicit in the manifest
(e.g., Whiteside °C → K). Simulator datasets that violate the real-data
discipline (a marine-engine and two wave-energy numerical models) are
excluded and listed, as are WebBook fluid pages below the 10-row
rejection-detection floor; one gauge whose instantaneous series was never
published is logged as an empty skip. ERA5 cells are marked as reanalysis,
not direct measurement.

**Table E1 — Composition by source.**

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
| **Total** | **91** | **150** | **761,515** | — |

**Protocol.** Each cell applies the §3.2 gates unchanged at α = 0.05 — A
with its B = 999 permutation p, B/C/D as registered — on the
positivity-filtered relation, subsampled to n ≤ 800 (seed 42) to bound
permutation cost. The run is a smoke pass: one realization per cell,
compatibility plus first verdicts, not a rate estimate. Of 150 cells, 147
run; three are registered `insufficient_n` by design (Bennett5 and
Roszman1 have no positive rows after the log-space positivity filter;
Danwood has n = 6) — the bench's honest-failure probes working as
intended.

**Table E2 — Smoke verdicts on the 147 run cells.**

| Gate | Rejections | Registered failures |
|---|---|---|
| A — permutation curvature | 97/147 (.660) | 10 declines, all `checked_false_small_n` (n ≤ 11) |
| B — RESET fitted F | 111/147 (.755) | 0 |
| C — RESET fitted HC3 | 108/147 (.735) | 0 |
| D — RESET exog HC3 | 94/147 (.639) | 3: `insufficient_df` (n = 13, four features), `hc3_leverage_degenerate` (n = 546), `rank_deficient` |

A and B agree on 123/147 cells (92 both reject, 31 both quiet; 5 A-only,
19 B-only). The failure columns are themselves informative: D's three
declines are exactly the three HC3 failure modes of Table 2 and Appendix
B, now on observational data; A's ten are uniform declines-to-run on
n ≤ 11 cells; B and C never fail. Nothing fails silently — every declined
evaluation is registered with a typed reason.

**Designed regimes.** The bench's design document assigns expected regimes
to 66 of the 150 cells before the smoke pass (Table E3). The assignment is
a reading aid, not a scored test: the H0-side laws are approximate by
design, and at n up to 800 an approximate law is a detectable null. The
12 rating-curve cells — whose registered form Q ~ h omits the documented
stage offset h₀ of the true rating law Q = C(h − h₀)ⁿ, a mismatch noted
in the manifest — are rejected unanimously by A and B alike (12/12):
detection of a documented misspecification, not a false alarm. The five
H0-side cells A retains are those with no designed departure mechanism —
the braking-distance v² law, the tree allometry, the road-load relation —
plus the two shear sites with the smallest fitted curvature.

**Table E3 — Designed-regime cells, rejections at α = 0.05.**

| Designed regime | Cells | A | B | C | D |
|---|---|---|---|---|---|
| Approximate power laws (H0-side) | 30 | 25 | 27 | 23 | 23 |
| Known non-power-law physics (H1-side) | 36 | 20 | 24 | 24 | 14 |

(A's H1-side count includes its nine WebBook declines at n ≤ 11 as
non-rejections.) On the H1 side no gate sweeps the board: the misses
concentrate in WebBook enthalpy cells whose measured temperature window
keeps the local power-law approximation good — e.g., toluene Δ_vapH ~ T:
all four gates quiet at n = 26 — a mutual-miss cell recorded as measured.

**Table E4 — Named cells: mechanism echoes and divergences; p-values at α = 0.05 in A / B / C / D order, — = declined with registered reason.**

| Cell | A / B / C / D | Reading |
|---|---|---|
| Wind shear v₁₀₀ ~ v₁₀, 10 ERA5 sites (n = 800 ea.) | 8/10 / 10/10 / 7/10 / 7/10 | B flags every site (p ≤ 2×10⁻⁴); A retains exactly the two smallest-curvature sites (.23, .41) — a graded response, not a threshold artifact |
| Gust factor ~ v₁₀, same sites (n = 800 ea.) | 10/10 ×4 | unanimous rejection at every site; no truth to adjudicate |
| puromycin rate ~ conc (Michaelis–Menten), n = 23 | **.011** / .069 / .169 / .169 | A detects the saturation curvature, B misses at α = .05 — Table 3's small-sample power gap, live |
| stormer Time ~ Viscosity + Wt (calibration law), n = 23 | .001 / 7×10⁻⁴ / 9×10⁻³ / — | single-feature cells quiet (A .94, B .99); the second feature carries the interaction the law contains |
| EPA CO₂ ~ MPG (stoichiometric inverse, exponent −1), n = 800 | .030 / 1.5×10⁻⁴ / .727 / .727 | A and B reject on a near-exact physical law bent by fleet drive-mix curvature |
| EPA combined MPG ~ displacement (road load), n = 800 | .386 / .596 / .599 / .599 | designed power law passes |
| CCPP PE ~ V, n = 800 | .642 / 8×10⁻⁴ / ≈0 / .023 | proxy curvature: V correlates .81 with the strongly curved AT, so B/C's fitted-power columns inherit AT's curvature while A's V-only direction set cannot — the mirror of dilution |
| automobile horsepower ~ engine size, n = 203 | .878 / 1.2×10⁻⁴ / .094 / .094 | divergence under heteroskedastic body-style scatter |
| rock perm ~ area + peri, n = 48 | .835 / 1×10⁻⁵ / 5.6×10⁻⁶ / .596 | A/D quiet, B/C reject; no truth to adjudicate |
| housing price ~ 4 features, n = 546 | .723 / .683 / ≈0 / — | D's leverage-degeneracy failure mode, live at n = 546, while A/B/C return verdicts |

**Scope.** What the bench supports is what S6 supports: compatibility at
scale (147/150 cells run across nine sources, zero silent failures),
failure registration on observational data, and mechanism-echo cells
whose directions match the synthetic map. Its per-gate rejection totals
(97/111/108/94 of 147) are compatibility statistics on single
realizations without ground truth — they rank nothing, and reading them
as performance would repeat exactly the leaderboard confusion the
benchmark exists to correct — A's calibrated conservatism (the
shear sites, the road-load law, puromycin at n = 23), A's interaction and
dilution behavior (stormer, CCPP), D's structural fragility (all three
failure modes live). What it cannot support: rates, calibration, or
rankings — single realizations, no ground truth, and time-ordered
NWIS/Open-Meteo series whose permutation-independence assumption is an
approximation marked in the manifest. All four p-values per cell are
released (`real_data/smoke_results.csv`) with the manifest
(`real_data/clean/manifest.json`, SHA-256 per source) and replayable
fetch/prepare/smoke scripts.

## Appendix G — Change notes (v2 → v5)

1. **Suite S5 added (gate comparison).** The question "why not just
   compare against a nested rival?" is promoted from an anticipated
   objection to a measured result: AICc-quad leads all four cells of
   Table 12; the permutation gate is repositioned from best-calibrated
   test to self-calibrating portfolio member. The repositioning is
   reflected in the abstract ("Fifth"), Q6, §8, and the title's framing
   ("gates" rather than "tests").
2. **Suite S6 added (real data).** World Bank macro relations with a
   preregistered list and a preregistered dilution probe; unanimous
   strong-signal cells and mechanism-consistent divergence, including a
   live reproduction of the DILUTION signature (Table 14). No rate claims.
3. **Totals updated.** 1,345 → 1,480 cells; 488,260 → 508,825 fitted
   datasets.
4. **S3 raw table dropped from body** (Table 5 in v2) to keep the table
   count within workshop limits; the raw family map remains in released
   per-cell data and the calibrated Table 5 carries the claim. v2's Q2
   and the "external validity" sentence in Q2 are superseded by the
   three-level dilution result (constructed → routed → observed).
5. **References extended** [12–14] (AICc, CV survey, World Bank data).
6. **Post-v3 addition: Appendix E (real-engineering bench).** A
   nine-source engineering companion to S6 — 91 datasets, 150
   preregistered relations, 761,515 rows — run as a single-realization
   smoke pass whose cells stand outside the seeded totals (unchanged):
   unanimous detection of a documented rating-curve misspecification,
   echoes of the small-sample power gap and the
   interaction-by-second-feature structure, all three registered HC3
   failure modes, and mutual misses on locally approximate laws.
7. **Major revision (v4) — reviewer response.** (a) The OR(A, E)
   combination's recall advantage is re-analyzed at the combination
   level: its naive joint false-positive budget is .065–.102, and
   recalibrated to a 5% joint budget it matches but never beats F
   (§4.7(d)); "best cell" claims withdrawn. (b) Fit-count Wilson
   intervals are supplemented with equation-level bootstrap intervals,
   leave-one-equation-out ranges, and per-equation rates (§4.1, Appendix
   C); the HC3 failure is shown to concentrate in 5 of 13 equations.
   (c) Permutation-null claims conditionalized throughout: residual
   permutation is an approximation calibrated empirically on our pools,
   not exact by construction (§3.2, §6). (d) **The HC3 claim is
   retracted and re-attributed**: an independent statsmodels check shows
   our pinv-based HC3 path — not the HC3 correction — produces the
   42–56% false alarms (1.00 vs .067 on the worst tasks; ≥3-digit
   agreement elsewhere); Tables 2/3/6 C/D rows are now labeled
   path-specific (Q3). (e) S5 calibration is validated out-of-sample by
   leave-one-equation-out threshold transfer (§4.7(c)). (f) Real-data
   claims downgraded to signature-consistency with confounds stated
   (§3.6, §4.8, Appendix E); S6's five seeds reframed as Monte Carlo
   replicates with measured spread. (g) Battery set-level errors and
   multiplicity reported (§4.6). (h) Related work extended to the wider
   specification-testing literature; references [15]–[22] split with
   access dates and per-source licenses; [8] description corrected;
   Table numbering restored to consecutive (v3 left a gap at 5).
8. **Second-round revision.** (a) Independent verification extended to
   exog-HC3 (D): all 13 H0 equations, stratified sampling, three-path
   comparison — our SVD path matches statsmodels (|Δp| ≤ .0023,
   covariance ≤ 3×10⁻⁴) while the X′X path fails only at measured
   cond(X′X) ≥ 4×10¹⁷; the numerics paragraph now states the exact
   call (NumPy 2.5.2 pinv, rcond = 10⁻¹⁵ relative) and replaces
   informal precision arithmetic with measured condition numbers (Q3).
   (b) The .067 re-verification restated as 8/120, Wilson [.034,.126] —
   "compatible with nominal", not precision calibration (Q3, §4.1).
   (c) Task-level paired bootstrap added for calibrated power: Δ(D−A),
   Δ(D−B), Δ(D−C) intervals with win fractions; "D strongest" softened
   to "highest average power in the current task pool" (§4.2,
   Appendix C). (d) Glyph coverage audited across the compiled PDF (zero
   uncovered characters; the reported missing-symbol boxes could not
   be reproduced). (e) S2 counts unchanged (32 tasks: 13 H0 + 19 H1; no task
   was removed; 1,480 cells and 508,825 datasets stand as audited).

9. **Goal-switch revision (post-v4): the entailment claim promoted to
   the paper's thesis.** The title, abstract, introduction, and
   conclusion now center the claim the evidence was built to test:
   fit quality and functional correctness are not connected by
   necessity, in either direction. New Section 4.10 (Table 15)
   quantifies both directions on the quality benches (57,000
   synthetic + 240,000 real-formula fits, C/D on the SVD path):
   within R² ≥ 0.99, true forms extrapolate at 0.1–0.2% error while
   wrong forms sit at 100% >10% (synthetic) and 9.7% (real); wrong
   forms below R² = 0.5 are jointly missed by all four gates at 57.5%
   (n=50), and ii.11.28 (R² ≈ 0.002) at 83%. New Q9 states the
   mechanism: detection tracks direction alignment
   corr(ŷ², residual), not fit magnitude. The gates are repositioned
   throughout as measuring instruments for the entailment question,
   not as competing detectors; all prior suite results, retractions,
   and audit chains stand unchanged. The analysis is regenerated by
   `no_causal_link.py` (twelve anchor checks, build-failing on drift).
10. **Review round 3 (six-item): bounded claims, joint release, SVD
   everywhere, conflict fixes.** (a) The thesis is re-bounded per
   review: "fit quality cannot by itself certify functional form;
   the evaluated specification gates remain jointly blind to some
   form mismatches" — the stronger "no entailment / no information /
   logical independence" wording is withdrawn (title, abstract,
   introduction, Q9, conclusion), and the analysis script's name is
   superseded by `svd_main_results.py` reporting the joint-release
   accounting. (b) The headline quantity is now the *joint event on
   the same models* (Section 4.10, Table 15): dangerous wrong forms
   (extrapolation error >10%) passing the R² ≥ 0.99 screen and
   jointly missed by all four gates — 0.67% of screened wrong forms,
   6.9% of dangerous ones, release rate of dangerous forms 7.8% →
   0.2% with n; certificate LR reported (dangerous wrong forms are
   65× less likely to hold a full-green certificate than true
   forms — the certificate is informative but non-excluding).
   Low-fit silences (ii.11.28) are repositioned as mechanism, not
   epidemiology — a model at R² ≈ 0.002 would not survive any sane
   fit screen. (c) All main-text C/D method comparisons now run on
   the verified SVD path (Tables 2, 3, 5, 15); the pinv(X′X) numbers
   survive only in Q3 as the numerical-implementation case study.
   This changes results substantively: D (exog-HC3, SVD) has the
   highest calibrated pooled power at every n on real formulas
   (cross-equation bootstrap confirming at n ≥ 200); the previous
   draft's "A leads at small and medium n" was an artifact of
   thresholds derived from artifact-contaminated p-mass. Table 5's
   trigonometric C/D entries similarly un-suppressed (.58/.86 under
   SVD). (d) Conflicts fixed: §4.6 empty-set wording now states the
   two-pool semantics (V never empty on member-truth cells; V = ∅
   is the correct outcome on out-of-family cells, at Table 10's
   rates); Table 12 excludes the OR row (own-5% budget ≠ per-gate
   5%) and reports it separately; §4.7(d) marks synth n=50 OR as
   uncalibratable (A alone fires at .080); the 805,825 total is
   replaced by the split accounting (508,825 fitted + 297,000
   re-evaluated). (e) Direction geometry is offered as a diagnostic,
   not an identified mechanism (circularity acknowledged; the
   independent experiment — vary alignment at fixed misfit — is
   named as future work).
9. **Third-round revision (reviewer response).** (a) **Denominators
   are explicit everywhere.** Table 15 previously printed C/B in its
   last column ("1.5% of screened / 13.2% of dangerous") while §4.10(b)
   and the abstract reported C/D ("7.8% → 0.2%") under the same label;
   the two ratios are now separate columns with their denominators
   named, and the n = 50 divergence (350 released of 2,655
   screened-dangerous, against 966 gate-silent of 12,369
   all-dangerous) is stated in the table note. (b) **The n-gradient
   claim is weakened.** §4.10(b) no longer reads "sample size closes
   the release channel": it reports both denominators, states that at
   n ≥ 200 every released dangerous model had already passed the
   screen, and names the three mixed effects (rising gate power, a
   changing numerator and denominator, the contracting best L₂
   approximation within the power-law class), with the decomposition
   left outstanding. (c) **Threshold conventions flagged.** The 10%
   danger bound and the R² ≥ 0.99 screen are labelled conventions in
   §4.10, and the screen-versus-gate budget asymmetry — an absolute
   goodness screen against a significance-level gate — is stated
   rather than left implicit. (d) **Certificate likelihood ratio
   corrected** to 65× (.845/.013; the earlier 50× came from rounding
   LR ≈ .02), in the abstract, §1, §4.10(c), §5, §8 and here. (e)
   **Unlike quantities are no longer juxtaposed**: the synthetic
   bench's 100% is labelled a design property of that bench, and 9.7%
   is named as the dangerous share of *screened* wrong forms. (f)
   **Appendix D is regenerated on the SVD path.** The previous
   draft's captions had drifted from the published images (D.1 was a
   CDF panel, not the retracted 42–56% claim; D.2 was the
   *uncalibrated* map; D.3 plotted S1 scenarios; D.4's ROC carried
   the pinv artifact), so the captions were first corrected to match
   the images and the panels were then re-rendered from the released
   cells (`paper/review3_analysis.py`, §3.7), with §7's bit-identical
   claim scoped to table cells and §6 recording the provenance. (g) **Appendix C
   re-derived on the SVD path** (C .053 [.051–.056], D .056
   [.054–.059]; task-level C [.045, .068]; leave-one-out
   [.048, .062]), with the pinv figures kept as the path-contrast
   case; the §4.1 sentence that attributed a pinv interval to the SVD
   path is corrected. (h) **Table 12's "matched FPR" claim scoped**:
   it holds in three of four columns; synth n = 50 is reported as
   measured because A fires at .080 there, the same standard the OR
   row is held to. (i) **§4.3 sensitivity readings downgraded** to
   directional — 30 seeds per cell cannot carry a 15-point trigger.
   (j) **The SR-method gap is stated as a gap** (Section 6): no
   symbolic-regression method is run anywhere in the benchmark, so the
   measured object is gate behavior on a pre-specified power-law
   candidate, not the pass rate of method-emitted wrong forms; that
   experiment is named as the next suite.
7. **Suite S7 added (joint-silence rate).** The panel's silence is now
   measured with the fit screen removed, on an independent re-sampling
   of S2's 18 non-monomial equations (2,000 datasets per task per n;
   36,000 per n). Joint silence — no gate rejects and no gate declines
   — is 37.8% / 29.8% / 23.2% naive and 38.6% / 30.2% / 23.6%
   H0-calibrated at n = 50 / 200 / 1000 (Table 16). This is a different
   quantity from §4.10's pipeline yield and the two are never pooled:
   §4.10 conditions on the R² ≥ 0.99 screen, S7 does not. Against the
   synthetic bench's 0.15–0.22% on the same unscreened basis the gap is
   170–250×. Adding the S5 model-comparison gates at matched
   calibration moves the n = 50 figure from 38.6% to 42.2%, against
   90.5% silence on true power laws (LR ≈ 2.2×; posterior ≈ 32% under a
   50/50 prior). **Provenance:** `runs/joint_silence/`; every number is
   regenerated from `joint_silence_srsd_per_task.csv`, whose values
   disagree with that directory's `REPORT.md` summary table (32.0% vs
   37.8%; 70% vs 90.5% for the H0 baseline — annotated in-file, October
   2026); the summary is not used anywhere, and no Wilson intervals are
   given because equations, not datasets, are the independent unit.
   Abstract, §4's table index, §6, and §8 are updated accordingly.
