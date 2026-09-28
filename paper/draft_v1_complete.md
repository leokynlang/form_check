# Silent Acceptance: A Falsification-Side Benchmark of Form-Specification Tests for Log-Linear Physical Model Discovery

**Complete draft v1 — working manuscript, not submitted.**

## Abstract

Symbolic regression and physical model discovery pipelines certify candidate models with goodness-of-fit metrics — R², cross-validation error, information criteria, coefficient standard errors — that measure how well a model fits, never whether an alternative model family would fit strictly better. A model can score near-perfectly on all of these metrics while having the wrong functional form, and fail catastrophically the moment it extrapolates. Classical form-specification tests exist — Ramsey's RESET and its built-regressor and exog-power variants, with heteroskedasticity-robust inference — but the model discovery community has no benchmark measuring whether these tests work on the data distributions such pipelines actually consume: physical laws sampled over physically realistic ranges with multiplicative noise, evaluated in log-transformed space.

We introduce a falsification-side benchmark filling this gap. Three task suites total 1,157 grid cells and 394,260 fitted datasets, with a further 472,500 held-out evaluations (240,000 same-distribution test sets and 232,500 extrapolation probes): (S1) 19 controlled scenarios with known ground truths; (S2) 32 tasks derived from the SRSD-Feynman registry of physical equations under the official per-variable sampling configurations; (S3) 30 Feynman equations stratified by structural family (exact power law, fractional power, exponential, rational, trigonometric, multivariable-additive). Four specification tests are compared under a single protocol: an identical significance level α = 0.05; a common permutation null for all tests; power compared only at a matched empirical false-positive rate of 5% estimated on true-monomial tasks; and explicit registration of every degenerate, rank-deficient, or non-finite evaluation.

Three findings stand out. First, aggregate rankings do not transfer across task families: the classical-F RESET variant that dominates synthetic quadratic benchmarks (power ≥ 0.99) is nearly blind to rational-function ground truths (0.13–0.19), because its fitted-value powers locally absorb exactly the curvature by which rational laws deviate from power laws. Second, heteroskedasticity-robust (HC3) Wald variants, which appear strongest in raw tables, carry false-positive rates of 42–56% on true monomial physical laws — flat in sample size from n = 50 to n = 1000 — so their apparent advantage disappears once methods are compared at matched empirical false-positive rates; the same signature (0.36–0.52) appears on the independent S3 pool. Third, no test is uniformly best: the permutation curvature test is the best-calibrated on both pools (5% ± 1%) but misses exponential-family laws (calibrated power 0.13) and has constructible blind spots (opposing curvature, power 0.02); each method's failure profile is complementary to the others'. Alarms enrich the probability of severe extrapolation error (precision 0.39–0.56 against a 0.20 base rate) but do not predict its magnitude, and the two worst extrapolators in S2 are missed by all four tests — no rejection is evidence of extrapolation safety. We release all generators, per-cell results, failure registrations, and content-hashed configurations under fixed seed rules.

## 1 Introduction

A regression model can be simultaneously excellent on every certification metric a discovery pipeline computes and wrong in its functional form. A power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes k-fold cross-validation, wins its information-criterion comparison, and reports small standard errors — while extrapolating with −20% error one decade outside the training domain. This is not a hypothetical: the failure was reproduced as a controlled scenario in our laboratory notebook and appears in this benchmark as scenario S1 (Section 3.3), where all four specification tests detect it at power ≥ 0.99. Fit quality and form correctness are different axes, and current discovery practice measures only the first.

The statistics literature has owned the second axis for over fifty years. Ramsey's RESET [1] tests whether powers of the fitted values carry additional signal; built-regressor variants use powers of the regressors; White [2] and MacKinnon–White [3] provide heteroskedasticity-robust inference. Yet these tests are absent from symbolic regression (SR) benchmarking. SRBench [5] measures whether discovery methods find the right expression across 252 problems; SRSD [6] redesigns Feynman-derived tasks with realistic sampling ranges and discovery-oriented metrics. In neither framework does anything measure whether an attached falsification layer would reject a wrong form. The asymmetry has practical consequences: a discovery pipeline with no working falsification layer silently accepts misfit forms, and its downstream use — extrapolation, engineering design, scientific claim-making — inherits the error silently.

Measuring the falsification side faithfully requires resolving a protocol problem that, to our knowledge, no prior work addresses for this test family: the tests reject on incommensurable bases. A bootstrap-t threshold ("alarm when |t| > 3") and an asymptotic F p-value ("reject when p < 0.05") measure different things; comparing them directly measures threshold placement, not statistical power. Our first pilot did exactly this and produced inverted rankings that vanished under a common basis. The protocol below (Section 3.1) therefore rests on three commitments: one significance level for all tests; one permutation null shared by all tests; and power compared only at a common empirical false-positive rate estimated on true-monomial tasks.

The question this benchmark answers is narrow and falsifiable: **when the ground truth is known and the candidate form is wrong, which specification tests detect the misspecification — and at what false-alarm cost, on data distributions that physical model discovery pipelines actually consume?** We contribute: (i) the benchmark itself — three suites, 394,260 fitted datasets plus 472,500 held-out evaluations, released generators with seed rules that make any single cell independently reproducible; (ii) the protocol — matched-basis comparison with empirical-FPR calibration and explicit failure registration; and (iii) the resulting failure-mode map — including two previously undocumented findings: the flat-in-n 42–56% false-positive rates of HC3-Wald RESET variants on real physical monomials, and the near-blindness of fitted-F RESET to rational laws.

## 2 Related Work

**Specification testing.** RESET [1] augments the regression with powers of the fitted values and tests the joint significance of the added coefficients; the exog variant uses powers of the regressors instead. Heteroskedasticity-consistent covariance estimation [2] and its HC3 leverage-weighted refinement [3, 4] are the standard small-sample corrections and are the default in several mainstream econometrics and statistics packages. These results are derived for linear models under fixed or quasi-fixed designs; none of the classical treatments contemplates log-transformed physical laws whose regressors span only a few orders of magnitude — precisely the regime our benchmark shows breaks HC3's small-sample behavior (Section 4.1).

**Symbolic regression benchmarks.** SRBench [5] established the discovery-side evaluation: 14 SR methods and 7 ML baselines on 252 problems, scored by exact recovery and held-out accuracy. SRSD [6] reconceived the task distribution itself — realistic per-variable ranges reviewed against the physics, rather than toy domains — and is the source of our task configurations. Both measure the discovery side; a pipeline that returns a wrong expression and one that flags it score identically. Our generators deliberately reuse SRSD's official sampling objects so that the falsification-side pools are the same distributions the discovery side already standardizes on.

**Physics-informed and dimensionally-aware discovery.** AI Feynman [7] and data-driven dimensional analysis [8] use physical units to constrain the hypothesis space; PySR [9] represents the modern search side. Dimensional correctness is necessary but not sufficient for form correctness: every ground-truth law in our suites is dimensionally self-consistent, yet each non-monomial law is a form misspecification for a power-law model. The falsification side complements these constraints — units narrow what to search; specification tests check what was found.

**Honest-failure evaluation.** Existing SR benchmarks report no analogue of a false-positive rate because they contain no explicit rejection mechanism to calibrate. Work on selective prediction and abstention addresses related epistemics for classification; we are not aware of a symbolic regression benchmark that measures attached falsification behavior, per family, with preregistered labels and failure registration. This is the gap the present benchmark occupies.

## 3 Benchmark Design

### 3.1 Design principles

**P1 — Common rejection basis.** All four tests decide at the same significance level α = 0.05 on p-values. The permutation curvature test's empirical p-value is computed from B = 999 permutations of the residual-permutation null that its native bootstrap already uses; no test is judged by a heuristic threshold.

**P2 — Power only at matched empirical FPR.** Headline power comparisons calibrate each method's decision threshold to a 5% empirical false-positive rate, estimated on pooled true-monomial tasks within the same suite, sample-size stratum, and noise regime. Raw-α results are reported alongside, because raw tables are what current practice would see.

**P3 — Preregistered labels.** H0 (strict monomial: the log-log mean function is exactly linear) versus H1 (any non-monomial structure) is fixed from symbolic structure of the ground truth before any test runs. No label is revised after seeing results; the family stratification of S3 is likewise fixed from expression structure.

**P4 — Explicit failure registration.** Every evaluation that cannot produce a valid decision is registered with a typed reason — zero-residual degeneracy, design rank deficiency, insufficient residual degrees of freedom, leverage degeneracy (1 − h → 0), singular covariance, statistic mismatch between wrapper and native implementation — and excluded from rate numerators while remaining in the denominators' audit trail. The zero-noise grid cells exist to exercise this path: all four tests register degeneracy there rather than emitting arbitrary p-values.

**P5 — Official data configurations.** S2 and S3 tasks reuse the official per-variable sampling objects of the source registry: logarithmic domains for DefaultSampling variables, linear domains for SimpleSampling variables, integer domains where specified, and the official positive/negative sign flags. Tasks whose inputs or outputs are not strictly positive under the official configuration are marked not-applicable and excluded — we do not shift or rescale data to make them log-fittable, since that changes the hypothesis under test. Sensitivity to the official ranges is measured explicitly (Section 4.3) rather than assumed.

### 3.2 Methods under test

All tests operate in log-log space: given positive data, let ly = log y and X = [1, lx₁, …, lxₖ] be the design of the candidate power-law model. The tests differ only in what they add to this design and how they test the addition.

**A. Permutation curvature (form_check).** Adds a single mean-centered quadratic column q = meanⱼ((lxⱼ − lx̄ⱼ)²) — for one feature, the squared centered log-feature itself; for several, their row-wise mean — and takes the added coefficient c₂ as statistic. The native implementation reports a t-statistic from 200 residual-permutation draws; for the benchmark, the wrapper recomputes the identical statistic and derives an empirical p-value from B = 999 permutations, using a vectorized formulation verified to machine precision against the loop implementation. A statistic-mismatch sentinel compares wrapper and native values on every call.

**B. RESET fitted, classical F.** Adds ŷ² and ŷ³, where ŷ is the fitted value from the restricted model, and tests their joint significance with the classical F statistic [(RSS_r − RSS_u)/q] / [RSS_u/(n − k)], powers (2, 3) throughout.

**C. RESET fitted, HC3-Wald.** Same augmented design; Wald statistic W = β_u' [V_HC3]⁻¹ β_u on the added coefficients, where V_HC3 = (X'X)⁻¹ [Σᵢ eᵢ²/(1−hᵢ)² xᵢxᵢ'] (X'X)⁻¹, referred to F(q, n − k).

**D. RESET exog, HC3-Wald.** Built regressors are the powers of each log-regressor (lxⱼ², lxⱼ³ for every feature j) rather than of the fitted values; HC3-Wald as in C.

Configurations were frozen before the main runs — powers (2, 3); B = 999; α = 0.05; the grids of Section 3.3 — and are identical across all three suites. The same harness invokes all four tests on every dataset, so method comparisons share data, seeds, and failure handling exactly.

### 3.3 Task suites

**S1 — Synthetic control (19 scenarios × 500 seeds).** Three monomial H0 scenarios (univariate; two-variable; two-variable with ρ ≈ 0.9 correlated design), each with strictly positive log-uniform features on [0.1, 10] × [0.05, 2]. Eleven structural H1 scenarios: quadratic curvature (c = 0.10); cubic (standardized odd term); multiplicative interaction (log-space cross term); threshold (exponent 0.4 → 1.0 at the domain midpoint); exponential law; logarithmic law; rational (saturating) law; omitted variable (a correlated, quadratically-entering covariate absent from the fitting surface); dilution (one of two variables curved); opposing curvature (+0.12, −0.12 on an exchangeable design, so the mean-square column's population coefficient is exactly zero); collinear-quadratic (ρ ≈ 0.95 design with curvature). Five noise-layer H1 scenarios have strict-monomial mean structure but non-Gaussian noise — Student-t(3); standard deviation growing linearly with the log prediction (0.3× to 1.7× the nominal); 2% point outliers (+8σ); 6% group outliers (+5σ contiguous); AR(1) correlation (ρ = 0.6) — rejections there are false alarms by construction and are reported as such, never pooled into structural power. Grid: n ∈ {50, 200, 1000}; relative Gaussian noise σ ∈ {0, 0.5%, 1%, 3%, 5%} injected additively in log space (equivalently multiplicatively in data space); 500 seeds per cell via a seed rule that embeds the grid position, making any cell independently reproducible.

**S2 — SRSD-Feynman tasks (32 tasks × 500 seeds).** From the 120-equation registry, 32 tasks pass the positivity screen under official configurations: 13 strict-monomial H0 (including physical-constant domains such as Planck constants at 10⁻³⁴–10⁻²³ scale) and 19 non-monomial H1. Noise is injected in log space at the S2 grid (excluding the zero-noise level for extrapolation probes). Each training fit (n ∈ {50, 200, 1000}) is accompanied by a same-distribution test set of 2,000 samples and a far-field probe: the first regressor extrapolated ×3 beyond both training domain edges (log-space shift ln 3), other regressors at training medians, scored against the closed-form ground truth in a numerically stable form (exp(log-ratio) − 1).

**S3 — Family-stratified Feynman equations (30 equations × 30 seeds).** Structural family classification from symbolic form: exact power law (12 H0); fractional power (1 H0, exponent 1/3); exponential (3 H1); rational (6 H1); trigonometric (5 H1); plus multivariable-additive, logarithmic, and other families represented only by not-applicable members (4 variables or non-positive under official ranges), which appear in the selection table with explicit reasons but are not tested. Each H1 equation runs over relative noise {0.5%, 1%, 3%, 5%} × {homoskedastic, heteroskedastic (SD growing 0.3×→1.7× with log prediction)} × sampling-range scale {official, ×0.5, ×2}; H0 equations run the noise grid at official range. Total across suites: 1,157 cells, 394,260 fitted datasets, plus 472,500 held-out evaluations — 240,000 same-distribution test sets (one per training fit) and 232,500 extrapolation probes (the remainder exclude one equation, i.26.2, whose ground truth is not finite on the ×3 probe domain; those probes are registered as not-a-number rather than silently scored).

## 4 Results

Table numbering is consecutive; Tables 1–2 report size, Tables 3–5 power, Table 6 noise robustness, Table 7 the family map, Table 8 extrapolation. All rates are rejection fractions; Wilson 95% intervals for headline rates appear in Appendix C. Seed counts per cell: 500 (S1, S2), 30 (S3).
### 4.1 Size: false-positive rates on true monomials

**Table 1 — S1 synthetic H0, false-positive rate at α = 0.05 (Gaussian noise pooled over σ ∈ {0.5%…5%}; 4,500 fits per cell).**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A — permutation curvature | .047 [.041–.053] | .037 [.032–.043] | .041 [.036–.047] |
| B — RESET fitted F | .045 [.039–.051] | .041 [.036–.047] | .038 [.033–.044] |
| C — RESET fitted HC3 | .048 [.042–.054] | .045 [.039–.051] | .038 [.033–.044] |
| D — RESET exog HC3 | .054 [.048–.061] | .042 [.036–.048] | .040 [.034–.046] |

**Table 2 — S2 true physical monomials (13 tasks), false-positive rate at α = 0.05 (26,000 fits per cell).**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .049 [.046–.052] | .053 [.050–.055] | .052 [.049–.054] |
| B | .050 [.047–.053] | .050 [.047–.052] | .053 [.050–.055] |
| C | **.422 [.416–.428]** | **.417 [.411–.423]** | **.418 [.412–.424]** |
| D | **.564 [.558–.570]** | **.563 [.557–.569]** | **.565 [.559–.571]** |

On synthetic Gaussian data, all four tests are calibrated (Table 1): every rate falls within or near its nominal band, and no method is distinguishable from another. On real physical laws the picture splits (Table 2). A and B remain calibrated — rates .049–.053, intervals overlapping 0.05, flat in n. Both HC3 variants break, and break *identically at every sample size*: C stabilizes at .417–.422, D at .563–.565, with 95% intervals that do not move from n = 50 to n = 1000. This is not a small-sample artifact that more data cures; it is a structural miscalibration.

The mechanism is collinearity of the built columns. Physical constants and engineered quantities span few log-decades; on such designs ŷ² and ŷ³ are nearly linear functions of the fitted values themselves (for C), and the built power columns for D are nearly constant or near-collinear with the design — S3 records the limiting case, a ⅓-power monomial whose built columns are effectively constant and where D's false-positive rate reaches 1.00. The sandwich estimator's asymptotic guarantee presumes the design matrix is well-behaved as n grows; here the design degeneracy is exact and does not wash out, so the Wald statistic's reference distribution is simply wrong. H0 exact-power-law cells of S3 reproduce the same signature at independent scale: C .52, D .36.

Noise-layer scenarios (Table 6, Section 4.5) confirm what does *not* cause this: heavy tails, growing variance, outliers, and autocorrelation leave all methods' sizes at .03–.07. Violations of the noise assumption are not what breaks specification tests on real data — the design is.

### 4.2 Power: raw, then calibrated

Raw power at α = 0.05 on S2 H1 (n = 200): A .52, B .28, C .63, D .94. A naive reading crowns D. Table 2 already priced that crown: D pays a .563 false-positive rate for it. Calibrating every method to a common 5% empirical FPR on the H0 pool (per sample size) reorders the field:

**Table 3 — S2 power at matched empirical FPR 5%.**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | **.41** | **.51** | .59 |
| B | .19 | .28 | .42 |
| C | .05 | .28 | .41 |
| D | .19 | .37 | **.65** |

At deployment-relevant sample sizes (n ≤ 200), A delivers roughly twice B's calibrated power and 2–10× C's and D's. At n = 1000, D overtakes (.65): with enough data, exog powers resolve the very structures that the mean-square average compresses, and its calibration failure becomes less damaging relative to its information gain. Even so, no method exceeds .65 — most of S2's H1 mass is *hard* for every test, which Section 4.3 decomposes by family.

On S1 structural H1 (raw α, rel = 1%, n = 200), the same pattern appears in sharper form. Quadratic, threshold, exponential-shape, logarithmic, rational-single, omitted-variable, and collinear-quadratic scenarios are caught by all methods at ≥ .96. The informative cells are the blind spots:

**Table 4 — S1 blind-spot cartography (raw α, n = 200, rel = 1%).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| cubic (standardized odd term) | .40 | 1.00 | 1.00 | 1.00 |
| interaction (lx₁·lx₂ cross term) | .21 | 1.00 | 1.00 | .10 |
| opposing curvature (+c, −c, exchangeable) | **.02** | .98 | .98 | 1.00 |
| dilution (one curved variable of two) | 1.00 | .88 | .87 | 1.00 |

Three mechanisms, each verified by construction. (i) A's statistic is driven by the mean of squared centered columns; when two variables carry equal and opposite curvature on an exchangeable design (OPPOSING), that mean has exactly zero population signal — A's .02 at n = 200 and .00 at n = 1000 is the statistic doing precisely what it is defined to do. (ii) A quadratic-projection statistic loses power on odd-symmetric deviations: A's cubic power (.40–.47 across n) never approaches B's saturation, because the cubic's second-order projection is small relative to its third-order content. (iii) B's fitted-value powers are a local Taylor surrogate for *any* smooth alternative — which is why B saturates on interaction (1.00) where A (.21) and D (.10) struggle — and, as S3 will show, the same surrogate property near-aborts B on rational laws.

### 4.3 Family-stratified results and sensitivity

**Table 5 — S3 family map, official sampling domains, homoskedastic noise pooled over rel (raw α).**

| Family (class) | A | B | C | D |
|---|---|---|---|---|
| exact power law (H0, n=13 eq. → size) | .04 | .04 | .52 | .36 |
| fractional power ⅓ (H0 → size) | .06 | .09 | .10 | 1.00 |
| exponential (H1, 3 eq.) | .11 | .05 | .70 | .68 |
| rational (H1, 6 eq.) | .66 | .13 | .68 | .95 |
| trigonometric (H1, 5 eq.) | .69 | .72 | .81 | 1.00 |

Calibrated to 5% empirical FPR (S3 H0 pool), the family map becomes the paper's sharpest statement:

**Table 5b — S3 calibrated power by family.**

| Family | A | B | C | D |
|---|---|---|---|---|
| exponential | .11 | .05 | **.63** | .04 |
| rational | **.66** | .14 | .50 | .66 |
| trigonometric | **.70** | **.73** | .08 | .05 |

Read jointly, Tables 5 and 5b are the benchmark's central finding. B — the classical test, near-perfect on synthetic quadratics — is nearly blind to rational laws (.13 raw, .14 calibrated): the fitted-value Taylor surrogate absorbs precisely the low-order curvature by which a saturating law departs from a power law. D's raw dominance on rational (.95) and trigonometric (1.00) is false-alarm budget, not power: calibrated, its trigonometric advantage collapses to .05. C is the exponential family's only reliable detector (.63 calibrated) while being worthless on trigonometric (.08). A misses exponentials (.11) — an exponential law's log-log curvature is modest at official scales — but is otherwise the most even performer. Every method's strength sits on another method's blind spot; the four columns of Table 5b have no row-major winner.

**Sensitivity.** Scaling every sampling domain ×2 moves exactly one family-method cell by more than 15 percentage points (rational × A: .66 → .84 — larger domains make curvature easier to see, the expected direction); ×0.5 moves none. The heteroskedastic replication of every S3 cell agrees with the homoskedastic one within .02 on every family-method cell. Sampling-range choice and variance model — the two design knobs a practitioner might suspect of driving the conclusions — do not drive them.

### 4.4 Noise-layer robustness

**Table 6 — S1 noise-layer H1 (strict-monomial mean structure) rejection rates at α = 0.05; these are false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| Student-t(3) noise | .034 | .036 | .030 | .030 |
| variance growing with prediction | .048 | .058 | .050 | .050 |
| AR(1) errors, ρ = 0.6 | .038 | .040 | .042 | .042 |
| 2% point outliers (+8σ) | .040 | .066 | .036 | .036 |
| 6% group outliers (+5σ) | .038 | .060 | .054 | .054 |

Every method, every non-Gaussian regime: .030–.066. Heavy tails, heteroskedasticity, serial correlation, and outliers — the standard rogues' gallery of regression diagnostics — do not push any of these tests above 7% false alarms. This table's function is exclusionary: it rules out the convenient explanation that the miscalibrations of Table 2 come from noisy, ugly, real-world data. The data in Table 2 are clean by construction; what differs from Table 1 is the *design* — physical constants and few-decade quantities. The design, not the noise, is what specification tests are fragile to.

### 4.5 Do alarms predict extrapolation failure?

Each S2 H1 fit is scored on 2,000 same-distribution test samples and on a ×3-domain far-field probe against the closed-form truth. Base rate of far-field relative error above 10%: .198. Specification-test alarms enrich this substantially:

**Table 7 — Alarm precision and recall for far-field error > 10% (S2, all n, rel ≥ 0.5%).**

| Method | precision | recall |
|---|---|---|
| A | .39 | .63 |
| B | .56 | .49 |

Two properties matter more than the margins. First, enrichment: an A alarm raises the probability of severe extrapolation failure from .198 to .39 (2.0×); a B alarm to .56 (2.8×). Second, complementarity: of the true-positive catches, 10,560 are A-only and 5,447 B-only, against 12,575 shared — the two tests' alarm sets overlap barely more than chance, so a portfolio of both captures substantially more failures than either. Alarms are therefore useful as *ordinal evidence* about form.

They are not useful as predictions of magnitude. P-value rank and error rank correlate only weakly within H1 (Spearman ρ = 0.25 for A across tasks), because rejection probability saturates long before extrapolation error does: a law can be decisively rejected in-domain yet extrapolate benignly (bonus.15: power 1.00, far-field error 9%), and the suite's two worst extrapolators — i.39.11 and i.43.43, both with far-field error 98% — carry near-zero in-domain curvature and are missed by all four tests (A power .06, B .27–.29). An alarm means form evidence exists; a non-alarm means nothing about the field beyond the training domain. No rejection is evidence of extrapolation safety.

### 4.6 Runtime

Per-dataset wall time (S1, all scenarios, means over 500-seed cells): A 12.0 ms (n = 50) to 25.3 ms (n = 1000) — the cost of 999 permutation evaluations, vectorized; B/C/D 1.0–1.8 ms across all n. All four tests are effectively free at the scale of a discovery pipeline that already fits thousands of candidates; the bottleneck resource is not compute but the *decision discipline* — the calibration set of Section P2 requires only true-monomial tasks that a physics repository already contains.

## 5 Analysis

**Q1. Does the low-df curvature test beat RESET on univariate quadratic deviations?** No — it ties. On S1's quadratic truth both saturate at power 1.00 from n = 50. What low degrees of freedom buy is calibration, not power: A stays within ±1 point of nominal size on every H0 pool in this paper, which is what makes its calibrated power numbers *trustworthy* rather than large. The asymmetries that do appear are structural: odd-symmetric cubic deviations (A .40–.47 vs B 1.00) and interaction terms (A .21 vs B 1.00, D .10), where A's quadratic projection is mathematically the wrong instrument.

**Q2. Does averaging squared columns miss feature-specific curvature and interactions?** Yes, by construction — and the benchmark quantifies when that matters. Cancellation is real: opposing curvature drives A to .02/.00. Dilution is not: one curved variable among two still yields A 1.00 versus B .88. Interactions blind A (.21) but saturate B. The decisive observation is external validity: neither the cancellation nor the interaction structure appears among the 32 real physical laws, where A's worst family is exponential (.11) and its blind spots are empty. Synthetic blind spots are true properties of the statistic; only a real-equation pool can reveal whether they bind in practice. Both suites are needed for that sentence to be sayable.

**Q3. Does HC3 fix false alarms under heteroskedasticity?** Heteroskedasticity was never the disease here. Synthetic growing-variance noise leaves all methods at .048–.058 (Table 6), and every S3 heteroskedastic cell replicates its homoskedastic counterpart within .02. What breaks on real physical monomials is HC3 itself (Table 2): the robustness correction presumes a design whose leverage structure is informative; on near-collinear built columns the leverage values are uninformative, and the correction multiplies noise into a statistic whose reference distribution is wrong at every n measured. The robust estimator is the failure mode — an inversion of its textbook reputation, documented at 26,000 fits per cell.

**Q4. Do alarms predict extrapolation error?** As binary evidence, usefully: 2.0–2.8× enrichment of severe-error probability, with complementary alarm sets (Section 4.5). As magnitudes or guarantees, no: rank correlation within H1 is ρ = 0.25, and the worst extrapolators in the suite carry the least in-domain curvature. The correct reading is asymmetric: alarm ⟹ form evidence (and elevated extrapolation risk); silence ⟹ nothing. Pipelines that treat a passed form check as an extrapolation certificate are making an error this benchmark can now price: in S2, that error rate is the .198 base rate itself.

**Portfolio reading.** The four columns of Table 5b have no dominant row; the alarm sets of Table 7 barely overlap; A and B disagree on cubic, interaction, and rational in opposite directions. The failure profiles are complementary enough that portfolio deployment — any alarm triggers model-family review — strictly dominates any single test, at a combined runtime under 30 ms per dataset. This is the deployment recommendation the benchmark supports, and it is a *falsification-side* recommendation: it concerns when to distrust a candidate, not how to find a better one.

## 6 Limitations

**Single candidate family.** All four tests evaluate a power law against alternatives; nothing here benchmarks model *selection* among families, only rejection of one. **Fixed augmentation.** Built-regressor powers are fixed at (2, 3); adaptive or higher-order constructions may shrink B's rational blind spot at calibration cost we have not measured. **Log-space scope.** All tests operate in log-log space by design; additive-error models are out of scope and marked not-applicable rather than log-shifted, which defines the benchmark's boundary rather than a finding. **Dimensional manipulation absent.** Units were recorded but not experimentally manipulated; the observation that a degenerate ⅓-power monomial drives D's false alarms to 1.00 suggests structure-aware screening could pre-empty some failures, but this is untested here. **Statistical precision.** Headline S1/S2 rates rest on 500 seeds (Wilson intervals ±.003–.006 at p ≈ .05, ±.02 near p ≈ .5); S3 family cells rest on 30 seeds and support directional claims only. **One noise draw per seed.** Datasets are not paired across methods by noise realization beyond shared seeds — shared seeds do pair them — but no cross-seed variance decomposition is attempted.

## 7 Reproducibility statement

Every generator is seeded from a documented rule embedding the grid position (suite, scenario, n, noise level, seed index), so any single cell of any table regenerates bit-identically. The harness records content hashes of the code under test — including the exact `form_check` implementation — and of the harness itself, in every run's environment file; a statistic-mismatch sentinel verifies on every call that the wrapper's recomputation agrees with the native implementation, so drift in the tested code cannot pass silently. All 1,157 cells' raw results (parquet and CSV), per-equation tables, failure registrations with typed reasons, selection and exclusion tables with per-equation justifications, analysis scripts, and figure sources are released. Zero-noise grid cells deliberately exercise the failure-registration path: all four tests register degeneracy there rather than emitting p-values, and those registrations are part of the released record.

## 8 Conclusion

We built the falsification-side benchmark that model discovery lacks: three suites, 394,260 fitted datasets and 472,500 held-out evaluations, one protocol that puts four classical specification tests on a common basis, and full registration of every failure. Its findings are specific and actionable. Robust-covariance RESET variants — the field's default "safe" choice — are structurally miscalibrated on real physical laws (42–56% false alarms, flat in n), an inversion of their textbook reputation documented here at a scale no single-dataset study could reach. Raw power rankings invert under matched false-positive rates, and invert again across task families: the test that dominates synthetic benchmarks is blind to rational laws; the test that is blind to interactions saturates on them synthetically. No test wins everywhere; the failure profiles are complementary; a portfolio of two tests at under 30 ms per dataset dominates any single one. And the epistemic boundary is now priced: alarms enrich severe-extrapolation probability 2–3×, silence certifies nothing. Benchmark papers are usually measured by the leaderboards they produce. This one is measured by the failures it makes visible — because in scientific discovery, the model you silently accept is the one that breaks your extrapolation.

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

## Appendix A — Suite composition

| Suite | Cells | Fits per cell | Total fits | Extra evaluations |
|---|---|---|---|---|
| S1 synthetic | 19 × 3n × 5σ = 285 | 500 | 142,500 | — |
| S2 SRSD-Feynman | 32 × 3n × 5σ = 480 | 500 | 240,000 | 240,000 same-distribution test sets; 232,500 extrapolation probes (7,500 not finite on probe domain, registered) |
| S3 family-stratified | 14 H1 × 4σ × 2 noise × 3 scales + 7 H0 × 4σ × 2 noise = 392 | 30 | 11,760 | — |
| **Total** | **1,157** | — | **394,260** | **472,500** |

*Correction note (v1): an earlier draft claimed "643,250 fitted datasets" and "1,167 cells"; the audited totals are 394,260 fits over 1,157 cells, plus 472,500 held-out evaluations (240,000 test sets; 232,500 valid extrapolation probes — 7,500 probes on equation i.26.2 return non-finite ground truth on the ×3 domain and are registered, not scored). This version supersedes the earlier figure everywhere; the annotated v0 draft retains it as an erratum example.*

## Appendix B — Failure registration taxonomy

| Reason | Meaning | Exercised by |
|---|---|---|
| `degenerate_zero_residual` | Zero-noise cell on exact ground truth; no residual variation → test undefined | S1/S2 σ = 0 cells |
| `rank_deficient` | Design matrix numerically rank-deficient | collinear designs |
| `insufficient_df` | n − k < 1 after augmentation | n = 50 multi-feature cells |
| `hc3_leverage_degenerate` | 1 − hᵢ ≤ 10⁻¹²; HC3 undefined | degenerate built columns |
| `hc3_singular_covariance` | V_HC3 numerically singular | degenerate built columns |
| `statistic_mismatch` | Wrapper statistic ≠ native statistic (drift sentinel) | continuous audit |
| `checked_false_*` | Native test declines to run (non-positive data, n < 12, non-finite) | not-applicable tasks |

All registered failures remain in the audit trail; rates' denominators are documented per table. Zero-noise cells intentionally trigger `degenerate_zero_residual` for all four methods — the registration path is tested by the benchmark itself.

## Appendix C — Headline Wilson 95% intervals

S2 H0, n = 200 (26,000 fits per cell): A .053 [.050–.055]; B .050 [.047–.052]; C .417 [.411–.423]; D .563 [.557–.569]. S1 H0, n = 200 (4,500 fits per cell): A .037 [.032–.043]; B .041 [.036–.047]; C .045 [.039–.051]; D .042 [.036–.048]. Intervals for all other cells are computable from released per-cell counts.
