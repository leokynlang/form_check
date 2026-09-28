# Silent Acceptance: A Falsification-Side Benchmark of Form-Specification Tests for Physical Model Discovery

**Complete draft v2 — working manuscript, not submitted.** Supersedes `draft_v1_complete.md`. Changes from v1: the family-viability battery is merged (§3.4, §4.6, Q5) with numbers refreshed from the formal 500-seed run (`runs/battery_full`; 94,000 datasets, zero data-layer failures); two errors in v1's extrapolation section are corrected (§4.5 and the change note in Appendix D); the table-numbering preamble is fixed.

## Abstract

Symbolic regression and physical model discovery pipelines certify candidate models with goodness-of-fit metrics — R², cross-validation error, information criteria, coefficient standard errors — that measure how well a model fits, never whether an alternative model family would fit strictly better. A model can score near-perfectly on all of these metrics while having the wrong functional form, and fail catastrophically the moment it extrapolates. Classical form-specification tests exist — Ramsey's RESET and its built-regressor and exog-power variants, with heteroskedasticity-robust inference — but the model discovery community has no benchmark measuring whether these tests work on the data distributions such pipelines actually consume: physical laws sampled over physically realistic ranges with multiplicative noise, evaluated in log-transformed space.

We introduce a falsification-side benchmark filling this gap. Three task suites total 1,157 grid cells and 394,260 fitted datasets, with a further 472,500 held-out evaluations (240,000 same-distribution test sets and 232,500 extrapolation probes): (S1) 19 controlled scenarios with known ground truths; (S2) 32 tasks derived from the SRSD-Feynman registry of physical equations under the official per-variable sampling configurations; (S3) 30 Feynman equations stratified by structural family. Four specification tests are compared under a single protocol: an identical significance level α = 0.05; a common permutation null for all tests; power compared only at a matched empirical false-positive rate of 5% estimated on true-monomial tasks; and explicit registration of every degenerate, rank-deficient, or non-finite evaluation. A fourth suite (S4, 188 cells, 94,000 datasets) extends the same machinery from rejection to routing: a four-family viability battery that returns, per dataset, the set of linearizable families not rejected.

Four findings stand out. First, aggregate rankings do not transfer across task families: the classical-F RESET variant that dominates synthetic quadratic benchmarks (power ≥ 0.99) is nearly blind to rational-function ground truths (0.13–0.19), because its fitted-value powers locally absorb exactly the curvature by which rational laws deviate from power laws. Second, heteroskedasticity-robust (HC3) Wald variants, which appear strongest in raw tables, carry false-positive rates of 42–56% on true monomial physical laws — flat in sample size from n = 50 to n = 1000 — so their apparent advantage disappears once methods are compared at matched empirical false-positive rates. Third, no test is uniformly best: the permutation curvature test is the best-calibrated on both pools (5% ± 1%) but misses exponential-family laws and has constructible blind spots; each method's failure profile is complementary to the others'. Alarms enrich the probability of severe extrapolation error (precision 0.39–0.56 against a 0.20 base rate) but do not predict its magnitude, and no rejection is evidence of extrapolation safety. Fourth, the battery shows that the same statistic can route between families — but only after its null distribution is rebuilt per family, because null validity is a coordinate-dependent property: the permutation null that is exact in log–log coordinates produces 30–34% false escalations in the rational family's (1/x, 1/y) coordinates, and a wild bootstrap repairs it to 4.5–9%, delivering containment of 93–95% on in-family truths and 1.00 correct escalation on the structures furthest outside all four families (quadratic, threshold, omitted-variable). We release all generators, per-cell results, failure registrations, and content-hashed configurations under fixed seed rules.

## 1 Introduction

A regression model can be simultaneously excellent on every certification metric a discovery pipeline computes and wrong in its functional form. A power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes k-fold cross-validation, wins its information-criterion comparison, and reports small standard errors — while extrapolating with −20% error one decade outside the training domain. This is not a hypothetical: the failure was reproduced as a controlled scenario in our laboratory notebook and appears in this benchmark as scenario S1 (Section 3.3), where all four specification tests detect it at power ≥ 0.99. Fit quality and form correctness are different axes, and current discovery practice measures only the first.

The statistics literature has owned the second axis for over fifty years. Ramsey's RESET [1] tests whether powers of the fitted values carry additional signal; built-regressor variants use powers of the regressors; White [2] and MacKinnon–White [3] provide heteroskedasticity-robust inference. Yet these tests are absent from symbolic regression (SR) benchmarking. SRBench [5] measures whether discovery methods find the right expression across 252 problems; SRSD [6] redesigns Feynman-derived tasks with realistic sampling ranges and discovery-oriented metrics. In neither framework does anything measure whether an attached falsification layer would reject a wrong form. The asymmetry has practical consequences: a discovery pipeline with no working falsification layer silently accepts misfit forms, and its downstream use — extrapolation, engineering design, scientific claim-making — inherits the error silently.

Measuring the falsification side faithfully requires resolving a protocol problem that, to our knowledge, no prior work addresses for this test family: the tests reject on incommensurable bases. A bootstrap-t threshold ("alarm when |t| > 3") and an asymptotic F p-value ("reject when p < 0.05") measure different things; comparing them directly measures threshold placement, not statistical power. Our first pilot did exactly this and produced inverted rankings that vanished under a common basis. The protocol below (Section 3.1) therefore rests on three commitments: one significance level for all tests; one permutation null shared by all tests; and power compared only at a common empirical false-positive rate estimated on true-monomial tasks.

One question the benchmark's own results force into the open: every test evaluated so far answers a single-form question — "is this a power law?" — and a rejection terminates the inquiry. What should the pipeline do next? Section 3.4 introduces a constructive answer evaluated as a fourth suite: a family-viability battery that runs the same curvature statistic in the linearizing coordinates of four families (power, exponential, logarithmic, rational) and returns the viability set — the families not rejected — instead of a single verdict. The battery inherits every protocol commitment of the main benchmark, and adds one finding of its own: the validity of a null distribution is a coordinate-dependent property, and a null that is exact in one family's coordinates can silently break in another's.

The question this benchmark answers is narrow and falsifiable: **when the ground truth is known and the candidate form is wrong, which specification tests detect the misspecification — and at what false-alarm cost, on data distributions that physical model discovery pipelines actually consume?** We contribute: (i) the benchmark itself — three suites, 394,260 fitted datasets plus 472,500 held-out evaluations, released generators with seed rules that make any single cell independently reproducible; (ii) the protocol — matched-basis comparison with empirical-FPR calibration and explicit failure registration; (iii) the resulting failure-mode map — including two previously undocumented findings: the flat-in-n 42–56% false-positive rates of HC3-Wald RESET variants on real physical monomials, and the near-blindness of fitted-F RESET to rational laws; and (iv) the family-viability battery — 94,000 further datasets showing the rejection machinery extends to routing between families once — and only once — the null is rebuilt per coordinates.

## 2 Related Work

**Specification testing.** RESET [1] augments the regression with powers of the fitted values and tests the joint significance of the added coefficients; the exog variant uses powers of the regressors instead. Heteroskedasticity-consistent covariance estimation [2] and its HC3 leverage-weighted refinement [3, 4] are the standard small-sample corrections and are the default in several mainstream econometrics and statistics packages. These results are derived for linear models under fixed or quasi-fixed designs; none of the classical treatments contemplates log-transformed physical laws whose regressors span only a few orders of magnitude — precisely the regime our benchmark shows breaks HC3's small-sample behavior (Section 4.1).

**Resampling nulls under non-i.i.d. errors.** The wild bootstrap [10, 11] was developed for inference under heteroskedasticity, and its permutation cousin requires exchangeability. The battery's mechanism finding (Section 4.6) is a compact new instance of this classical dichotomy: a log transform homoskedasticizes multiplicative noise, so residual permutation is valid in coordinates ending in log y and invalid in coordinates ending in y or 1/y — the null's validity inherits the noise geometry of the coordinates the test runs in.

**Symbolic regression benchmarks.** SRBench [5] established the discovery-side evaluation: 14 SR methods and 7 ML baselines on 252 problems, scored by exact recovery and held-out accuracy. SRSD [6] reconceived the task distribution itself — realistic per-variable ranges reviewed against the physics, rather than toy domains — and is the source of our task configurations. Both measure the discovery side; a pipeline that returns a wrong expression and one that flags it score identically. Our generators deliberately reuse SRSD's official sampling objects so that the falsification-side pools are the same distributions the discovery side already standardizes on.

**Physics-informed and dimensionally-aware discovery.** AI Feynman [7] and data-driven dimensional analysis [8] use physical units to constrain the hypothesis space; PySR [9] represents the modern search side. Dimensional correctness is necessary but not sufficient for form correctness: every ground-truth law in our suites is dimensionally self-consistent, yet each non-monomial law is a form misspecification for a power-law model. The falsification side complements these constraints — units narrow what to search; specification tests check what was found.

**Honest-failure evaluation.** Existing SR benchmarks report no analogue of a false-positive rate because they contain no explicit rejection mechanism to calibrate. Work on selective prediction and abstention addresses related epistemics for classification; we are not aware of a symbolic regression benchmark that measures attached falsification behavior, per family, with preregistered labels and failure registration. This is the gap the present benchmark occupies.

## 3 Benchmark Design

### 3.1 Design principles

**P1 — Common rejection basis.** All four tests decide at the same significance level α = 0.05 on p-values. The permutation curvature test's empirical p-value is computed from B = 999 permutations of the residual-permutation null that its native bootstrap already uses; no test is judged by a heuristic threshold.

**P2 — Power only at matched empirical FPR.** Headline power comparisons calibrate each method's decision threshold to a 5% empirical false-positive rate, estimated on pooled true-monomial tasks within the same suite, sample-size stratum, and noise regime. Raw-α results are reported alongside, because raw tables are what current practice would see.

**P3 — Preregistered labels.** H0 (strict monomial: the log-log mean function is exactly linear) versus H1 (any non-monomial structure) is fixed from symbolic structure of the ground truth before any test runs. No label is revised after seeing results; the family stratification of S3, the battery's synthetic truth registry, and the battery's SRSD selection are likewise fixed from expression structure before any battery run.

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

### 3.4 Suite S4 — the family-viability battery

The four tests of §3.2 share one property: each evaluates a *single* null form. A rejection is actionable only as an escalation, but the family-stratified results of S3 hint the machinery could do more: each test is a curvature check *in some linearizing coordinates*, and several physical families linearize in known coordinates.

**Families and members.** The battery registers four families with their linearizing coordinate pairs (u(x), v(y)), u applied column-wise: power (log x, log y); exponential (x, log y); logarithmic (log x, y); rational (1/x, 1/y). For each family F, member test H0,F is "the conditional mean is linear in F's coordinates": the design is [1, U, q] where U are the standardized transformed features and q = meanⱼ(Uⱼ²) is the single quadratic column; the statistic is q's coefficient. Two construction differences from method A are deliberate: columns are standardized before squaring so that features with different log-variances (physical constants at 10⁻³⁴ next to meter-scale quantities) contribute equally, and the null is a wild bootstrap rather than a permutation — for the reason the next paragraph documents.

**Null validity is a coordinate-dependent property.** A residual-permutation null is exact only under exchangeable errors in the transformed coordinates. Log transforms homoskedasticize multiplicative noise — which is why the log-log design was chosen for the main benchmark — but the same noise is heteroskedastic in (x, log y), (log x, y), and (1/x, 1/y) coordinates. This prediction was preregistered in the probe script's docstring before any run: permutation should be valid for families ending in v = log y (power, exponential) and break for families ending in v = y or v = 1/y (logarithmic, rational). The probe measured exactly that (Table 9): under multiplicative noise, the permutation-null rational member rejects its *own* true family at 30–34% — false escalations caused purely by null misspecification — while power and exponential members stay calibrated. Replacing the permutation draws with a wild bootstrap (Rademacher signs on restricted residuals [10, 11]; B = 999, vectorized through the pseudo-inverse row, so the coefficient is recomputed as a linear functional without refitting) drops the rational false escalation to 4.5–9% and leaves the log-coordinate members' size unchanged. One line of null construction, no change to statistic or design. The formal run uses the wild null throughout.

**Viability semantics.** The battery's output is the viability set V = {F : p_F ≥ α}, not a winner. Three errors are defined and reported separately, never merged: *false exclusion* (the true family absent from V, the limiting case being false escalation V = ∅) — costs routing the dataset away from a usable form; *false inclusion* (a wrong family present in V) — costs an inconsistent shortlist downstream; and for structures outside all four families (quadratic, threshold, cubic, interaction, …), *escalation is the correct outcome*, reported as such, never as a win. Because V is a non-rejection set, each family's exclusion error is a per-family property and needs no cross-family multiplicity correction; what multiplicity does buy — near-neighbour families co-occurring in V — is reported as a finding rather than suppressed.

**Grid (frozen before the run).** n = 200; rel ∈ {0.5%, 1%, 3%, 5%}; α = 0.05 per family; 500 seeds per cell; seed rule 900,000,000 + cell_index × 10,000,000 + seed × 97. Synthetic side: all 19 S1 scenarios, each with a truth-family label frozen in the code's registry before the run — power (the 3 structural H0 scenarios and the 5 noise-layer scenarios), exponential, logarithmic, rational, and `none` for the 8 structures no battery family can linearize. SRSD side: the 14 H1 equations of the preregistered S3 family-stratified selection (a stratified-representative subset of S2's 19 H1 tasks, drawn by a frozen per-family cap before any battery result) × {homoskedastic, heteroskedastic} noise. Total: 188 cells, 94,000 datasets. The battery is implemented at the experiment layer over the project's core code path without modifying it; a core-integration patch with tests is released but deliberately not applied, so the audited core-code hashes of the main benchmark remain valid.

## 4 Results

Table numbering is consecutive: Tables 1–2 report size, Tables 3–4 power, Tables 5–6 the family map, Table 7 noise robustness, Table 8 extrapolation alarms, Tables 9–11 the battery. All rates are rejection fractions; Wilson 95% intervals for headline rates appear in Appendix C. Seed counts per cell: 500 (S1, S2, S4), 30 (S3).

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

On synthetic Gaussian data, all four tests are calibrated (Table 1): every rate falls within or near its nominal band, and no method is distinguishable from another (Figure: `figures/full_500_calibration.png`). On real physical laws the picture splits (Table 2). A and B remain calibrated — rates .049–.053, intervals overlapping 0.05, flat in n. Both HC3 variants break, and break *identically at every sample size*: C stabilizes at .417–.422, D at .563–.565, with 95% intervals that do not move from n = 50 to n = 1000. This is not a small-sample artifact that more data cures; it is a structural miscalibration.

The mechanism is collinearity of the built columns. Physical constants and engineered quantities span few log-decades; on such designs ŷ² and ŷ³ are nearly linear functions of the fitted values themselves (for C), and the built power columns for D are nearly constant or near-collinear with the design — S3 records the limiting case, a ⅓-power monomial whose built columns are effectively constant and where D's false-positive rate reaches 1.00. The sandwich estimator's asymptotic guarantee presumes the design matrix is well-behaved as n grows; here the design degeneracy is exact and does not wash out, so the Wald statistic's reference distribution is simply wrong. H0 exact-power-law cells of S3 reproduce the same signature at independent scale: C .52, D .36.

Noise-layer scenarios (Table 7, Section 4.4) confirm what does *not* cause this: heavy tails, growing variance, outliers, and autocorrelation leave all methods' sizes at .03–.07. Violations of the noise assumption are not what breaks specification tests on real data — the design is.

### 4.2 Power: raw, then calibrated

Raw power at α = 0.05 on S2 H1 (n = 200, pooled over rel ≥ 0.5%): A .52, B .28, C .63, D .94. A naive reading crowns D. Table 2 already priced that crown: D pays a .563 false-positive rate for it. Calibrating every method to a common 5% empirical FPR on the H0 pool (per sample size) reorders the field:

**Table 3 — S2 power at matched empirical FPR 5%.**

| Method | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | **.41** | **.51** | .59 |
| B | .19 | .28 | .42 |
| C | .05 | .28 | .41 |
| D | .19 | .37 | **.65** |

At deployment-relevant sample sizes (n ≤ 200), A delivers roughly twice B's calibrated power and 2–10× C's and D's. At n = 1000, D overtakes (.65): with enough data, exog powers resolve the very structures that the mean-square average compresses, and its calibration failure becomes less damaging relative to its information gain. Even so, no method exceeds .65 — most of S2's H1 mass is *hard* for every test, which Section 4.3 decomposes by family (Figures: `figures/full_500_power_heat.png`, `full_500_n_scaling.png`).

On S1 structural H1 (raw α, rel = 1%, n = 200), the same pattern appears in sharper form. Quadratic, threshold, exponential-shape, logarithmic, rational-single, omitted-variable, and collinear-quadratic scenarios are caught by all methods at ≥ .96. The informative cells are the blind spots:

**Table 4 — S1 blind-spot cartography (raw α, n = 200, rel = 1%).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| cubic (standardized odd term) | .40 | 1.00 | 1.00 | 1.00 |
| interaction (lx₁·lx₂ cross term) | .21 | 1.00 | 1.00 | .10 |
| opposing curvature (+c, −c, exchangeable) | **.02** | .98 | .98 | 1.00 |
| dilution (one curved variable of two) | 1.00 | .88 | .87 | 1.00 |

Three mechanisms, each verified by construction. (i) A's statistic is driven by the mean of squared centered columns; when two variables carry equal and opposite curvature on an exchangeable design (OPPOSING), that mean has exactly zero population signal — A's .02 at n = 200 and .00 at n = 1000 is the statistic doing precisely what it is defined to do. (ii) A quadratic-projection statistic loses power on odd-symmetric deviations: A's cubic power (.40–.47 across n) never approaches B's saturation, because the cubic's second-order projection is small relative to its third-order content. (iii) B's fitted-value powers are a local Taylor surrogate for *any* smooth alternative — which is why B saturates on interaction (1.00) where A (.21) and D (.10) struggle — and, as S3 will show, the same surrogate property near-blinds B on rational laws.

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

**Table 6 — S3 calibrated power by family.**

| Family | A | B | C | D |
|---|---|---|---|---|
| exponential | .11 | .05 | **.63** | .04 |
| rational | **.66** | .14 | .50 | .66 |
| trigonometric | **.70** | **.73** | .08 | .05 |

Read jointly, Tables 5 and 6 are the benchmark's central finding. B — the classical test, near-perfect on synthetic quadratics — is nearly blind to rational laws (.13 raw, .14 calibrated): the fitted-value Taylor surrogate absorbs precisely the low-order curvature by which a saturating law departs from a power law. D's raw dominance on rational (.95) and trigonometric (1.00) is false-alarm budget, not power: calibrated, its trigonometric advantage collapses to .05. C is the exponential family's only reliable detector (.63 calibrated) while being worthless on trigonometric (.08). A misses exponentials (.11) — an exponential law's log-log curvature is modest at official scales — but is otherwise the most even performer. Every method's strength sits on another method's blind spot; the four columns of Table 6 have no row-major winner.

**Sensitivity.** Scaling every sampling domain ×2 moves exactly one family-method cell by more than 15 percentage points (rational × A: .66 → .84 — larger domains make curvature easier to see, the expected direction); ×0.5 moves none. The heteroskedastic replication of every S3 cell agrees with the homoskedastic one within .02 on every family-method cell. Sampling-range choice and variance model — the two design knobs a practitioner might suspect of driving the conclusions — do not drive them.

### 4.4 Noise-layer robustness

**Table 7 — S1 noise-layer H1 (strict-monomial mean structure) rejection rates at α = 0.05; these are false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds).**

| Scenario | A | B | C | D |
|---|---|---|---|---|
| Student-t(3) noise | .034 | .036 | .030 | .030 |
| variance growing with prediction | .048 | .058 | .050 | .050 |
| AR(1) errors, ρ = 0.6 | .038 | .040 | .042 | .042 |
| 2% point outliers (+8σ) | .040 | .066 | .036 | .036 |
| 6% group outliers (+5σ) | .038 | .060 | .054 | .054 |

Every method, every non-Gaussian regime: .030–.066. Heavy tails, heteroskedasticity, serial correlation, and outliers — the standard rogues' gallery of regression diagnostics — do not push any of these tests above 7% false alarms. This table's function is exclusionary: it rules out the convenient explanation that the miscalibrations of Table 2 come from noisy, ugly, real-world data. The data in Table 2 are clean by construction; what differs from Table 1 is the *design* — physical constants and few-decade quantities. The design, not the noise, is what specification tests are fragile to.

### 4.5 Do alarms predict extrapolation failure?

Each S2 fit (H0 and H1 tasks alike — the deployment situation, since a pipeline does not know which candidates are right-formed) is scored on 2,000 same-distribution test samples and on a ×3-domain far-field probe against the closed-form truth. Base rate of far-field relative error above 10%: .198. Specification-test alarms enrich this substantially:

**Table 8 — Alarm precision and recall for far-field error > 10% (S2, all tasks, all n, rel ≥ 0.5%).**

| Method | precision | recall |
|---|---|---|
| A | .39 | .63 |
| B | .56 | .49 |

Two properties matter more than the margins. First, enrichment: an A alarm raises the probability of severe extrapolation failure from .198 to .39 (2.0×); a B alarm to .56 (2.8×). Second, complementarity: of the true-positive catches, 10,560 are A-only and 5,447 B-only, against 12,575 shared — the two tests' alarm sets overlap barely more than chance, so a portfolio of both captures substantially more failures than either. Alarms are therefore useful as *ordinal evidence* about form.

They are not certificates of magnitude. At the task level, the rank correlation between median in-domain p-value and median far-field error is moderate for A (Spearman ρ = −0.33) and stronger for B (ρ = −0.68; negative = smaller p, stronger alarm, goes with larger error), yet both orderings fail exactly where it hurts most: the suite's two worst extrapolators — i.39.11 and i.43.43, both at 98% median far-field error — sit at the weak end of every method's evidence (A rejection .06, B .27–.29 at n = 200, rel = 1%), because their in-domain curvature is near-zero while the failure lives beyond the domain edge. Conversely, ii.38.14 is rejected by A essentially always (rejection 1.00) yet extrapolates benignly (median far-field error 7.7%): an alarm means form evidence exists; a non-alarm means nothing about the field beyond the training domain, and an alarm does not price the error either. No rejection is evidence of extrapolation safety (Figure: `figures/full_500_extrap.png`).

### 4.6 The battery: from rejection to routing

**Null legality, probe evidence (Table 9).** The preregistered prediction — permutation valid only in coordinates ending in log y — was measured at 200 seeds per cell before the formal run:

**Table 9 — Probe: permutation vs wild null by family (n = 200, 200 seeds; entries rel = 1% / 3%).**

| Truth | perm null | wild null |
|---|---|---|
| power (size of power member) | .070 / .075 | .075 / .060 |
| rational: false escalation (V = ∅) | **.34 / .31** | .045 / .090 |
| rational: containment (rational ∈ V) | .66 / .70 | **.96 / .91** |
| exponential: containment | .94 / .95 | .92 / .95 |
| logarithmic: containment | .87 / .93 | .89 / .95 |
| threshold (no battery family): escalation — *correct outcome* | 1.00 / 1.00 | 1.00 / 1.00 |

The rational row is the mechanism in one line: with the permutation null, the battery's own rational member rejects the best case for the rational test — data generated as y = x/(0.2 + 0.1x) with multiplicative noise — a third of the time. Nothing about rational laws is at fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The wild bootstrap repairs exactly that and nothing else: log-coordinate members (power, exponential) are untouched, confirming the coordinate-dependence reading rather than a blanket "wild beats permutation."

**Formal run, synthetic side (Tables 10–11).** At 500 seeds per cell, n = 200, wild null throughout, zero data-layer failures in 94,000 datasets:

**Table 10 — Battery on synthetic truths (pooled over rel ∈ {0.5–5%}; 2,000 datasets per truth family per cell).**

| Truth family | Containment (truth ∈ V) | False escalation (V = ∅) |
|---|---|---|
| power — size, 8 scenario groups incl. all noise layers | — | .0496 [.046–.053] pooled; per-scenario .043–.054 |
| exponential | .95 | .05 |
| logarithmic | .94 | .005 |
| rational | .93 | .066 pooled (≤ .082 worst rel) |

Under logarithmic truth the modal viable set is {logarithmic, rational} — the rational member co-survives in 91% of logarithmic-truth datasets. This is false inclusion in the strict sense, and it is the multiplicity property §3.4 promised to report rather than suppress: over a ×100 domain, a logarithmic and a saturating rational law are near-neighbours the battery cannot and should not arbitrate between. Arbitration belongs to the discovery layer (compare fitted candidates); the battery's contract is "never silently route a truth out of V," and after the wild fix that contract holds: containment .93–.95 against a 5% per-family budget, false escalation at or below 8.2% in the worst cell.

**Table 11 — Battery on out-of-family structures (truth ∈ no battery family; escalation is the *correct* outcome; 2,000 datasets per scenario).**

| Scenario | Correct escalation | False-inclusion rate | Modal false V |
|---|---|---|---|
| quadratic | 1.00 | .00 | — |
| threshold | 1.00 | .00 | — |
| omitted variable | 1.00 | .00 | — |
| dilution (1 of 2 features curved) | .92 | .08 | rational |
| cubic (odd term) | .24 | .76 | power (.68) |
| interaction (cross term) | .05 | .95 | power (.95) |
| opposing curvature | **.00** | 1.00 | power+rational (.84) |
| collinear-quadratic | .03 | .97 | rational (.61) |

Table 11 is the blind-spot cartography of Table 4 re-expressed as routing errors, and the correspondence is exact: the power member inherits A's documented blind spots (opposing curvature, odd terms, interactions — each orthogonal or aliased to the mean-square quadratic column), so those structures stay in V under the power label instead of escalating. Two rows are new information. Collinear-quadratic — where A alone among the four tests had full power — fails to escalate because the curvature hides behind the degenerate design for every member simultaneously. And dilution, A's synthetic win, escalates correctly 92% of the time: the battery converts A's best-case power into a correct routing decision.

**Formal run, SRSD side.** On the 14 preregistered H1 equations (homo- and heteroskedastic noise agreeing within .02 on escalation for 12 of 14; worst deviation .09 on i.48.2), the battery's routing is descriptive — no containment score is computed, for a reason worth stating. The selection table's family labels are surface-form classifications, and surface form is not linearizable form: the two equations labelled "exponential" that are Planck-type x³/(e^{c/x} − 1) laws are *correctly* rejected by the exponential member (rejection 1.00) — exp of a reciprocal argument is not y = a·e^{bx} in any coordinate the battery owns. Scoring containment against surface labels would mis-grade the battery for doing its job. What the descriptive record shows: median |V| is 1–2 for 11 of 14 equations (the exceptions being the escalator ii.38.14 at 0 and the under-informative i.26.2 and ii.11.28 at 3–4); i.30.5 routes to the singleton {exponential} in 78% of datasets; ii.38.14 — the main benchmark's most-rejected benign extrapolator (A rejection 1.00, far-field error 7.7%) — escalates on 67%, the battery correctly refusing to house it anywhere; and ii.11.28 keeps all four families viable over its official domain (per-member rejection ≈ .05), a near-power-law regime where no in-domain evidence exists to route on. A symbolically exact linearizability registry — checking candidate families against the closed form rather than a surface label — is the missing piece for fair SRSD-side scoring, and is future work.

### 4.7 Runtime

Per-dataset wall time (S1, all scenarios, means over 500-seed cells): A 12.0 ms (n = 50) to 25.3 ms (n = 1000) — the cost of 999 permutation evaluations, vectorized; B/C/D 1.0–1.8 ms across all n (Figure: `figures/full_500_timing.png`). The full four-family battery at n = 200 runs in 24 ms per dataset (94,000 datasets in 37 minutes single-threaded). All four tests and the whole battery are effectively free at the scale of a discovery pipeline that already fits thousands of candidates; the bottleneck resource is not compute but the *decision discipline* — the calibration set of Section P2 requires only true-monomial tasks that a physics repository already contains.

## 5 Analysis

**Q1. Does the low-df curvature test beat RESET on univariate quadratic deviations?** No — it ties. On S1's quadratic truth both saturate at power 1.00 from n = 50. What low degrees of freedom buy is calibration, not power: A stays within ±1 point of nominal size on every H0 pool in this paper, which is what makes its calibrated power numbers *trustworthy* rather than large. The asymmetries that do appear are structural: odd-symmetric cubic deviations (A .40–.47 vs B 1.00) and interaction terms (A .21 vs B 1.00, D .10), where A's quadratic projection is mathematically the wrong instrument.

**Q2. Does averaging squared columns miss feature-specific curvature and interactions?** Yes, by construction — and the benchmark quantifies when that matters. Cancellation is real: opposing curvature drives A to .02/.00. Dilution is not: one curved variable among two still yields A 1.00 versus B .88. Interactions blind A (.21) but saturate B. The decisive observation is external validity: neither the cancellation nor the interaction structure appears among the 32 real physical laws, where A's worst family is exponential (.11) and its blind spots are empty. Synthetic blind spots are true properties of the statistic; only a real-equation pool can reveal whether they bind in practice. Both suites are needed for that sentence to be sayable — and the battery's Table 11 adds the corollary: a synthetic blind spot stops being a power deficit and becomes a *routing* error (wrong family kept viable) the moment the same statistic is used to route.

**Q3. Does HC3 fix false alarms under heteroskedasticity?** Heteroskedasticity was never the disease here. Synthetic growing-variance noise leaves all methods at .048–.058 (Table 7), and every S3 heteroskedastic cell replicates its homoskedastic counterpart within .02. What breaks on real physical monomials is HC3 itself (Table 2): the robustness correction presumes a design whose leverage structure is informative; on near-collinear built columns the leverage values are uninformative, and the correction multiplies noise into a statistic whose reference distribution is wrong at every n measured. The robust estimator is the failure mode — an inversion of its textbook reputation, documented at 26,000 fits per cell.

**Q4. Do alarms predict extrapolation error?** As binary evidence, usefully: 2.0–2.8× enrichment of severe-error probability, with complementary alarm sets (Section 4.5). As magnitudes or guarantees, no: task-level rank correlations are ρ = −0.33 (A) and −0.68 (B), the two worst extrapolators in the suite carry the least in-domain curvature and are missed by every method, and the most-rejected benign extrapolator shows the converse. The correct reading is asymmetric: alarm ⟹ form evidence (and elevated extrapolation risk); silence ⟹ nothing. Pipelines that treat a passed form check as an extrapolation certificate are making an error this benchmark can now price: in S2, that error rate is the .198 base rate itself.

**Q5. Can the same machinery answer "not a power law — then what?"** Yes, as a viability set — after one repair the main benchmark could not have predicted. The battery is not a fifth test and claims no new statistic: each member is A's construction transplanted into a family's linearizing coordinates (with column standardization so unequal feature scales weight equally). Its first finding is about nulls, not families: the permutation null that calibrates A in log-log space is invalid in coordinates ending in y or 1/y, and the resulting 30–34% false-escalation rate on rational truths (probe, Table 9) is a *coordinate* property — log transforms homoskedasticize multiplicative noise, other transforms do not. The wild bootstrap, chosen precisely because it is valid under heteroskedastic errors of unknown form [10, 11], repairs the rational member to ≤ 9% and the formal run then delivers the contract: containment .93–.95 on all three battery families, pooled power-family size .0496 flat across structural and noise-layer scenarios, correct escalation at 1.00 on quadratic, threshold, and omitted-variable structures. The costs are three, each located: near-neighbour overlap (logarithmic ↔ rational, 91% co-viability) is inherent to non-arbitration semantics; the carry-over blind spots of Table 11 are A's blind spots wearing routing clothes; and on real equations the surface-form family labels dissolve under examination — the battery correctly rejects Planck-type "exponentials" that no linearizable family owns, which is why SRSD-side results are reported descriptively. The battery thus extends the falsification side one step up the pipeline — from "reject this form" to "hand the discovery layer a defensible shortlist" — without claiming the discovery layer's job: ranking within V is model selection, explicitly out of scope.

## 6 Limitations

**Single candidate family (main benchmark).** All four tests of §3.2 evaluate a power law against alternatives; nothing there benchmarks model *selection* among families, only rejection of one. The battery addresses routing, not selection: it returns V and does not rank its members; selection between viable families (e.g., by information criteria on family-fitted candidates) is the discovery layer's job and is not benchmarked. **Fixed augmentation.** Built-regressor powers are fixed at (2, 3); adaptive or higher-order constructions may shrink B's rational blind spot at calibration cost we have not measured. **Battery scope.** The battery runs at n = 200 only, single-decision per dataset; its synthetic truths are the S1 generators (single- and two-feature instances); it has no member for interaction or trigonometric structures, and the carry-over blind spots of Table 11 are quantified, not fixed. Its SRSD side is descriptive pending a symbolically exact linearizability registry; the collinear-quadratic false-inclusion row (.97 under rational) is registered as an open failure mode of member behavior on degenerate designs. **Log-space scope.** All main-benchmark tests operate in log-log space by design; additive-error models are out of scope and marked not-applicable rather than log-shifted, which defines the benchmark's boundary rather than a finding. **Dimensional manipulation absent.** Units were recorded but not experimentally manipulated; the observation that a degenerate ⅓-power monomial drives D's false alarms to 1.00 suggests structure-aware screening could pre-empty some failures, but this is untested here. **Statistical precision.** Headline S1/S2 rates rest on 500 seeds (Wilson intervals ±.003–.006 at p ≈ .05, ±.02 near p ≈ .5); S3 family cells rest on 30 seeds and support directional claims only; battery probe comparisons (Table 9) rest on 200 seeds per cell and carry the mechanism claim, with the formal run (500 seeds) carrying the headline rates.

## 7 Reproducibility statement

Every generator is seeded from a documented rule embedding the grid position (suite, scenario, n, noise level, seed index; the battery's rule is given in §3.4), so any single cell of any table regenerates bit-identically. The harness records content hashes of the code under test — including the exact `form_check` implementation — and of the harness itself, in every run's environment file; a statistic-mismatch sentinel verifies on every call that the wrapper's recomputation agrees with the native implementation, so drift in the tested code cannot pass silently. The battery is implemented at the experiment layer over the unmodified core code path, so the audited core-code hashes remain valid; a core-integration patch with its own test suite (153 tests green in a verification worktree) is released separately and deliberately not applied. All cells' raw results (parquet and CSV), per-equation tables, failure registrations with typed reasons, selection and exclusion tables with per-equation justifications, analysis scripts, and figure sources are released. Zero-noise grid cells deliberately exercise the failure-registration path: all four tests register degeneracy there rather than emitting p-values, and those registrations are part of the released record; the battery's 94,000 datasets contain zero data-layer failures, which is itself a registered outcome, not an omission.

## 8 Conclusion

We built the falsification-side benchmark that model discovery lacks: three suites, 394,260 fitted datasets and 472,500 held-out evaluations, one protocol that puts four classical specification tests on a common basis, and full registration of every failure; then extended it one pipeline step with a fourth suite — 94,000 datasets — that turns the same machinery from rejection into family routing. Its findings are specific and actionable. Robust-covariance RESET variants — the field's default "safe" choice — are structurally miscalibrated on real physical laws (42–56% false alarms, flat in n), an inversion of their textbook reputation documented here at a scale no single-dataset study could reach. Raw power rankings invert under matched false-positive rates, and invert again across task families: the test that dominates synthetic benchmarks is blind to rational laws; the test that is blind to interactions saturates on them synthetically. No test wins everywhere; the failure profiles are complementary; a portfolio of two tests at under 30 ms per dataset dominates any single one. The epistemic boundary is now priced: alarms enrich severe-extrapolation probability 2–3×, silence certifies nothing. And the battery shows the way past the dead end that a rejection is: route, don't just refuse — but rebuild the null per family first, because a null that is exact in one coordinate system is silently wrong in another. Benchmark papers are usually measured by the leaderboards they produce. This one is measured by the failures it makes visible — because in scientific discovery, the model you silently accept is the one that breaks your extrapolation.

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

## Appendix A — Suite composition

| Suite | Cells | Fits per cell | Total fits | Extra evaluations |
|---|---|---|---|---|
| S1 synthetic | 19 × 3n × 5σ = 285 | 500 | 142,500 | — |
| S2 SRSD-Feynman | 32 × 3n × 5σ = 480 | 500 | 240,000 | 240,000 same-distribution test sets; 232,500 extrapolation probes (7,500 not finite on probe domain, registered) |
| S3 family-stratified | 392 | 30 | 11,760 | — |
| S4 battery | 19 scenarios × 4σ = 76 synth; 14 eq × 4σ × 2 noise = 112 SRSD; 188 total | 500 | 94,000 | — |
| **Total** | **1,345** | — | **488,260** | **472,500** |

*Correction note (v1, retained): an earlier draft claimed "643,250 fitted datasets" and "1,167 cells" for the main benchmark; the audited totals are 394,260 fits over 1,157 cells, plus 472,500 held-out evaluations. This version supersedes the earlier figure everywhere; the annotated v0 draft retains it as an erratum example.*

## Appendix B — Failure registration taxonomy

| Reason | Meaning | Exercised by |
|---|---|---|
| `degenerate_zero_residual` | Zero-noise cell on exact ground truth; no residual variation → test undefined | S1/S2 σ = 0 cells; battery degenerate draws |
| `rank_deficient` | Design matrix numerically rank-deficient | collinear designs |
| `insufficient_df` | n − k < 1 after augmentation | n = 50 multi-feature cells |
| `hc3_leverage_degenerate` | 1 − hᵢ ≤ 10⁻¹²; HC3 undefined | degenerate built columns |
| `hc3_singular_covariance` | V_HC3 numerically singular | degenerate built columns |
| `statistic_mismatch` | Wrapper statistic ≠ native statistic (drift sentinel) | continuous audit |
| `checked_false_*` | Native test declines to run (non-positive data, n < 12, non-finite) | not-applicable tasks |
| `nonpositive` / `transform_failed` / `transform_nonfinite` / `degenerate_null` | Battery member precondition failures: data not positive in a family's coordinates; u/v transform overflow; transformed values non-finite; null draws degenerate | battery members on edge-case data |

All registered failures remain in the audit trail; rates' denominators are documented per table. Zero-noise cells intentionally trigger `degenerate_zero_residual` for all four methods — the registration path is tested by the benchmark itself. The battery's formal run registered zero data-layer failures across 94,000 datasets.

## Appendix C — Headline Wilson 95% intervals

S2 H0, n = 200 (26,000 fits per cell): A .053 [.050–.055]; B .050 [.047–.052]; C .417 [.411–.423]; D .563 [.557–.569]. S1 H0, n = 200 (4,500 fits per cell): A .037 [.032–.043]; B .041 [.036–.047]; C .045 [.039–.051]; D .042 [.036–.048]. Battery power-family size (S4, 16,000 member tests): .0496 [.046–.053]. Intervals for all other cells are computable from released per-cell counts.

## Appendix D — Change note v1 → v2

1. **Battery merged.** The extension previously maintained as `battery_extension_section.md` (probe-scale, 200 seeds) is now §3.4/§4.6/Q5, with all headline numbers refreshed from the formal 500-seed run (`runs/battery_full`): pooled power size .0496; containment .95/.94/.93 (exponential/logarithmic/rational); false escalation ≤ .082; out-of-family escalation and false-inclusion rows per Table 11; SRSD-side routing reported descriptively with the surface-label vs linearizability dissociation (new finding, replaces the probe draft's "SRSD validation in progress").
2. **Extrapolation section corrected.** v1's example "bonus.15: power 1.00, far-field error 9%" misattributed exog-HC3's raw rate (1.00, achieved while paying a 56% false-positive rate) to the benchmark's power measure; form_check's rate on bonus.15 is .18–.25. The example is replaced by ii.38.14 (A rejection 1.00, median far-field error 7.7%). The base rate .198 is now correctly described as pooled over H0+H1 tasks (deployment semantics; v1's text implied H1-only, for which the base rate would be .34). Rank correlations are refreshed from the 500-seed data: task-level Spearman −0.33 (A) / −0.68 (B); v1's "ρ = 0.25" was a pilot figure.
3. **Numbering fixed.** v1's §4 preamble misstated the table map (it announced "Table 7 the family map, Table 8 extrapolation" while the family map was Tables 5/5b and extrapolation was Table 7); v2 renumbers the family map as Tables 5–6, noise as Table 7, extrapolation as Table 8, battery as Tables 9–11.
4. **Totals updated.** 1,157 → 1,345 cells and 394,260 → 488,260 fitted datasets with the battery suite included; the v1 correction note about the main-benchmark totals is retained above.



