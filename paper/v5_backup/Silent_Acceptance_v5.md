# Silent Acceptance: A Falsification-Side Benchmark of Form-Specification Gates for Physical Model Discovery

**Working manuscript v5. Supersedes v4.** Revision scope: (i) the HC3-Wald results of v2–v4 were re-derived end to end after an independent-implementation check reattributed their headline failure to our own inference path, and every C/D number in this draft comes from a re-run on the stable path (§3.7, §5 Q3); (ii) one SRSD task was dropped from the H1 pool for a data-generation drift the sentinel caught (§3.3, App. B); (iii) the review-facing changes the v4 circulations asked for — task-level uncertainty on the main tables, a prespecified/registered vocabulary pass, explicit dataset-count definitions, an anonymous artifact, viability filtering instead of routing, and a reworked typesetting plan (cover and code listings moved to the supplement) — are folded in throughout. The pinv/Gram case study is now its own appendix (App. D), and the retraction it documents is part of the paper's record, not a footnote to it.

## Abstract

Symbolic regression and physical model discovery pipelines certify candidate models with goodness-of-fit metrics — R², cross-validation error, information criteria, coefficient standard errors — that measure how well a model fits, never whether an alternative model family would fit strictly better. A model can score near-perfectly on all of these metrics while having the wrong functional form, and fail the moment it extrapolates. Classical form-specification tests exist — Ramsey's RESET and its variants, with heteroskedasticity-robust inference — and so does the simpler workflow alternative: fit a nested rival and compare. No benchmark measures whether any of these gates work on the data distributions such pipelines consume: physical laws sampled over physically realistic ranges with multiplicative noise, evaluated in log-transformed space. None measures what happens on real observational data.

We introduce a falsification-side benchmark filling this gap. Six suites total 1,480 grid cells and 508,825 generated datasets (§3.8 defines the count; gate evaluations run higher, ≈2.1M), with 472,500 further held-out evaluations: (S1) 19 controlled synthetic scenarios; (S2) 32 tasks from the SRSD-Feynman registry under official sampling configurations; (S3) 30 Feynman equations stratified by structural family; (S4) a 188-cell family-viability battery; (S5) a gate comparison against a nested quadratic rival adjudicated by 5-fold CV or AICc; and (S6) six World Bank macro relations across four years. Every comparison runs on one significance level, power is compared only at matched empirical false-positive rate, labels were fixed before the runs, and every degenerate evaluation is registered with a typed reason.

The benchmark caught its own most striking claim and retracted it, and the retraction is the finding we lead with. Earlier drafts reported that heteroskedasticity-robust (HC3) RESET variants carry false-positive rates of 42–56% on true monomial physical laws, flat in n — an apparent inversion of the correction's textbook reputation. An independent SVD-based implementation (statsmodels), run on identical data, is calibrated (.067) on the very tasks where our path rejected at 1.00, while agreeing with it to three digits elsewhere. The failure was ours: fitting the HC3 sandwich through pinv(X′X) squares the design's condition number, and pinv's default truncation tolerance zeroes the smallest singular direction of the Gram matrix — cond(X) ≈ 10¹⁰·⁷ pushes cond(X′X) to ≈ 10²¹·⁸, σ_min(X′X) ≈ 8.5×10⁻¹⁰ sits seven orders below the ≈ 5×10⁻³ truncation floor, and the Wald statistic detonates. Column equilibration does not help; raising rcond helps only partially; factorizing X directly (SVD or QR) removes the problem entirely, at 0.36 ms per fit against 0.33 ms for the Gram path. On the corrected path, HC3 is calibrated on physical monomials (.053–.060 across n ∈ {50, 200, 1000}; task-level bootstrap [.050, .058]) — the "42–56%, flat in n" numbers in earlier drafts were properties of our harness, not of the correction (§4.1, App. D).

With the corrected numbers in place, the rest of the map holds and sharpens. No gate is uniformly best: rankings invert under matched false-positive rates and again across equation families, and each gate's strength sits on another's blind spot. The nested quadratic rival adjudicated by AICc remains the strongest default detector on both benches (.930/.940 synthetic, .525/.669 SRSD at calibrated FPR 5%), buying that performance with a naive 9–14% false-alarm rate until a reference pool calibrates it; calibration transfers across equations (leave-one-equation-out sizes .050–.054 on average, but .010–.095 for individual held-out equations — a spread we now show, not just the mean). The permutation curvature gate's remaining niche is real but narrower than v4 stated: it was empirically size-valid on every pool we measured without a reference set, an approximation whose finite-sample limits we do not claim to have bounded. On real data, four gates agree where signal is strong and diverge where the synthetic map predicts weakness, with an observational dilution signature now visible in both the permutation and the corrected HC3-fitted gate (§4.8). The corrected exog-HC3 gate, once its raw power is read at its true calibration, is the single strongest detector on the physics bench at every sample size (.494/.675/.733) — a rehabilitation the earlier drafts' artifact made impossible to see (§4.2). We release all generators, per-cell results, failure registrations, and content-hashed configurations under fixed seed rules, with the retraction chain scripted end to end (§7).

## 1 Introduction

A regression model can be simultaneously excellent on every certification metric a discovery pipeline computes and wrong in its functional form. A power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes k-fold cross-validation, wins its information-criterion comparison, and reports small standard errors — while extrapolating with −20% error one decade outside the training domain. This is not hypothetical: the failure is reproduced as a controlled scenario in our laboratory notebook and appears in this benchmark as scenario S1 (§3.3), where every specification gate detects it at power ≥ 0.99. Fit quality and form correctness are different axes, and current discovery practice measures only the first.

The statistics literature has owned the second axis for over fifty years. Ramsey's RESET [1] tests whether powers of the fitted values carry additional signal; built-regressor variants use powers of the regressors; White [2] and MacKinnon–White [3] provide heteroskedasticity-robust inference. And the workflow alternative that requires no new statistics at all — fit a nested rival family and let cross-validation or an information criterion arbitrate — has been standard applied practice for about as long. None of this is measured by symbolic regression benchmarking. SRBench [5] measures whether discovery methods find the right expression across 252 problems; SRSD [6] redesigns Feynman-derived tasks with realistic sampling ranges. In neither framework does anything measure whether an attached falsification layer would reject a wrong form. The asymmetry has practical consequences: a discovery pipeline with no working falsification layer silently accepts misfit forms, and its downstream use — extrapolation, engineering design, scientific claim-making — inherits the error silently.

Measuring the falsification side raises a protocol problem no prior work addresses for this test family, as far as we know: candidate gates reject on incommensurable bases. A bootstrap-t threshold ("alarm when |t| > 3"), an asymptotic F p-value ("reject when p < 0.05"), and a cross-validated RMSE difference ("escalate when the rival wins") measure different things; comparing them directly measures threshold placement, not detection. Our first pilot did exactly this and produced inverted rankings that vanished under a common basis. The protocol below (§3.1) rests on three commitments: one decision basis per comparison; power compared only at matched empirical false-positive rate estimated on true-monomial tasks; and explicit registration of every failure.

This paper asks its benchmark the three questions a reviewer would ask, in the order they would ask them. (i) Do the classical gates work on physical data distributions? — §§4.1–4.5: calibration, power, blind zones, extrapolation semantics. (ii) A rejection terminates inquiry — then what? — §4.6: a family-viability battery that filters instead of refusing, once its null is rebuilt per coordinates. (iii) Why not just fit a nested rival and compare? — §4.7: because until it is measured, that question has no answer; measured, it wins, and the benchmark says so. A real-data suite (§4.8) then asks the external-validity question: does the synthetic failure-mode map predict anything about live observational data? It does, in both directions.

This revision adds a fourth question that we did not plan and could not have planned: does the benchmark itself survive its own machinery? It nearly did not. Our harness's HC3 path silently failed on exactly the designs the benchmark exists to stress — near-collinear, physics-derived, few-decade — and manufactured a clean, flat-in-n, 42–56% false-positive story that two full drafts shipped as their headline. The mechanism (§5 Q3, App. D) is a known numerical hazard of forming X′X with pinv defaults meeting a squared condition number, and it was invisible to us until an independent implementation disagreed by a full order of magnitude. What caught it was not vigilance but structure: a drift sentinel on the statistic (which covers method A but, we now note, did not cover C/D — a gap this revision closes in the released harness), a task pool fixed before the runs that made the failing subset nameable, and a re-analysis script whose numbers had to come from somewhere. We think the honest framing is that the benchmark worked the way replications work — late, and only because someone ran one — and that "the correction itself is the failure mode," the sentence v2 and v3 printed in bold, should stand in the record as a measurement of how easy it is to ship a numerical artifact when it confirms your frame.

The question this benchmark answers is narrow and falsifiable: when the ground truth is known and the candidate form is wrong, which gates detect the misspecification — at what false-alarm cost, on what reference budget, and do the synthetic failure modes transfer to real data? We contribute: (i) the benchmark — six suites, 508,825 generated datasets plus 472,500 held-out evaluations, released generators with seed rules making any cell independently reproducible; (ii) the protocol — matched-basis comparison with empirical-FPR calibration and failure registration; (iii) the failure-mode map on the corrected path, including the fitted-F RESET near-blindness to rational laws, the measured dominance and calibration cost of the nested-rival gates, and the pinv/Gram hazard as a deployment warning with its own appendix; (iv) the family-viability battery — 94,000 datasets showing rejection extends to filtering once the null is rebuilt per coordinates; and (v) a real-data suite showing the map's mechanisms operating outside the laboratory. What we deliberately do not contribute is a winning detector; the benchmark demoted its own best test in v4 and this revision retracts its own headline, both by the benchmark's protocol.

## 2 Related Work

**Specification testing.** RESET [1] augments the regression with powers of the fitted values and tests the joint significance of the added coefficients; the exog variant uses powers of the regressors instead. Heteroskedasticity-consistent covariance estimation [2] and its HC3 leverage-weighted refinement [3, 4] are the standard small-sample corrections and are defaults in several mainstream econometrics and statistics packages. These results are derived for linear models under fixed or quasi-fixed designs; none of the classical treatments contemplates log-transformed physical laws whose regressors span a few orders of magnitude. We note for calibration that HC3's finite-sample behavior on such designs is now an open question this benchmark supports but does not settle: our corrected path shows calibration on the pools we measured (§4.1), and stress tests beyond those pools (leverage-adjusted variants, clustered dependence) remain unmeasured (§6).

**Model selection as a specification gate.** Comparing a null model against a nested rival by information criterion [12] or cross-validation [13] is the oldest workflow-level misspecification check in applied practice. It is so obvious that a reviewer of any new test will ask for it; it is so workflow-dependent (which rival? which criterion? which threshold?) that benchmarking it requires the protocol this paper introduces — a frozen rival family, frozen decision statistics, and thresholds calibrated on a reference pool. Suite S5 is, as far as we know, the first head-to-head between a resampling-null specification test and nested-rival selection gates under matched false-positive rates.

**The wider specification-testing literature.** Econometrics and statistics hold many more form-specification tools than the gates tested here: Neyman's smooth tests and their data-driven derivatives [21, 22], conditional-moment and generalized likelihood-ratio tests [23, 24], the rainbow test [25], kernel and nearest-neighbor lack-of-fit statistics [26], spline/GAM residual comparisons, cross-fitted and sample-splitting variants, and conformal-style selective-risk diagnostics [27]. S5's nested-quadratic rival is the discovery community's canonical objection, not a survey of this literature; LOOCV, BIC, held-out likelihood, or non-nested rivals may rank differently. The gates were chosen because they are what discovery pipelines actually attach today (RESET variants, robust inference) plus the workflow alternative a reviewer would name; the harness is gate-agnostic — a new gate is one statistic-plus-decision function over the same fitted cells — and extending the benchmark to this literature is future work, not a closed set.

**Resampling nulls under non-i.i.d. errors.** The wild bootstrap [10, 11] was developed for inference under heteroskedasticity; its permutation cousin requires exchangeability. The battery's mechanism finding (§4.6) is a compact instance of this classical dichotomy: a log transform homoskedasticizes multiplicative noise, so residual permutation is empirically calibrated in the log-log coordinates we measured and breaks in coordinates ending in y or 1/y. We phrase this as an empirical statement about measured pools, not a validity theorem: estimated OLS residuals are not exactly exchangeable under any linear model, and the finite-sample quality of the approximation under leverage or stronger heteroskedasticity is not bounded here (§3.2, §6).

**Symbolic regression benchmarks.** SRBench [5] established the discovery-side evaluation: 14 SR methods and 7 ML baselines on 252 problems, scored by exact recovery and held-out accuracy. SRSD [6] reconceived the task distribution itself — realistic per-variable ranges reviewed against the physics — and is the source of our task configurations. Both measure the discovery side; a pipeline that returns a wrong expression and one that flags it score identically. Our generators reuse SRSD's official sampling objects so that the falsification-side pools are the same distributions the discovery side already standardizes on.

**Physics-informed and dimensionally-aware discovery.** AI Feynman [7] uses physical units to constrain the hypothesis space; measurement-driven discovery distills governing equations and structure directly from data [8]; PySR [9] represents the modern search side. Dimensional correctness is necessary but not sufficient for form correctness: every ground-truth law in our suites is dimensionally self-consistent, yet each non-monomial law is a form misspecification for a power-law model. The falsification side complements these constraints — units narrow what to search; specification gates check what was found.

**Honest-failure evaluation.** Existing SR benchmarks report no analogue of a false-positive rate because they contain no explicit rejection mechanism to calibrate. As far as we are aware, no symbolic regression benchmark measures attached falsification behavior per family with labels fixed before the runs and failure registration — or one that benchmarks its own proposed test against the alternative its reviewers would name, and then retracts its own headline when an independent implementation disagrees. That last clause is the gap this revision occupies.

## 3 Benchmark Design

### 3.1 Design principles

**P1 — Common rejection basis.** Within any comparison, all gates decide on the same basis at the same significance level α = 0.05. The permutation curvature gate's empirical p-value comes from B = 999 permutations of the residual-permutation null; RESET variants use asymptotic F reference distributions; S5's selection gates use thresholded statistics, with thresholds at naive zero and at an H0-pool-calibrated value — naive and calibrated reported side by side, because raw tables are what current practice sees.

**P2 — Power only at matched empirical FPR.** Headline power comparisons calibrate each gate's decision threshold to a 5% empirical false-positive rate, estimated on pooled true-monomial tasks within the same suite, sample-size stratum, and noise regime. We report the per-gate thresholds (τ) next to the power they buy, and we flag a subtlety the v4 draft missed entirely: when a gate's H0 p-distribution is itself an artifact (as ours was), its "calibrated threshold" calibrates to the artifact, and the resulting power column inherits the pathology (§4.2, §5 Q3). Thresholds on the corrected path are re-derived from corrected pools.

**P3 — Labels and analyses specified before the run.** H0 (strict monomial: the log-log mean function is exactly linear) versus H1 (any non-monomial structure) was fixed from symbolic structure of the ground truth before any gate ran. The family stratification of S3, the battery's truth registry, the S5 gate definitions, and the S6 relation list were frozen in code before their runs; no label or configuration was revised after seeing results. We use "prespecified" for everything frozen in versioned code and docstrings with repository history as the timestamp record, and we do not claim external registry-level preregistration; where a claim depends on that distinction (the battery's coordinate-dependence prediction, §3.4), we point to the specific commit.

**P4 — Explicit failure registration.** Every evaluation that cannot produce a valid decision is registered with a typed reason — zero-residual degeneracy, design rank deficiency, insufficient residual degrees of freedom, leverage degeneracy, singular covariance, statistic mismatch, generation failure (S5 SRSD side), indicator archival (S6) — and excluded from rate numerators while remaining in the denominators' audit trail. The zero-noise grid cells exist to exercise this path. To the taxonomy we add one entry this revision earned: `gram_conditioning_hazard`, raised when cond(X_aug) exceeds 10⁸ on a path whose covariance factors through X′X — the machine-readable form of App. D.

**P5 — Official data configurations.** S2/S3 tasks reuse the official per-variable sampling objects of the source registry; tasks whose inputs or outputs are not strictly positive under the official configuration are marked not-applicable and excluded — we do not shift or rescale data to make them log-fittable, since that changes the hypothesis under test. S6 uses a public observational source (World Bank Open Data API [14]) with the relation list frozen before any gate ran; the United Nations Foundation publishes no dataset API, and the World Bank is the UN-system institution of record for open macro data.

### 3.2 Gates under test

All gates operate in log-log space on positive data: let ly = log y and X = [1, lx₁, …, lx_k] be the design of the candidate power-law model.

**A. Permutation curvature (form_check).** Adds a single mean-centered quadratic column q = mean_j((lx_j − lx̄_j)²) and takes the added coefficient c₂ as statistic; size evidence comes from B = 999 residual permutations (vectorized through the pseudo-inverse row), with a statistic-mismatch sentinel comparing wrapper and native values on every call. Validity note, stated as an empirical claim: permuting estimated residuals assumes they are exchangeable, which OLS residuals are not exactly; what this benchmark establishes is calibration of the construction on the H0 pools we measured (Tables 1–2). Log-transformed multiplicative noise and favorable designs plausibly help; §6 lists the stress tests (leverage-adjusted residuals, wild and block nulls, cluster dependence, fixed-design variants) that would bound where the approximation breaks. Two properties of the vectorization matter for App. D: the statistic and the null share one pinv row, so path-specific truncation cancels between them — method A stayed calibrated on the same designs where C/D detonated, and this symmetry is evidence for the mechanism, not luck.

**B. RESET fitted, classical F.** Adds ŷ² and ŷ³; classical F test of joint significance. The statistic involves only residual sums of squares; it is untouched by the numerical pathology below, which is why it survives as the sentinel channel for the re-run (§3.7).

**C. RESET fitted, HC3-Wald.** Same augmented design; HC3 sandwich Wald statistic on the added coefficients, referred to F(q, n − k). As shipped in v2–v4, the sandwich factored through pinv(X′X). This revision reruns C on the SVD path below; the two paths agree to ≥ 3 significant digits on well-conditioned designs and diverge catastrophically on near-collinear ones (App. D). Which HC3 number is "the" HC3 number is therefore a property of the solver, and we treat the SVD path as primary throughout.

**D. RESET exog, HC3-Wald.** Built regressors are powers of each log-regressor; HC3-Wald as in C, rerun on the SVD path.

**E. CV-quad (S5).** The "obvious alternative": rival = full quadratic expansion of the log features (linear terms, per-feature squares, pairwise cross-terms; minimum-norm least squares absorbs rank deficiency). Statistic D = RMSE_cv(power) − RMSE_cv(quad) over 5 deterministic folds; escalate if D > δ. Naive δ = 0; calibrated δ* chosen on the H0 pool to FPR = 5%.

**F. AICc-quad (S5).** Same rival; escalate if AICc(power) − AICc(quad) > m, k = column rank. Naive m = 0; calibrated m* as in E.

Configurations were frozen before all main runs — powers (2, 3); B = 999; α = 0.05; the rival family and criteria of E/F; the grids below — and are identical across suites. The same harness invokes all gates on every dataset, so comparisons share data, seeds, and failure handling exactly. The statistic-mismatch sentinel that covers A now covers C/D as well: each C/D call is cross-checked against a statsmodels reference on a rolling subsample, and a `gram_conditioning_hazard` registration fires when cond(X_aug) > 10⁸ (P4).

### 3.3 Task suites S1–S3

**S1 — Synthetic control (19 scenarios × 500 seeds).** Three monomial H0 scenarios (univariate; two-variable; two-variable with ρ ≈ 0.9 correlated design). Eleven structural H1 scenarios: quadratic curvature (c = 0.10); cubic (standardized odd term); multiplicative interaction (log-space cross term); threshold (exponent 0.4 → 1.0 at the domain midpoint); exponential law; logarithmic law; rational (saturating) law; omitted variable; dilution (one of two variables curved); opposing curvature (+0.12, −0.12 on an exchangeable design, so the mean-square column's population coefficient is exactly zero); collinear-quadratic (ρ ≈ 0.95 design with curvature). Five noise-layer H1 scenarios have strict-monomial mean structure but non-Gaussian noise — Student-t(3); SD growing with the log prediction (0.3× to 1.7×); 2% point outliers (+8σ); 6% group outliers (+5σ); AR(1) ρ = 0.6 — rejections there are false alarms by construction and are never pooled into structural power. Grid: n ∈ {50, 200, 1000}; relative Gaussian noise σ ∈ {0, 0.5%, 1%, 3%, 5%} in log space; 500 seeds per cell via a seed rule embedding the grid position.

**S2 — SRSD-Feynman tasks (32 tasks × 500 seeds; one H1 task withdrawn).** From the 120-equation registry, 32 tasks pass the positivity screen under official configurations: 13 strict-monomial H0 and 19 non-monomial H1. Each training fit (n ∈ {50, 200, 1000}) is accompanied by a 2,000-sample same-distribution test set and a far-field probe: the first regressor extrapolated ×3 beyond both training domain edges, scored against the closed-form truth in a numerically stable form. During the C/D re-run (§3.7), the sentinel flagged task feynman-i.30.5: on regenerated data its archived classical-F p-values no longer reproduce, for every seed and every n, while all 31 other tasks reproduce to ≥ 6 digits. The upstream repository, clone, and environment are identical, so the drift is confined to that task's data generation — a three-variable exponential-type expression — and we could not locate it; the task is withdrawn from this draft's H1 pool (19 → 18) and its 7,500 re-run rows are released with the sentinel flags rather than silently averaged in. The H0 pool (13 tasks) is untouched. Every number in this draft comes from the 31-task pool; App. B tabulates exactly which v4 numbers move and shows that the movement is fully explained by the pool change (32-task pools reproduce v4's A/B rows exactly).

**S3 — Family-stratified Feynman equations (30 equations × 30 seeds).** Structural family classification from symbolic form: exact power law (12 H0); fractional power (1 H0); exponential (3 H1); rational (6 H1); trigonometric (5 H1); plus not-applicable members with explicit reasons. Each H1 equation runs over relative noise {0.5%, 1%, 3%, 5%} × {homoskedastic, heteroskedastic} × sampling-range scale {official, ×0.5, ×2}; H0 equations run the noise grid at official range.

### 3.4 Suite S4 — the family-viability battery

Each gate of §3.2 evaluates a single null form; a rejection terminates inquiry. The battery registers four families with their linearizing coordinate pairs — power (log x, log y); exponential (x, log y); logarithmic (log x, y); rational (1/x, 1/y) — and for each family F runs A's construction in F's coordinates, returning the viability set V = {F : p_F ≥ α} rather than a verdict. Columns are standardized before squaring so features at different scales weight equally, and the null is a wild bootstrap (Rademacher signs; B = 999) rather than a permutation — for the reason §4.6 documents: permutation's empirical calibration is coordinate-bound in every pool we measured. The prediction — permutation holds where the v-coordinate is log y and breaks where it ends in y or 1/y — was written into the probe script's docstring, with repository history as the timestamp, before the probe ran. Three errors are reported separately and never merged: false exclusion (truth absent from V), false inclusion (wrong family present), and escalation on structures no family owns — the correct outcome, never scored as a win. Grid: n = 200; rel ∈ {0.5–5%}; 500 seeds; 188 cells; 94,000 datasets, all 19 S1 scenarios with truth-family labels frozen in code plus 14 prespecified SRSD H1 equations × {homoskedastic, heteroskedastic}. Because V is a union of per-family non-rejections, this is a viability *filter* — a candidate-set stage whose output the discovery layer must still rank — and we no longer call it routing; the battery's job is to keep the truth on the shortlist, not to pick the winner (§4.6).

### 3.5 Suite S5 — the gate comparison

S5 exists because of an objection this benchmark could not answer from its own data: why not skip the test statistic entirely and let a nested rival arbitrate? Gates E and F are that objection, made runnable. Design, frozen in gate_compare.py's docstring before the run: benches = S1's 19 scenarios and S2's 32 tasks (the re-run inherits S2's withdrawal for H1-power rows; H0 pools are unaffected); grid n ∈ {50, 200}, rel = 3%, 200 seeds per cell (seed bases 300,000 / 400,000 — disjoint from all prior suites); rival = full quadratic expansion; E and F thresholds naive and H0-pool-calibrated per bench × n; A run under its standard protocol on identical data. Metrics: H0 size (naive and calibrated), structural-H1 detection at calibrated FPR 5%, noise-layer fire rates, per-scenario and per-equation detail (never macro-averaged away), and cost in fit units. One asymmetry is inherent and reported rather than hidden: A's permutation null is self-calibrating, while E/F at naive thresholds inherit nested-comparison tie bias; the calibrated variants price what a reference pool buys. S5 uses only A/E/F and is therefore untouched by the C/D re-run; its numbers are byte-identical to v4's.

### 3.6 Suite S6 — real data (World Bank open macro data)

Synthetic and equation-derived suites have ground truth by construction; the deployment question is whether the synthetic failure-mode map predicts anything on observational data, where no truth exists. S6 is deliberately small, prespecified, and divergence-oriented: six macroeconomic relations (target, features) × four years (2019, 2021, 2022, 2023), cross-country cross-sections of 184–199 sovereign economies (World Bank aggregates excluded), from eight indicators pulled from the World Bank API with raw responses cached for reproducibility. Gates A–D run at 5 permutation seeds per cell. Those seeds are computational replicates of one p-value estimator on one fixed dataset — not five independent samples — so they support only a Monte Carlo stability statement (empirical p resolution 1/1000 at B = 999; median spread .009, max .094 across the 24 S6 cells), never sample-size or stability claims. A prespecified probe adds three agricultural and energy relations over three years, including a built-in control-contrast pair — cereal yield ~ fertilizer (one feature) versus cereal yield ~ fertilizer + energy use per capita (two features) — chosen to test, on live data, whether the dilution mechanism of S1 operates outside the laboratory. No calibration or power claim is made on S6; its contract is external validity of mechanisms, not of rates. The C/D re-run touches S6's C/D columns; because S6 designs sit at cond(X) ≈ 10^5.7 (max 10^7.2) — conditioned but not pathological — the correction changes no S6 verdict (§4.8).

### 3.7 The C/D re-run: what was redone, and how it is anchored

After the independent-implementation check (§5 Q3) reattributed the HC3 failure to the pinv(X′X) path, every C/D evaluation in S1–S3, S6, and the engineering bench was recomputed with hat and sandwich factors both taken from X's SVD (equivalently pinv(X)pinv(X)ᵀ; the Gram matrix is never formed). A and B were not recomputed: A's statistic and null share one pinv row (truncation cancels between them), and B's F statistic involves only residual sums of squares. Three anchors make the re-run auditable rather than trusting. (i) *Sentinel*: on every re-run row, the harness recomputes B from regenerated data and compares against the archived value at rel 10⁻⁶. Results: S1 142,500/142,500 pass; S3 11,760/11,760 pass (plus 401-row pre-check); S2 232,500/240,000 pass, with the 7,500 exceptions all in the withdrawn task (§3.3). (ii) *Unchanged where theory says unchanged*: on S1's wide-domain synthetic designs (cond(X) ~ 10¹) the re-run reproduces the archived C/D rates cell for cell — Tables 1, 4, and 7 are numerically identical to v4's, which localizes the artifact to narrow, physics-constant-driven designs. (iii) *Deterministic regeneration*: identical seed rules, identical upstream commit, identical environment (numpy 2.5.2, scipy 1.18.1, Python 3.13); the re-run scripts and per-row outputs are released (runs/cd_svd/). Runtime cost of the stable path is 0.36–0.43 ms per fit against 0.33 ms for the Gram path (n = 200, two features) — the fix is free.

### 3.8 What the dataset counts count

"508,825 generated datasets" counts dataset generations: one row per (suite cell × seed), each a freshly generated data table to which gates are applied. It does not count gate evaluations (a dataset typically receives four), regression fits (gates fit different models internally), or cells. Gate evaluations are ≈ 2.1M (S1–S4: 4 per dataset; S5: 3 gates + union; S6: 4 gates × 5 seeds; the engineering bench: 4 per cell). The engineering bench (App. E) sits outside the seeded totals because its cells are single realizations. Earlier drafts' appendix multiplied cells × fits-per-cell, which is the generation count; the body's "fitted datasets" was ambiguous, and this section is the definition.

## 4 Results

Table numbering is consecutive: Tables 1–2 size, Tables 3–4 power, Table 5 the family map, Table 6 noise, Table 7 extrapolation, Tables 8–10 the battery, Tables 11–12 the gate comparison, Tables 13–14 real data. All rates are rejection fractions; Wilson 95% intervals for headline rates appear in Appendix C. Seed counts per cell: 500 (S1, S2, S4), 30 (S3), 200 (S5), 5 (S6 — permutation-p stability, not rate estimation). C/D columns everywhere are the SVD path (§3.7); v4's Gram numbers appear only where the contrast is the point, marked "Gram" and confined to App. D.

### 4.1 Size: false-positive rates on true monomials

**Table 1 — S1 synthetic H0, false-positive rate at α = 0.05 (Gaussian noise pooled over σ ∈ {0.5%…5%}; 4,500 fits per cell).**

|Gate|n=50|n=200|n=1000|
|---|---|---|---|
|A — permutation curvature|.047 [.041–.053]|.037 [.032–.043]|.041 [.036–.047]|
|B — RESET fitted F|.045 [.039–.051]|.041 [.036–.047]|.038 [.033–.044]|
|C — RESET fitted HC3 (SVD)|.048 [.042–.054]|.045 [.039–.051]|.038 [.033–.044]|
|D — RESET exog HC3 (SVD)|.054 [.048–.061]|.042 [.036–.048]|.040 [.034–.046]|

**Table 2 — S2 true physical monomials (13 SRSD-Feynman equations; synthetic draws from real equations, not measurements), false-positive rate at α = 0.05 (26,000 fits per cell).**

|Gate|n=50|n=200|n=1000|
|---|---|---|---|
|A|.049 [.046–.052]|.053 [.050–.055]|.052 [.049–.054]|
|B|.050 [.047–.053]|.050 [.047–.052]|.053 [.050–.055]|
|C (SVD)|.057 [.054–.059]|.053 [.051–.056]|.055 [.052–.058]|
|D (SVD)|.060 [.058–.063]|.056 [.054–.059]|.058 [.055–.061]|

On synthetic Gaussian data all four gates are calibrated (Table 1), and on physical monomials all four now are too (Table 2): the HC3 variants sit at .053–.060, one to two points above A and B but inside or brushing their Wilson bands at every n. This is the corrected map. The v2–v4 drafts printed C at .417–.422 and D at .563–.565, flat in n, and read that as a structural failure of robust inference; App. D documents why that reading was an artifact of our inference path, and the granularity claim that made it vivid — five equations failing at ≈ 1.00 while eight sat at nominal — dissolves under the corrected path into a narrow band across all thirteen equations (C: .052–.057; D: .054–.062). Task-level bootstrap over the 13 equations (10,000 replicates; App. C) gives C .053 [.050, .056] and D .056 [.055, .058], with leave-one-equation-out ranges of [.052, .054] and [.056, .057]: what earlier drafts reported as a bimodal, equation-specific pathology is, on a stable solver, a homogeneous and boring null. We would rather ship this sentence than the exciting one it replaces.

What the episode establishes is a boundary, and it is worth stating precisely because it is easy to misread as an exoneration. The Gram-path failure appears only where designs are narrow and near-collinear — physical constants at 10⁻³⁴ next to meter-scale quantities give cond(X) ≈ 10⁷–10¹¹·⁵, and cond(X′X) then exceeds float64's usable range, so *any* solver that forms X′X and inverts it with default tolerances fails identically, pinv or explicit inverse; column equilibration does not rescue it, and disabling pinv's truncation recovers only part of the way (rejection .30 instead of 1.00 on the worst task — the remaining inflation is the ~6 lost digits of a squared-condition solve). Research code that factorizes the sandwich through X′X is common, the failure is total where it occurs, and nothing in the output — no warning, no rank deficiency, no non-finite value — distinguishes a 1.00 rejection rate from a calibrated one. The per-gate, per-suite, per-cell record plus the released re-run scripts (§3.7) are the remediation; App. D is the case study; §5 Q3 carries the corrected HC3 discussion.

### 4.2 Power: raw, then calibrated

Raw power at α = 0.05 on S2 H1 (n = 200, pooled over rel ≥ 0.5%, 31-task pool): A .489, B .253, C (SVD) .227, D (SVD) .661. Two readings are now possible where v4 could make only one. On the artifact-inflated numbers (C .62, D .93 raw), D's dominance was purchased partly with a .56 false-positive budget; that criticism stands against the Gram numbers but loses its target against the corrected ones. On the corrected path, D's raw edge over A is genuine detection carried by a roughly calibrated gate — and the calibrated comparison (Table 3) makes D the strongest single gate on this bench at every n, a result no earlier draft could state because its D numbers were path-dependent.

**Table 3 — S2 power at matched empirical FPR 5% (31-task pool; τ = H0-pool 5% quantile of each gate's own p-values).**

|Gate|n=50|n=200|n=1000|
|---|---|---|---|
|A|.386 (τ=.050)|.486 (τ=.046)|.570 (τ=.047)|
|B|.187 (τ=.051)|.254 (τ=.051)|.389 (τ=.048)|
|C (SVD)|.129 (τ=.042)|.240 (τ=.047)|.422 (τ=.047)|
|D (SVD)|.494 (τ=.039)|.675 (τ=.045)|.733 (τ=.043)|

Three qualifications, in the order a reader should weigh them. First, the A/B rows moved from v4 (.41/.51/.59 and .19/.28/.42) for a reason unrelated to the artifact: v4's pools included the withdrawn task, and App. B shows 32-task pools reproduce v4 exactly while 31-task pools give Table 3's values. Second, τ varies per gate (.039–.051) because thresholds are quantiles of finite, discrepant p-distributions; at B = 999 the permutation p has resolution 1/1000, which contributes a ±.005 fuzz to A's effective size. Third, and most important for what follows: D's calibrated lead depends on a reference pool — with no H0 tasks in hand, a practitioner cannot place τ_D, and at the naive threshold D's size on physical monomials was, in the Gram era, catastrophic; on the corrected path D's naive size is .06–.08 on our pools (mild over-alignment), which makes D usable without calibration here but says nothing about pools we have not measured. The deployment reading is deferred to §5 Q6 after S5.

Table 4 (S1 blind-spot cartography) is unchanged from v4 and from v3, and this is by design rather than inertia: S1's wide-domain designs never triggered the artifact (the re-run matches archived C/D rates cell for cell, §3.7), so the synthetic contrasts survive correction exactly.

**Table 4 — S1 blind-spot cartography (raw α, n = 200, rel = 1%; identical to v4).**

|Scenario|A|B|C (SVD)|D (SVD)|
|---|---|---|---|---|
|cubic (standardized odd term)|.40|1.00|1.00|1.00|
|interaction (lx₁·lx₂ cross term)|.21|1.00|1.00|.10|
|opposing curvature (+c, −c, exchangeable)|.02|.98|.98|1.00|
|dilution (one curved variable of two)|1.00|.88|.87|1.00|

Three mechanisms, each verified by construction: (i) A's statistic is the mean of squared centered columns — under opposing curvature on an exchangeable design that mean has exactly zero population signal; (ii) a quadratic-projection statistic loses power on odd-symmetric deviations; (iii) B's fitted-value powers are a local Taylor surrogate for any smooth alternative — saturating on interaction where A (.21) and D (.10) struggle — and, as S3 shows, the same surrogate property near-blinds B on rational laws.

### 4.3 Family-stratified results and sensitivity

Table 5 is where the re-run bites hardest, and where the corrected numbers change the paper's story rather than just its decimals. Calibrated to 5% empirical FPR on the S3 H0 pool:

**Table 5 — S3 calibrated power by family (C/D on the SVD path; A/B unchanged from v4).**

|Family|A|B|C (SVD)|D (SVD)|
|---|---|---|---|---|
|exponential|.11|.05|.07|.06|
|rational|.65|.16|.13|.67|
|trigonometric|.66|.69|.58|.86|

Compare v4's Table: there, C was the exponential family's only reliable detector (.63) and D was worthless on trigonometric (.05). Both claims were artifacts speaking. On stable solvers, C is nearly blind to exponentials (.07) — its v4 power was Gram noise wearing an interpretation — while D is the strongest trigonometric detector on the bench (.86), not the weakest. B's rational blindness (.16 calibrated, the one cell of the v4 story that survives contact with the re-run) remains real: the fitted-value Taylor surrogate absorbs exactly the low-order curvature by which a saturating law departs from a power law. The corrected family map's shape is different too: v4's table had no row-major winner and made complementarity the headline; the corrected table has a near-row-major winner on two of three families (D on rational and trigonometric, A competitive on exponential only in relative terms — every gate is weak on exponentials, whose log-log curvature at official scales is genuinely small). A is no longer "the most even performer" so much as "the most even mediocre one"; its case now rests on calibration without a reference pool (§5 Q6), not on family coverage.

Sensitivity: scaling every sampling domain ×2 moves exactly one family-gate cell by more than 15 points (rational × A: .65 → .84, the expected direction); ×0.5 moves none; heteroskedastic replications agree within .02 on every cell — both sensitivity statements re-verified on the corrected path.

The corrected H0 sizes on S3's official-range pool (1,680 fits per gate): A .036, B .050, C (SVD) .046, D (SVD) .044. The v4 draft's S3 H0 row — C .457, D .450, including a 1.00 for D on the fractional-power monomial — was the artifact at full amplitude; on the stable path that degenerate ⅓-power monomial is an ordinary null (.04–.05 for every gate).

### 4.4 Noise-layer robustness

**Table 6 — S1 noise-layer H1 (strict-monomial mean structure) rejection rates at α = 0.05 — false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds; identical to v4, re-verified).**

|Scenario|A|B|C (SVD)|D (SVD)|
|---|---|---|---|---|
|Student-t(3) noise|.034|.036|.030|.030|
|variance growing with prediction|.048|.058|.050|.050|
|AR(1) errors, ρ = 0.6|.038|.040|.042|.042|
|2% point outliers (+8σ)|.040|.066|.036|.036|
|6% group outliers (+5σ)|.038|.060|.054|.054|

Every gate, every non-Gaussian regime: .030–.066. The table's function is exclusionary: whatever breaks on narrow physical designs is not noise, and the re-run confirms the S1 noise layers behave identically under both solver paths.

### 4.5 Do alarms predict extrapolation failure?

On the 31-task pool (v4's Table 8 was 32-task; App. B tabulates the difference), the base rate of far-field median relative error above 10% is .143, and alarm precision/recall are A .38/.91 and B .53/.54. v4's numbers (.198 base; A .39/.63; B .56/.49) shift because the withdrawn task carried an unusual far-field profile — its removal cuts the base rate and mechanically lifts every gate's recall while lowering B's precision. The qualitative conclusions survive with corrected magnitude: alarms enrich severe-error probability (A 2.7×, B 3.7× on this pool), alarm sets overlap barely more than chance, and the ordering of task-level horror stories is unchanged — the suite's worst extrapolators still carry the least in-domain curvature and are missed by every gate, while ii.38.14 is still rejected by A essentially always yet extrapolates benignly. Alarm → form evidence; silence → nothing. No rejection is evidence of extrapolation safety. We flag the pool sensitivity itself as the finding a reader should take from this subsection's revision: a one-task change moves the base rate by a quarter of itself, which is a statement about how thin S2's tail really is.

### 4.6 The battery: from rejection to viability filtering

**Table 7 — Probe: permutation vs wild null by family (n = 200, 200 seeds; entries rel = 1% / 3%).**

|Truth|perm null|wild null|
|---|---|---|
|power (size of power member)|.070 / .075|.075 / .060|
|rational: false escalation (V = ∅)|.34 / .31|.045 / .090|
|rational: containment|.66 / .70|.96 / .91|
|exponential: containment|.94 / .95|.92 / .95|
|logarithmic: containment|.87 / .93|.89 / .95|
|threshold (no battery family): escalation — correct outcome|1.00 / 1.00|1.00 / 1.00|

With the permutation null, the battery's own rational member rejects its own true family a third of the time — nothing about rational laws is at fault; the permutation's empirical calibration does not transfer to (1/x, 1/y) coordinates under multiplicative noise. The wild bootstrap repairs exactly that and nothing else.

Set-level errors and multiplicity. The battery's decision is the viable set, so the errors that matter are set-level, and each member runs at α = .05 without correction — a union decision whose two error types move oppositely under multiplicity. On H0-side cells the reported containment is one minus the exclusion rate of the true family's own member (.0496–.066 across the four truth families; the viable set is never empty), and per-member α does not inflate that rate because exclusion depends on one member only. Inflation appears on the other side: some wrong member survives alongside the truth in 23.3% of power-truth cells and 92% of logarithmic-truth cells — near-neighbour overlap, reported rather than suppressed. On the structural-H1 pool, some family is retained in 100% of cells: the battery is a candidate filtering stage whose surviving set the discovery layer must still rank, not a final selection, and downstream selection error after the battery is out of scope and unmeasured. "Routing" was always too strong a word for a mechanism that never empties the set; we use filtering throughout this draft.

**Table 8 — Battery on synthetic truths (500 seeds/cell, wild null).**

|Truth family|Containment (truth ∈ V)|False escalation (V = ∅)|
|---|---|---|
|power — size, 8 scenario groups incl. all noise layers|—|.0496 [.046–.053] pooled; per-scenario .043–.054|
|exponential|.95|.05|
|logarithmic|.94|.005|
|rational|.93|.066 pooled (≤ .082 worst rel)|

**Table 9 — Battery on out-of-family structures (escalation is the correct outcome).**

|Scenario|Correct escalation|Modal false V|
|---|---|---|
|quadratic / threshold / omitted variable|1.00 / 1.00 / 1.00|—|
|dilution|.92|rational|
|cubic (odd term)|.24|power (.68)|
|interaction|.05|power (.95)|
|opposing curvature|.00|power+rational (.84)|
|collinear-quadratic|.03|rational (.61)|

Table 9 is Table 4 re-expressed as filtering errors: the power member inherits A's documented blind spots, so those structures stay in V under the power label instead of escalating. On the SRSD side the filtering remains descriptive (median |V| 1–2 for 11 of 14 equations; i.30.5 → {exponential} in 78% of datasets in the v4-era run — a row now annotated as withdrawn-task-descriptive, since the task itself failed re-run reproducibility; ii.38.14 escalates on 67%), because surface-form family labels dissolve under examination: the battery correctly rejects Planck-type "exponentials" that no linearizable family owns.

### 4.7 The gate comparison: why not just compare against a nested rival?

**Table 10 — S5 size on H0 pools (200 seeds/cell; naive and calibrated thresholds for E/F; unchanged from v4 — S5 uses only A/E/F).**

|H0 pool|A|E naive|E cal.|F naive|F cal.|
|---|---|---|---|---|---|
|synth n=50 (600 fits)|.080|.212|.050|.132|.050|
|synth n=200 (600)|.040|.183|.050|.135|.050|
|SRSD n=50 (2,600)|.051|.167|.050|.086|.050|
|SRSD n=200 (2,600)|.038|.152|.050|.123|.050|

**Table 11 — S5 structural-H1 detection at calibrated FPR 5% (mean over tasks; unchanged from v4).**

| |synth n=50|synth n=200|SRSD n=50|SRSD n=200|
|---|---|---|---|---|
|A (permutation curvature)|.774|.782|.348|.493|
|E (CV-quad)|.929|.934|.462|.628|
|F (AICc-quad)|.930|.940|.525|.669|
|OR(A, E) at 5% joint budget|.774*|.940*|.547*|.675*|

\* Joint-budget rows re-read per §4.7(d): at n = 50 the recalibrated union collapses to A alone.

Four findings. **(a) The obvious alternative is the better default detector, and the margin is exactly A's blind zone.** The quadratic expansion carries per-feature squares and cross-terms, so INTERACT (A .24/.26) and OPPOSING (A .03/.02) are caught at 1.00 by both rival gates, and the two SRSD extrapolation-blind equations (i.39.11, i.43.43) are detected at .99–1.00. The direction set, not the null construction, was the binding constraint. CUBIC curvature defeats all three gates equally (~.3): odd curvature is near-orthogonal to every even-column design — the honest general form of the blind-zone lesson. One re-run footnote lands here: the corrected Table 5 shows D beating F on trigonometric families (.86 vs F's bench-average lead), a reminder that S5's rival is the canonical objection, not a ceiling; §2's survey paragraph lists what would have to be tested to close this.

**(b) The rival gates do not come calibrated.** Because the rival nests the null, ties and noise break toward escalation: naive thresholds false-alarm at .15–.21 (CV) and .09–.14 (AICc), and calibration to 5% required an H0 reference pool — knowledge a production gate does not have. A's residual-permutation null is an approximation (§3.2) but was empirically size-valid on every pool we measured (.038–.080; synth n = 50 sits at .080 and is reported as measured) without any reference set. A's remaining structural advantage has narrowed twice in this paper's history — v4 cut it from "valid by construction" to "valid as measured," and this revision notes that D, on the corrected path, is nearly as usable without calibration (.06–.08 naive size, §4.2) — but it survives: no other gate in the comparison produced a size-valid p-value on all four pools with zero reference data.

**(c) Does the calibrated threshold transfer?** Calibrating E and F on all but one H0 equation and evaluating on the held-out equation (13-fold per bench × n) gives held-out sizes of .053/.050 (E, n = 50/200) and .051/.052 (F) on average — but with per-equation spread [.010, .095]. Earlier drafts emphasized the mean; the spread is the honest statistic. A reference pool of thirteen monomial equations calibrates a gate to within ±.04 of nominal *per equation*, and which side of nominal a given deployment lands on is luck. Reference-pool calibration is not an in-sample artifact, but it is also not a guarantee; pools should be larger than thirteen, and the released harness treats pool size as a deployment parameter, not a constant.

**(d) The OR union's recall advantage is a joint-budget artifact.** OR(A, E) with each gate at its own 5% runs a joint false-positive budget of .065–.102 across the four benches — its TPR edge over E is purchased with that larger budget, not with information. Recalibrating the combination to a 5% joint budget collapses the union: at n = 50 A's own size consumes the budget and the union is A (TPR .348 SRSD / .774 synthetic); at n = 200 the recalibrated union reaches .664/.940 against F's .669/.940 — parity with F, never a gain. At any common error budget F is the frontier; the union is defensible only when a workflow explicitly budgets for its larger naive size. Noise-layer scenarios leave all three gates quiet (.04–.07): the gates measure form, not noise assumptions. Cost in fit units: A ≈ 2 full fits + 1 vectorized pseudo-inverse (999 permutations); F = 2 fits; E = 12 — all negligible against the candidate flood a discovery pipeline already fits, which is the actual status quo's cost point: the workflow runs none of them today.

### 4.8 Real data: does the map predict live behavior?

**Table 12 — S6 World Bank relations: verdicts (median permutation p over 5 seeds for A; deterministic p for B–D; α = .05). C/D verdicts are SVD-path; every C/D verdict agrees with the Gram-era verdict on all 120 cells, as expected at cond(X) ≈ 10^5.7 — the artifact never reached S6.**

|Relation (target ~ features)|Years all four gates agree|Divergent cells (which gate, which way)|
|---|---|---|
|infant mortality ~ GDP pc|4/4 reject|—|
|CO₂ pc ~ GDP pc|4/4 reject|—|
|internet penetration ~ GDP pc|4/4 reject (p < 10⁻¹⁶)|—|
|GDP ~ population|4/4 quiet (A .71–.91; B .19–.33)|—|
|life expectancy ~ GDP pc|0/4|A alone rejects 2019 (.035) & 2022 (.008); B alone rejects 2022 (.0002); C/D never|
|life expectancy ~ GDP pc + health exp pc|2/4 quiet|B rejects 2022 (.0002) while A quiet (.108)|

Where signal is strong the gates agree — three relations rejected unanimously every year, one true scaling relation cleared unanimously — and A's five p-values per cell are Monte Carlo replicates of one estimator whose spread is small (median .009, max .094 at resolution 1/1000; no verdict flips beyond boundary cells). Divergence concentrates where the synthetic map puts weak curvature: the life-expectancy cells flip year to year, and the two-feature cell shows B rejecting where A stays quiet. We label this pattern consistent with the dilution mechanism (§4.2, Table 4) rather than a reproduction of it: on observational data the p-value change after adding a second feature can also reflect feature correlation, omitted-variable changes, sample composition, measurement error, cross-country heterogeneity, or spatial dependence — no ground truth exists to separate these, so the claim is signature-consistency, not mechanism confirmation.

**Table 13 — S6 dilution probe: the same target with and without a second correlated curved feature (3 years × 5 seeds; detection = median gate verdict; C/D on the SVD path).**

|Relation|A|B|C (SVD)|D (SVD)|
|---|---|---|---|---|
|cereal yield ~ fertilizer (1 feature)|rejects 3/3 years (.011–.031)|rejects 3/3 (~.003)|rejects 3/3 (.049/.076/.136)|rejects 3/3 (.049/.076/.136)|
|cereal yield ~ fertilizer + energy pc (2 features)|quiet 3/3 years (.067–.149)|rejects 2/3 (.025–.092)|quiet 3/3 (.109/.176/.226)|rejects 3/3 (.018/.0001/.0002)|

The probe's content changed in the re-run, in a way that strengthens the mechanism reading. v4 reported A silenced by the second feature while RESET variants kept alarming. On the corrected path C is silenced too: adding the correlated curved feature moves C's median p from .049–.136 (rejecting at two of three years) to .109–.226 (quiet everywhere), while D holds at ≤ .018 across years and features. The synthetic DILUTION story predicted exactly this split between fitted-power gates (which absorb the diluting feature's curvature into ŷ-powers) and built-direction gates (whose explicit columns survive it) — and now two of the three fitted-power gates show it live, not one. No ground truth exists on S6 and none is claimed; the contract of §3.6 is signature consistency in both directions: unanimous cells where the map predicts signal, and map-consistent divergence where it predicts weakness.

A companion engineering bench (App. E) pushes the same probe across nine physical-engineering sources (147 run cells): unanimous detection of a documented rating-curve misspecification, mechanism-echo divergences matching the synthetic map, and a numerical fingerprint — twenty-one C cells and ten D cells where the Gram path produced p-values below 10⁻²⁰ that the SVD path maps back to [0.07, 1.0] — that is now the bench's cleanest illustration of why App. D matters outside the laboratory.

### 4.9 Runtime

Per-dataset wall time (S1 means over 500-seed cells): A 12.0 ms (n = 50) to 25.3 ms (n = 1000); B 1.0–1.8 ms; C/D on the SVD path 0.36–0.43 ms per call at n = 200 against 0.33 ms for the Gram path (§3.7) — the correctness fix costs single-digit milliseconds and nothing else; the battery at n = 200 runs in 24 ms; S5 gates E/F 1.5–3.5 ms. All gates are effectively free at the scale of a discovery pipeline that already fits thousands of candidates; the bottleneck resource is not compute but the decision discipline — and, this revision adds, the solver discipline.

## 5 Analysis

**Q1. Does the low-df curvature gate beat RESET on univariate quadratic deviations?** No — it ties; both saturate from n = 50. What low degrees of freedom buy is calibration, not power — which is what makes A's calibrated numbers trustworthy rather than large. The asymmetries that do appear are structural (odd terms, interactions), where A's quadratic projection is mathematically the wrong instrument.

**Q2. Does averaging squared columns miss feature-specific curvature and interactions?** Yes, by construction — and the benchmark quantifies when that matters, in three escalating venues: synthetic power cells (Table 4), filtering errors (Table 9), and live observational data (Table 13). The dilution mechanism in particular is now a three-level result: constructed (S1), filtered (S4), and signature-consistent on observational data (S6) — with the re-run upgrading the third level from one gate showing the signature to two (A and the corrected C), while D's resistance to it is itself predicted by the same geometry (built columns keep the curved direction the ŷ-powers absorb). A synthetic blind spot stops being a laboratory artifact the moment it anticipates which real relation goes quiet — though on S6 that correspondence remains pattern-level, not mechanism-proof.

**Q3. Does HC3 fix false alarms under heteroskedasticity — and what did our original "42–56%" actually measure?** Heteroskedasticity was never the disease (Table 6): growing-variance noise leaves every gate at .048–.058, and every S3 heteroskedastic cell replicates its homoskedastic twin within .02. What Tables 1–2 of v2–v4 actually measured is now clear, and the measurement has three layers worth keeping separate. *Layer 1, the observation*: our harness's HC3-Wald path rejected at rates up to 1.00 on 5 of 13 physics-derived monomial tasks, flat in n. *Layer 2, the retraction*: a standard SVD-based implementation (statsmodels HC3-Wald, same F reference) is calibrated on the very tasks where our path rejected at 1.00 (.067 over 60 seeds each) and agrees with ours to ≥ 3 significant digits on well-conditioned tasks (e.g., p = .0906 vs .0906) — so the data was never diseased, and "the robust estimator is the failure mode," the sentence two drafts printed in bold, was false. *Layer 3, the localization* (App. D): the divergence tracks design conditioning (point-biserial r = .84 between task-median log₁₀ cond(X_aug) and rejection rate in our scan; the paper's earlier r = .78 used a coarser aggregation), and a 2×2 carrier decomposition isolates the sandwich's inverse factor as the sole culprit — swapping the hat vector between Gram and SVD sources changes nothing, swapping the covariance factor flips 1.00 to .00. The mechanism: forming X′X squares cond(X) (10^10.7 → 10^21.05 on the worst task), pinv's default truncation tolerance (max(n,k)·eps·σ_max ≈ 5.1×10⁻³) then zeroes a direction whose σ_min(X′X) ≈ 8.5×10⁻¹⁰ sits seven orders below the floor, and the Wald statistic, built through that truncated inverse twice, explodes. Equilibration fails because the direction is already dead in the Gram matrix's geometry, not in its column scales; rcond = 0 fails partially (rejection .30, not .00) because a 10²² condition solve keeps only ~6 digits regardless of truncation. Corrected conclusion, scoped honestly: on the pools we measured, the SVD path's HC3 is calibrated on physical monomials (.053–.060, task-bootstrap [.050, .058]); whether any residual HC3 miscalibration exists under conditions we did not sample — narrower domains, stronger collinearity, heteroskedastic errors on top of degenerate designs — is open, and a full re-derived sweep is exactly what §3.7's harness now makes cheap.

**Q4. Do alarms predict extrapolation error?** As binary evidence, usefully (2.7–3.7× enrichment on the 31-task pool, complementary alarm sets); as magnitudes or guarantees, no. The correct reading is asymmetric: alarm → form evidence; silence → nothing. The pool-sensitivity finding (§4.5) is itself part of the answer: one task moving in or out of the pool shifts the base rate by a quarter, so anyone quoting an enrichment factor should quote the pool it came from.

**Q5. Can the same machinery answer "not a power law — then what?"** Yes, as a viability set, after one repair the main benchmark could not have predicted: the permutation null's empirical calibration is coordinate-bound, and the wild bootstrap is the repair. The formal run delivers the contract (containment .93–.95, pooled size .0496, correct escalation 1.00 on quadratic/threshold/omitted-variable), with near-neighbour overlap reported rather than suppressed and A's blind spots wearing filtering clothes (Table 9). What the battery does not do is rank, decide, or empty the set — a viability filter hands the discovery layer a shortlist, and "routing" was our over-reach in earlier drafts.

**Q6. Why not just compare against a nested rival — and does the re-run change the answer?** Measured, it still wins, and the re-run strengthens the case. F (AICc-quad) led every S5 cell in v4 and those numbers are untouched (S5 never touched the broken path). What changed is the rest of the menu. The re-run rehabilitated D: on corrected numbers D is the strongest single gate on the physics bench at every n (.494/.675/.733 calibrated; §4.2), which v4 could not see because its D rows mixed detection with artifact. A's niche narrowed again, and honesty requires tracking the narrowing: v3 said A was the best-calibrated test; v4 said its niche was calibration without a reference pool; this draft notes D now nearly matches that niche (.06–.08 naive size on our pools, §4.2) while beating A on power everywhere. What keeps A's niche from closing entirely is that D's un-calibrated usability is measured on exactly the pools where we fixed the solver, and A's self-calibration mechanism (statistic and null share one pinv row) has a structural argument behind it, not just four pools of evidence. The deployment ordering this revision supports: a pipeline with a reference pool should default to F; one without should treat A as the only gate whose small-sample size has a mechanism-level argument and measured evidence on every pool, while watching the corrected D as the likely better gate the moment a modest H0 set exists; recall-critical deployments should budget for F, not for an uncalibrated union. The benchmark's own proposed test remains demoted from default to portfolio member — and the re-run's rehabilitation of D is the demotion's sharpest edge, because the gate our benchmark was built around is now, on corrected numbers, not even the best classical test on its home bench.

**Q7. Does any of this transfer to real data?** Where the synthetic map predicts strong signal, four gates built for physical laws agree on macroeconomic relations they were never designed for (Table 12) — cheap evidence that the map is not a laboratory artifact. Where the map predicts weakness, observational divergence appears with the predicted signature (dilution-consistency, Table 13) — consistent-with, not confirmation-of: too many observational confounds and no truth to adjudicate. The re-run's contribution to Q7 is negative and useful: at S6's conditioning (median cond(X) ≈ 10^5.7) the artifact never fired, so the Gram-era S6 verdicts stand unconditionally — a boundary drawn around the hazard from the observational side. What S6 and App. E cannot supply — rates, calibration, extrapolation truth, mechanism identification — they do not attempt; real data here is a compatibility probe, not a leaderboard, and their rejection counts are compatibility statistics, not performance rankings.

## 6 Limitations

**Single candidate family (main benchmark).** All gates of §3.2–3.5 evaluate a power law; the battery addresses filtering, not selection — ranking within V is the discovery layer's job and is not benchmarked. **One rival, one criterion each (S5).** The rival is the quadratic expansion and the criteria are 5-fold CV and AICc — the reviewer's canonical objection, not a survey of selection gates; LOOCV, held-out likelihood, or non-nested rivals may rank differently, and §2's survey paragraph names the literature S5 does not cover. **Battery scope.** n = 200 only; no member for interaction or trigonometric structures; SRSD side descriptive pending a symbolically exact linearizability registry. **Real-data scope (S6).** Six relations, one source, four years, no ground truth: no rate claims are possible, the relation list is prespecified to prevent cherry-picking but remains six choices, and year-to-year verdict flips in weak-curvature cells signal genuinely boundary-level signal rather than gate defect. **Withdrawn task.** feynman-i.30.5's non-reproducibility is bounded (sentinel flags on all 7,500 rows) but unexplained; if its data generation drifted in a way that also affected other tasks before the v2 run, we have no independent channel to detect it — the sentinel covers regeneration from the current clone, not archaeology of the original one. The 32-task re-analysis in App. B bounds the blast radius numerically (v4's A/B rows reproduce exactly from 32-task pools). **Fixed augmentation.** Built powers fixed at (2, 3). **Log-space scope.** All gates operate in log-log space; additive-error models are out of scope and marked not-applicable. **Statistical precision and clustering.** S1/S2/S4 rates rest on 500 seeds; S3 on 30 (directional); S5 on 200; S6's 5 seeds are Monte Carlo replicates of one estimator, not samples. Seeds repeat within equations, so fit-count intervals understate cross-equation uncertainty; Table 2 carries task-level bootstrap intervals and LOEO ranges (App. C), and the power tables (3, 5) now report per-gate thresholds — but power-table task-level intervals remain wider than the differences they would test (18 H1 tasks is a thin meta-pool), so we report LOEO sensitivity instead of pretending precision. **No multiple-comparison correction** across the 1,480 cells; rates are descriptive, and family-level set decisions are analyzed as sets (§4.6) rather than as four independent tests. **HC3 re-scope.** The SVD re-run covers C/D across S1–S3, S6, and the engineering bench; the statsmodels cross-check covers the fitted variant on two failing tasks (60 seeds) plus spot checks — a full independent-implementation sweep of every cell is released as a script but not executed end to end, and C/D numbers should be read as SVD-path properties until then. **Permutation validity.** The residual-permutation null is an approximation whose finite-sample behavior under leverage, dependence, or stronger heteroskedasticity is unmeasured; leverage-adjusted, wild, block, cluster, and fixed-design null variants are the required stress set. **OR budget.** Union gates are reported at both naive and joint-5% calibration; any deployment comparison must fix the budget first (§4.7d).

## 7 Reproducibility statement

Every generator is seeded from a documented rule embedding the grid position, so any cell of any table regenerates bit-identically — and this revision required that promise to be kept under stress, which is exactly how the one place it broke became visible (§3.3). The harness records content hashes of the code under test — including the exact form_check implementation — and of the harness itself; a statistic-mismatch sentinel verifies on every call that the wrapper's recomputation agrees with the native implementation. The sentinel's C/D gap is closed: each HC3 call is now cross-checked against a statsmodels reference on a rolling subsample, and a `gram_conditioning_hazard` registration fires above cond(X) = 10⁸ (P4). The battery is implemented at the experiment layer over the unmodified core code path, so audited core-code hashes remain valid. S5's gate definitions and grid are frozen in gate_compare.py's prespecification docstring (seed bases disjoint from all prior suites). S6 pulls from the World Bank API with every raw response cached under the run directory (runs/wb_real/raw/), reproducible against the cached snapshot and refreshable against the live API; the relation and probe lists are frozen in wb_run.py / wb_blindspot.py before any gate ran. The revision's re-analyses are scripted and released with their raw inputs: the C/D re-run driver (paper/revision_loc/rerun_cd_svd.py and replay_s6_smoke_svd.py) regenerates every SVD-path number in §4 from the released per-cell parquet/CSV files, with the sentinel statistics of §3.7 printed per suite; the calibration re-derivation (calibrated_power_v5.py), task bootstrap (task_bootstrap_v5.py), S3 family tables (s3_family_table_v5.py), and the i.30.5 diagnosis (diag_i305.py) regenerate the rest; the independent HC3 check uses statsmodels 0.15 on identical data frames. The engineering companion bench of App. E is released the same way (real_data/): raw snapshots cached with SHA-256 hashes and access dates (2026-09-06), a manifest of prespecified relations with per-relation notes and per-source license fields, and replayable fetch/prepare/smoke scripts whose smoke output path is now a parameter (an earlier hard-coded path briefly let a re-run overwrite the audited Gram-era file — itself fixed and disclosed in the released incident note); every declined evaluation carries a typed reason. The benchmark artifact (repository archive with all generators, per-cell results, and content-hashed configurations) is deposited with the publication; a permanent DOI is assigned at release, and the anonymous review copy ships as an anonymized archive link on the title page (§ Typesetting note in the supplement). All cells' raw results (parquet and CSV), failure registrations with typed reasons, selection and exclusion tables, analysis scripts, and figure sources are released. Zero-noise cells deliberately exercise the failure-registration path; the battery's 94,000 and S5's 20,400 datasets contain zero data-layer failures, which is itself a registered outcome.

## 8 Conclusion

We built the falsification-side benchmark that model discovery lacks: six suites, 508,825 generated datasets and 472,500 held-out evaluations, one protocol that puts every gate on a common basis, and full registration of every failure. This revision's deliverable is not a winning detector — it is a corrected map and a record. The record first: the benchmark's most striking claim, the one two drafts led with, did not survive contact with an independent implementation. Our 42–56% HC3 false-alarm story was a silent numerical failure of our own inference path — pinv's default truncation meeting a condition number that forming X′X had squared — and the retraction chain, from the statsmodels disagreement through the carrier decomposition to the boundary the re-run drew (synthetic wide-domain designs unaffected, physical monomials recalibrated, engineering bench fingerprinted by 31 p-value underflows, macro data untouched), is released scripted end to end. We print the retracted numbers nowhere in the main text except as contrasts, and the case study (App. D) exists because the next harness that forms X′X on a 10¹⁰-condition design deserves the warning we did not get.

The corrected map, in three sentences. Robust-covariance RESET variants, on a solver that does not self-destruct, are calibrated on real physical monomials and — in their exog form — are the strongest single detector on the physics bench (.494/.675/.733 calibrated at every n), a rehabilitation the artifact had hidden. The fitted-F gate stays nearly blind to rational laws, the one structural finding three drafts have now confirmed; no gate sees exponentials well; and the permutation gate, whose calibration without a reference pool once looked like a superpower, ends the paper where v4 left it — demoted to portfolio member, its niche narrowed by a rehabilitated rival, its value an argument about mechanisms rather than a leaderboard position. The nested quadratic rival adjudicated by AICc remains the best default where a reference pool exists (its calibration transfers across equations, .050–.054 on average — with per-equation luck spread of .010–.095 that we now show rather than smooth), the battery turns rejection into a viability shortlist once its null is rebuilt per coordinates, and on live data the map predicted where four gates would agree, which gate would go quiet on which relation, and — added by the re-run — that a second fitted-power gate joins the permutation gate in the dilution signature. Benchmark papers are usually measured by the leaderboards they produce. This one is measured by the failures it makes visible, including its own — and the one this revision made visible, in the spirit of its title: the artifact you silently accept, the one that flatters your framing and ships in bold, breaks exactly the extrapolation you built the benchmark to protect.

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

[21] Neyman, J. (1937). "Smooth test" for goodness of fit. *Skandinavisk Aktuarietidskrift*, 20, 149–199.

[22] Kallenberg, W. C. M., & Ledwina, T. (1997). Data-driven smooth tests when the hypothesis is composite. *Journal of the American Statistical Association*, 92(439), 1094–1104.

[23] Bierens, H. J., & Ploberger, W. (1997). Asymptotic theory of integrated conditional moment tests. *Econometrica*, 65(5), 1129–1151.

[24] Fan, J., Zhang, C., & Zhang, J. (2001). Generalized likelihood ratio statistics and Wilks phenomenon. *Annals of Statistics*, 29(1), 153–193.

[25] Utts, J. M. (1982). The rainbow test for lack of fit in regression. *Communications in Statistics — Theory and Methods*, 11(24), 2801–2815.

[26] Stute, W. (1997). Nonparametric model checks for regression. *Annals of Statistics*, 25(2), 613–641.

[27] Vovk, V., Gammerman, A., & Shafer, G. (2005). *Algorithmic Learning in a Random World*. Springer.

[28] Hart, J. D. (1997). *Nonparametric Smoothing and Lack-of-Fit Tests*. Springer.

## Appendix A — Suite composition

|Suite|Cells|Datasets per cell|Generated datasets|Extra evaluations|
|---|---|---|---|---|
|S1 synthetic|285|500|142,500|—|
|S2 SRSD-Feynman|480|500|240,000|240,000 test sets; 232,500 extrapolation probes|
|S3 family-stratified|392|30|11,760|—|
|S4 battery|188|500|94,000|—|
|S5 gate comparison|102|200|20,400|—|
|S6 real data + probe|33|5|165|—|
|Total|1,480|—|508,825|472,500|

"Generated datasets" = dataset generations (one per cell × seed); each typically receives four gate evaluations, ≈ 2.1M total (§3.8). Outside the seeded totals, the App. E engineering bench contributes 150 single-realization cells (147 run), one fit per gate per cell. C/D columns throughout this draft are SVD-path (§3.7); the Gram-era C/D numbers appear only inside App. D's case study and the App. B withdrawal analysis.

## Appendix B — The withdrawn task and what moved

**What happened.** During the C/D re-run, the B-sentinel (§3.7) flagged feynman-i.30.5 on all 7,500 rows: archived classical-F p-values do not reproduce from regenerated data at any seed or n, while every other S2 task reproduces to ≥ 6 significant digits. The upstream SRSD clone is the same commit (7d00b45, cloned 2026-09-04), the environment is identical (numpy 2.5.2, scipy 1.18.1, Python 3.13), and the task's generator is deterministic — so the drift is inside that task's data path, is bounded by the sentinel, and resisted our diagnosis. The task is withdrawn (H1 pool 19 → 18); its re-run rows ship with sentinel flags.

**What moved.** Everything quantified, nothing structural. H0 pools exclude the task, so Tables 1–2 are untouched. H1 power rows lose one of 19 tasks: raw n=200 pooled power moves A .516→.489, B .280→.253 (32- vs 31-task pools); calibrated power (Table 3) moves A .513→.486, B .281→.254 at n=200; the extrapolation base rate moves .198→.143 and Table 7's precision/recall shifts with it (§4.5). The withdrawn task's own C/D contrast matched its synthetic twin (S1 EXP_MODEL, both-solver difference zero), so no SVD-vs-Gram conclusion depends on it. What we cannot rule out: a drift in the *original* run's data generation for this task that also touched other tasks — the sentinel covers regeneration from the current clone, not archaeology of the original; the exact reproduction of v4's A/B rows from 32-task pools is the available bound.

## Appendix C — Task-level and headline intervals

**Headline Wilson 95% intervals.** S2 H0, n = 200 (26,000 fits per cell): A .053 [.050–.055]; B .050 [.047–.052]; C (SVD) .053 [.051–.056]; D (SVD) .056 [.054–.059]. S1 H0, n = 200 (4,500 fits per cell): A .037 [.032–.043]; B .041 [.036–.047]; C (SVD) .045 [.039–.051]; D (SVD) .042 [.036–.048]. Battery power-family size: .0496 [.046–.053]. S5 H0 pools (600–2,600 fits per cell): intervals within ±.013 of every Table 10 entry. Intervals for all other cells are computable from released per-cell counts.

**Task-level uncertainty (S2, n = 200, 13 H0 equations).** Seeds repeat within an equation, so fit-count intervals understate cross-equation uncertainty; the equation is the honest unit for generalization over physical laws. Resampling the 13 equations (10,000 bootstrap replicates): A .053 [.051–.056]; B .050 [.048–.052]; C (SVD) .053 [.050–.056]; D (SVD) .056 [.055–.058]. Leave-one-equation-out pooled rates: A [.052, .054]; B [.049, .050]; C (SVD) [.052, .054]; D (SVD) [.056, .057]. On the corrected path the per-equation spread has collapsed to a band — C ranges .052–.057 across equations and D .054–.062 — where the Gram-era table showed five equations at .998–1.00 against eight at nominal. For the power tables, the honest statement is LOEO sensitivity: the H1 meta-pool (18 tasks) is too thin for bootstrap intervals narrower than the effects at issue, and per-gate thresholds τ are reported alongside power for that reason (Table 3). Where a single held-out H0 equation calibrates E or F (§4.7c), held-out sizes spread [.010, .095] around a .050–.054 mean — reported as a spread, not smoothed to the mean.

## Appendix D — Case study: a silent numerical failure of pinv-squared-conditioning inference

**The symptom.** v2–v4's Table 2: our HC3-Wald path rejecting at .998–1.00 on five physics-derived monomial tasks (i.34.27, i.38.12, ii.27.18, iii.12.43, iii.15.14) while sitting at .048–.060 on the other eight, flat in n from 50 to 1000, with clean Gaussian noise and clean data by construction. The pattern looked like a discovery: "the robustness correction is itself the failure mode." It was our solver.

**The disagreement.** On the two worst tasks, 60 seeds each at n = 200: our path rejects at 1.00; statsmodels' HC3-Wald (SVD on X, same F reference) rejects at .067. On well-conditioned tasks the implementations agree to ≥ 3 significant digits (p = .0906 vs .0906), so both implement the same definition and the divergence is numerical.

**The localization.** Three experiments (paper/revision_loc/locate_pinv_failure.py; all on the harness's own code path, official seed rules, unmodified generators). *E1, reproduction*: same data, four paths — harness HC3 1.000/1.000/0.050 (two failing tasks + one control), classical F .050/.033/.067, statsmodels .033/.033/.050, our own SVD-native rewrite .033/.033/.050. The failure follows the implementation, not the data. *E2, carrier decomposition*: a 2×2 swap of hat-vector and covariance-factor sources (Gram vs SVD) shows the covariance factor is the sole carrier — gram_hat + gram_cov rejects 1.00, svd_hat + gram_cov rejects 1.00, gram_hat + svd_cov rejects .00, svd_hat + svd_cov rejects .00. The hat vector, computed through the same broken pinv(X′X), is innocent. *E2b, mechanism adjudication*: on one task and ten seeds — default pinv(G): reject 1.00 (median p ≈ 0); pinv(G, rcond=0): reject .30; column-equilibrated pinv(G): reject 1.00; SVD on X: reject .00 (median p .63). The numbers, on the worst task: cond(X) = 10^10.7, cond(X′X) = 10^21.8, σ_min(X′X) = 8.5×10⁻¹⁰, pinv's default truncation floor = max(n,k)·eps·σ_max ≈ 5.1×10⁻³. The smallest direction sits seven orders below the floor and is silently annihilated; the Wald statistic then travels through the truncated inverse twice and detonates. *E3, the boundary*: across all 13 H0 tasks, rejection rate tracks conditioning (point-biserial r = .84 against task-median log₁₀ cond(X_aug)); every task above 10^7.6 fails, everything below is nominal, and S1's wide-domain synthetic designs (cond ~ 10¹) never fire at all — which is why this draft's Tables 1, 4, and 6 are byte-identical to v4's.

**Why each tempting fix fails.** Column equilibration: the direction is annihilated in the Gram matrix's spectrum, not by column scaling — equilibrating changes nothing. rcond = 0: recovers partially (.30 rejection, not .00) because inverting a 10²²-condition matrix in float64 keeps only ~6 digits no matter the tolerance; the solve itself is the hazard. Blaming pinv: wrong on both ends — numpy's pinv is itself an SVD (gelsd), and an explicit inverse of X′X dies identically; the hazard is forming X′X at cond(X) ≳ 10⁸ with any default-tolerance inverse. A's immunity: method A's permutation vectorization also uses a pinv row, but its statistic and its null share that row, so truncation cancels between them — calibration survives on the same designs, which is evidence for the carrier mechanism, not luck.

**The fingerprint, and the remediation.** On the engineering bench (App. E), the Gram path left 31 p-value underflows (p < 10⁻²⁰: 21 C cells, 10 D cells) that the SVD path maps back to [0.07, 1.0] — a signature that can be screened for post hoc in any repository that still ships the Gram path. The fix is algorithmic, not parametric: factorize X directly (SVD or QR), never form X′X at high condition numbers; cost, 0.36 vs 0.33 ms per fit. The harness additions are a `gram_conditioning_hazard` registration above cond(X) = 10⁸ and a rolling statsmodels cross-check on C/D (§7) — closing the sentinel gap that let this artifact pass two review cycles in silence.

**What this case does and does not establish.** It establishes: a concrete, scripted, seeded instance of a known numerical hazard manufacturing a flat-in-n, headline-grade statistical claim; the carrier (the sandwich's inverse factor), the trigger (default truncation), the boundary (cond(X) ≳ 10⁷·⁶ in our pools), and two failed fixes. It does not establish that HC3 is universally calibrated on physical designs — our corrected evidence covers these pools and solvers; narrower domains, stronger collinearity, or heteroskedasticity layered on degenerate designs are unmeasured (§6). The claim "HC3 fails on physical laws" is retracted; the claim "research code that inverts X′X through default tolerances silently fails exactly where this benchmark lives" stands, and is the version worth citing.

## Appendix E — The real-engineering data bench (companion to S6)

S6 probes mechanism transfer on macroeconomic aggregates. The deployment domain these gates are built for is physical-engineering measurement — fluid, civil, hydraulic, energy, meteorological, seismic, chemical, and aerospace data, where power laws are the daily candidate models. This appendix documents a companion bench of that kind, run under the same contract as S6: real data has no ground truth, so the bench registers compatibility and divergence, not error rates; its cells are single realizations and stand outside Appendix A's seeded totals.

Composition and discipline. Nine open sources under their stated licenses (redistribution rights checked per dataset; per-source license fields in the released manifest): UCI, EPA fueleconomy, Open-Meteo/ERA5 (reanalysis, marked), OpenML, USGS NWIS rating curves, USGS FDSN, RDatasets, NIST StRD, NIST WebBook; 91 datasets; 150 prespecified candidate log-log relations; 761,515 raw rows. Raw bytes are cached (~115 MB) with SHA-256 hashes; cleaning is replayable and imputation-free — rows with missing or non-positive values are dropped and logged, and the few unit conversions are explicit in the manifest (e.g., Whiteside °C → K). Simulator datasets that violate the real-data discipline (a marine-engine and two wave-energy numerical models) are excluded and listed, as are WebBook fluid pages below the 10-row rejection-detection floor; one gauge whose instantaneous series was never published is logged as an empty skip.

**Table E1 — Composition by source.**

|Source|Datasets|Relations|Rows|
|---|---|---|---|
|USGS NWIS rating curves|12|12|442,269|
|UCI|11|38|145,920|
|EPA fueleconomy|1|4|48,624|
|Open-Meteo / ERA5|10|20|87,600|
|OpenML|3|8|21,997|
|USGS FDSN|1|3|12,334|
|RDatasets|14|26|1,349|
|NIST StRD|13|13|968|
|NIST WebBook|26|26|454|
|Total|91|150|761,515|

Protocol. Each cell applies the §3.2 gates unchanged at α = 0.05 — A with its B = 999 permutation p, B as registered, C/D on the SVD path — on the positivity-filtered relation, subsampled to n ≤ 800 (seed 42) to bound permutation cost. The run is a smoke pass: one realization per cell, compatibility plus first verdicts, not a rate estimate. Of 150 cells, 147 run; three are registered insufficient_n by design (Bennett5 and Roszman1 have no positive rows after the log-space positivity filter; Danwood has n = 6) — the bench's honest-failure probes working as intended.

**Table E2 — Smoke verdicts on the 147 run cells (C/D SVD-path).**

|Gate|Rejections|Registered failures|
|---|---|---|
|A — permutation curvature|97/147 (.660)|10 declines, all checked_false_small_n (n ≤ 11)|
|B — RESET fitted F|111/147 (.755)|0|
|C — RESET fitted HC3 (SVD)|93/147 (.633)|0|
|D — RESET exog HC3 (SVD)|87/147 (.592)|3: insufficient_df (n = 13, four features), hc3_leverage_degenerate (n = 546), rank_deficient|

A and B agree on 123/147 cells (92 both reject, 31 both quiet; 5 A-only, 19 B-only). The failure columns are themselves informative: D's three declines are exactly the three HC3 failure modes of Table 2's taxonomy, now on observational data; A's ten are uniform declines-to-run on n ≤ 11 cells; B and C never fail. Nothing fails silently — every declined evaluation is registered with a typed reason.

The Gram-era bench overstated C and D (rejection .735/.639), and the mechanism of the overstatement is legible in the re-run's fingerprint: 31 cells (21 C, 10 D) where the Gram path returned p < 10⁻²⁰ — numbers that are not evidence of anything except a truncated inverse — map under SVD back into [0.07, 1.0]. Datasets named in v4's divergent-cell prose are re-read on corrected values: the CCPP proxy-curvature cell (PE ~ V, where V correlates .81 with the strongly curved AT) keeps its reading — C(SVD) still rejects at ≈ 0, B at 8×10⁻⁴ — because that cell's divergence was mechanism, not noise; the automobile horsepower cell (A .878 quiet, B 1.2×10⁻⁴ reject under heteroskedastic body-style scatter) is likewise unchanged. Cells whose v4 divergence *was* the artifact — the WebBook enthalpy family, housing (D's leverage-degeneracy cell now registered rather than rejected), petrol and whiteside under D — are re-read as ordinary nulls or registered failures, and their v4 prose is superseded.

Designed regimes. The bench's design document assigns expected regimes to 66 of the 150 cells before the smoke pass (Table E3). The assignment is a reading aid, not a scored test: the H0-side laws are approximate by design, and at n up to 800 an approximate law is a detectable null. The 12 rating-curve cells — whose registered form Q ~ h omits the documented stage offset h₀ of the true rating law Q = C(h − h₀)ⁿ, a mismatch noted in the manifest — are rejected unanimously by A, B, and C(SVD) (12/12) and by D on 11/12 (usgs_14105700 is the miss, p_D = 1.0 against p < 10⁻⁷ for every other gate on the same station; we read it as a genuine power boundary of D's built-column geometry on that station's curvature shape, not a numerical event — the station passes every other gate). The five H0-side cells A retains are those with no designed departure mechanism — the braking-distance v² law, the tree allometry, the road-load relation — plus the two wind-shear sites with the smallest fitted curvature (Altamont Pass p = .229, Baltic Rügen p = .411).

**Table E3 — Designed-regime cells, rejections at α = 0.05 (C/D SVD-path).**

|Designed regime|Cells|A|B|C (SVD)|D (SVD)|
|---|---|---|---|---|---|
|Approximate power laws (H0-side)|55|43|48|43|43|
|Known non-power-law physics (H1-side)|50|27|32|21|20|

(Counts follow the released design document's dataset-level labels; v4's 30/36 split used a different cell granularity and is superseded.) On the H1 side no gate sweeps the board: the misses concentrate in WebBook enthalpy cells whose measured temperature window keeps the local power-law approximation good — e.g., toluene Δ_vapH ~ T: all four gates quiet at n = 26 — a mutual-miss cell recorded as measured.

**Table E4 — Named cells: mechanism echoes and divergences; p-values at α = 0.05 in A / B / C(SVD) / D(SVD) order, — = declined with registered reason.**

|Cell|A / B / C / D|Reading|
|---|---|---|
|Wind shear v₁₀₀ ~ v₁₀, 10 ERA5 sites (n = 800 ea.)|8/10 / 10/10 / 9/10 / 7/10|B flags every site (p ≤ 2×10⁻⁴); A retains exactly the two smallest-curvature sites (.23, .41) — a graded response, not a threshold artifact; C(SVD) retains Baltic Rügen (.119)|
|Gust factor ~ v₁₀, same sites (n = 800 ea.)|10/10 ×4|unanimous rejection at every site; no truth to adjudicate|
|puromycin rate ~ conc (Michaelis–Menten), n = 23|.011 / .069 / .169 / .169|A detects the saturation curvature, B misses at α = .05 — Table 3's small-sample power gap, live|
|stormer Time ~ Viscosity + Wt (calibration law), n = 23|.001 / 7×10⁻⁴ / 9×10⁻³ / —|single-feature cells quiet (A .94, B .99); the second feature carries the interaction the law contains; D registers rank_deficient|
|EPA CO₂ ~ MPG (stoichiometric inverse, exponent −1), n = 800|.030 / 1.5×10⁻⁴ / .727 / .727|A and B reject on a near-exact physical law bent by fleet drive-mix curvature|
|EPA combined MPG ~ displacement (road load), n = 800|.386 / .596 / .599 / .599|designed power law passes|
|CCPP PE ~ V, n = 800|.642 / 8×10⁻⁴ / ≈0 / .023|proxy curvature: V correlates .81 with the strongly curved AT, so B/C's fitted-power columns inherit AT's curvature while A's V-only direction set cannot — the mirror of dilution|
|automobile horsepower ~ engine size, n = 203|.878 / 1.2×10⁻⁴ / .094 / .094|divergence under heteroskedastic body-style scatter|
|rock perm ~ area + peri, n = 48|.835 / 1×10⁻⁵ / 5.6×10⁻⁶ / .596|A/D quiet, B/C reject; no truth to adjudicate|
|housing price ~ 4 features, n = 546|.723 / .683 / ≈0 / —|D's leverage-degeneracy failure mode, registered live at n = 546, while A/B/C return verdicts|

Scope. What the bench supports is what S6 supports: compatibility at scale (147/150 cells run across nine sources, zero silent failures), failure registration on observational data, and mechanism-echo cells whose directions match the synthetic map. Its per-gate rejection totals (97/111/93/87 of 147) are compatibility statistics on single realizations without ground truth — they rank nothing, and reading them as performance would repeat exactly the leaderboard confusion the benchmark exists to correct. What it cannot support: rates, calibration, or rankings — single realizations, no ground truth, and time-ordered NWIS/Open-Meteo series whose permutation-independence assumption is an approximation marked in the manifest. All four p-values per cell are released (real_data/ with the Gram-era and SVD-era result files distinguished by suffix) with the manifest (real_data/clean/manifest.json, SHA-256 per source) and replayable fetch/prepare/smoke scripts.
