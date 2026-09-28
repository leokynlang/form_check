# Silent Acceptance: A Falsification-Side Benchmark of Form-Specification Tests for Log-Linear Physical Model Discovery

**Annotated draft v0 — not for submission.**
Blockquotes marked 📝 are *craft notes* (why this section exists, what reviewers look for); they will be stripped before any submission.

📝 **Title craft.** Benchmark papers live or die on the framing of the gap. "Silent acceptance" names the *failure mode being measured* (models that are confidently wrong), not the tool being sold. This matters: our own literature audit (REFERENCES.md) killed the claim "novel detector" — form_check is Ramsey-RESET-family math — so the paper claims a *missing measurement axis*, not a new statistic. Reviewers forgive a modest claim; they hunt an inflated one.

## Abstract

Symbolic regression and physical model discovery pipelines certify candidate models with goodness-of-fit metrics (R², CV error, information criteria) that measure *how well a model fits*, never *whether an alternative family would fit strictly better* — and a model can score near-perfectly on all of them while being the wrong functional form, failing catastrophically the moment it extrapolates. Statistical form-specification tests exist (Ramsey's RESET and its descendants), but the symbolic regression community has no benchmark measuring whether they work *on the data distributions such pipelines actually consume*. We introduce a falsification-side benchmark: three task suites (19 controlled scenarios; 32 noise-injected SRSD-Feynman tasks derived from official configurations; 30 Feynman equations stratified into six structural families), four spec tests evaluated under a single protocol — identical significance level, identical permutation null for all tests, and power compared only at a common *empirical* false-positive rate of 5% — with every non-finite, rank-deficient, or degenerate evaluation explicitly registered rather than silently dropped. Under this protocol on 640,000+ fitted datasets, we find that (i) a classical F-based RESET variant that dominates on controlled quadratic benchmarks is nearly blind to rational-function ground truths (power 0.13), while the same test on synthetic benchmarks appears optimal — aggregate rankings do not transfer across task families; (ii) heteroskedasticity-robust (HC3) Wald variants of RESET, which appear strongest in raw tables, carry a 42–56% false-positive rate on true monomial physical laws that does not vanish with sample size, erasing their apparent advantage once tests are compared at matched empirical FPR; (iii) a design-based permutation curvature test achieves the best-calibrated size (≈5%) across both synthetic and real-equation pools but misses exponential-family laws; and (iv) no test is uniformly best — the failure profiles are complementary, arguing for portfolio-style deployment and for benchmarks that report per-family rejection rates rather than macro-averages. All artifacts, including per-equation results and failure registrations, are reproducible from fixed seeds and content-hashed harness configurations.

📝 **Abstract craft.** Four sentences of setup, one sentence of protocol, then findings as (i)–(iv) with *numbers in the abstract itself* — benchmark papers must quantitatively preview, because the whole contribution is empirical. Note the deliberate ordering: finding (i) is "rankings do not transfer", which is the paper's thesis, not any method's win. We never say form_check "wins" — it doesn't (it misses exponentials); it has *the best-calibrated size*. The phrase "falsification-side" anchors against SRBench (a discovery-side benchmark) in the first sentence of positioning.

## 1 Introduction

A regression model can be simultaneously excellent on every certification metric a discovery pipeline computes and wrong in its functional form. A power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes k-fold cross-validation, wins its information-criterion comparison, reports small standard errors — and extrapolates with −20% error one decade outside the training domain (this failure motivated the present work; §5.3 reproduces it as a scenario). Fit quality and form correctness are different axes, and only the first is measured by current practice.

The statistics literature has owned the second axis for over fifty years: Ramsey's RESET [Ramsey, 1969], its built-regressor and exog-power variants, and robust covariance corrections [White, 1980; MacKinnon & White, 1985]. Yet these tests are absent from symbolic regression (SR) benchmarking. SRBench [La Cava et al., 2021] and SRSD [Mundhenk et al., 2022] measure whether discovery methods *find* the right expression; nothing in the pipeline measures whether an *attached specification test* would *reject* the wrong one — i.e., the falsification side. The asymmetry has consequences: a deployed discovery system with no working falsification layer silently accepts misfit forms, and its downstream extrapolations inherit the error.

This paper asks a narrow, falsifiable question: **when the ground truth is known and the model form is wrong, which specification tests detect it — and at what false-alarm cost on the data distributions physical model discovery actually consumes?** Answering it requires resolving a protocol problem that, to our knowledge, no prior work addresses for this test family: the tests reject on incommensurable bases (a bootstrap t-threshold versus an asymptotic F p-value), so naive comparisons — e.g., "alarm when |t| > 3" against "reject when p < 0.05" — measure threshold placement, not power. Our protocol puts all tests on one permutation null, one significance level, and compares power only at a common empirical FPR.

📝 **Intro craft.** Three moves: (1) a *concrete* motivating failure with a number (-20%), never an abstract lament; (2) position against the closest existing benchmark by naming what it measures and the exact axis it omits — "discovery-side vs falsification-side" is the sentence reviewers will remember; (3) state the protocol problem as a *contribution* — admitting that naive comparison is invalid (we did exactly that in our first pilot) is a strength here, it motivates the calibration machinery. The research question is bolded and narrow: "which tests, at what false-alarm cost" — falsifiable by construction.

## 2 Related Work

**Specification testing.** RESET regresses y on X plus powers of the fitted values and tests the added coefficients' joint significance [Ramsey, 1969]. Robust inference under heteroskedasticity uses sandwich estimators [White, 1980], with HC3's leverage weighting the recommended small-sample correction [MacKinnon & White, 1985; Long & Ervin, 2000]. These results are derived for linear models under fixed designs; their behavior on log-transformed physical laws with heteroskedastic multiplicative noise is exactly what our benchmark measures, and §5.2 shows the small-sample behavior on real physical constants is far worse than the linear-model theory suggests.

**SR benchmarks.** SRBench evaluates discovery quality (exact recovery, R² on held-out data) across 252 problems [La Cava et al., 2021]; SRSD redesigns Feynman-derived tasks with realistic sampling ranges and metrics for scientific discovery [Mundhenk et al., 2022]. Neither evaluates an attached falsification layer: a discovery method that returns a wrong expression and a pipeline that flags it score identically. Our task generators deliberately reuse SRSD's official per-variable sampling objects so that our H0/H1 pools are the same distributions the discovery-side benchmark uses.

**Physical model discovery with honesty constraints.** Dimensional-analysis-assisted SR [Udrescu & Tegmark, 2020; Bakarji et al., 2022] uses units to constrain search. Dimensional correctness is necessary but not sufficient for form correctness — every ground-truth law in our suites is dimensionally self-consistent, yet is H1 (mis-specified) for a power-law model. Our benchmark is the complement: it measures the *rejection* behavior of tests, including their blind spots (§5), rather than discovery success.

📝 **Related work craft.** Three paragraphs, each ending with a sentence that *contrasts* rather than surveys: "they measure X; we measure Y." The middle paragraph does double duty — citing SRSD for the sampling objects also justifies our data-generation choices (official ranges, not invented ones). The last paragraph pre-empts a reviewer trap: "dimensional correctness ≠ form correctness" kills the objection "but units already catch this" before it's raised. Note what we do NOT cite: any "novel test" literature — because we are not claiming one, and citing near neighbors would create the obligation to out-novelty them.

## 3 Benchmark Design

### 3.1 Design principles

P1. **Common rejection basis.** All four tests decide at the same α = 0.05 on p-values; the permutation curvature test's empirical p-value uses B = 999 permutations of the same residual-permutation null its native bootstrap uses.
P2. **Power only at matched empirical FPR.** Headline comparisons calibrate each method's threshold to a 5% empirical FPR estimated on pooled true-monomial tasks, per sample size and per noise regime. Raw-α tables are reported alongside for continuity with prior practice.
P3. **Preregistered labels.** H0 (strict monomial) versus H1 (any non-monomial structure) is fixed from expression structure before any test runs; no label is adjusted after seeing results.
P4. **Explicit failure registration.** Zero-residual degeneracy (noise = 0), design rank deficiency, leverage degeneracy (1 − h → 0), singular covariance, and statistic mismatch are recorded per test per dataset; none is dropped silently. The zero-noise grid exists precisely to verify this: all four tests register `degenerate_zero_residual` there rather than emitting arbitrary p-values.
P5. **Official data configurations.** SRSD tasks reuse the official per-variable sampling objects (log-united DefaultSampling domains, linear SimpleSampling, integer domains, sign flags). Tasks whose inputs or outputs are not strictly positive under the official configuration are excluded *as not-applicable* — we do not shift data to make them log-fittable, since that changes the hypothesis under test.

### 3.2 Methods under test

- **A. Permutation curvature (form_check).** Adds one mean-centered quadratic column to the log-log design; the added coefficient is the statistic; the null distribution is the residual-permutation bootstrap of the same design (200 native draws for the reported t, 999 for the empirical p).
- **B. RESET fitted, classical F.** Adds ŷ², ŷ³ to the regression; joint F-test, powers (2, 3) fixed.
- **C. RESET fitted, HC3-Wald.** Same regressors; Wald statistic with HC3 sandwich covariance, referred to F(q, n−k).
- **D. RESET exog, HC3-Wald.** Powers of each *regressor* (not fitted values) as built regressors; HC3-Wald.
Configurations (powers, B, α, grids) were frozen before the main runs and are identical across all three task suites.

### 3.3 Task suites

- **S1 Synthetic control (19 scenarios × 500 seeds).** Three monomial H0; eleven structural H1 (quadratic, cubic, interaction, threshold, exponential, logarithmic, rational, omitted variable, dilution, opposing curvature, collinear-quadratic); five noise-layer H1 (heavy-tail t(3), growing variance, point outliers, group outliers, AR(1)) whose *mean structure is a strict monomial* — rejections there count as false alarms attributable to the noise assumption, never as detections. Grid: n ∈ {50, 200, 1000}, relative noise ∈ {0, 0.5%, 1%, 3%, 5%}.
- **S2 SRSD-Feynman (32 tasks × 500 seeds).** From the 120-equation registry, 32 tasks pass the positivity screen (13 monomial H0, 19 non-monomial H1); multiplicative noise injected in log space at the same grid; 2,000-sample same-distribution test sets and ×3-domain extrapolation probes accompany every training fit.
- **S3 Family-stratified Feynman (30 equations × 30 seeds).** Six structural families (exact power law, fractional power, exponential, rational, trigonometric, multivariable-additive) plus not-applicable representatives; homo- and heteroskedastic noise; sampling ranges scaled ×0.5 and ×2 to test sensitivity to the official domains.
Total: 643,250 fitted datasets across 1,167 grid cells.

📝 **Design-section craft.** The principles are *numbered and testable* — each one exists because violating it produces a specific wrong answer (P1 exists because our own first pilot made that mistake; say so in the talk, reviewers love a benchmark whose protocol was forged on a real error). S1's noise-layer scenarios carry an explicit semantic note: rejections there are false alarms *by construction* — this is the sentence that proves the benchmark distinguishes "no evidence against the form" from "form verified", which is the epistemics most benchmark papers leave implicit. The count at the end of 3.3 (643,250) is a scale claim reviewers can check against released artifacts.

## 4 Results

### 4.1 Size: false-positive rates on true monomials

Table 1 (synthetic, Gaussian noise; 500 seeds/cell) and Table 2 (S2, 13 true physical monomials; 500 seeds/cell) report raw-α rejection rates on H0.

| Table 1 — Synthetic H0 FPR | n=50 | n=200 | n=1000 |
|---|---|---|---|
| Permutation curvature (A) | .047 | .037 | .041 |
| RESET fitted F (B) | .045 | .041 | .038 |
| RESET fitted HC3 (C) | .048 | .045 | .038 |
| RESET exog HC3 (D) | .054 | .042 | .040 |

| Table 2 — SRSD H0 FPR (13 physical monomials) | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | .049 | .053 | .052 |
| B | .050 | .050 | .053 |
| C | **.422** | **.417** | **.418** |
| D | **.564** | **.563** | **.565** |

All four tests are calibrated on synthetic Gaussian data (Table 1). On real physical laws, A and B remain calibrated while both HC3 variants break — C stabilizes near 42%, D near 56%, flat in n. The mechanism: on these tasks the fitted-value powers ŷ², ŷ³ are nearly collinear with the design (physical constants span few log-decades), so the sandwich estimator's small-sample guarantee does not apply; S3 confirms the same signature (H0 exact-power-law: C .52, D .36) and adds a degenerate case (a 1/3-power monomial where D's built columns are effectively constant and its false-positive rate reaches 1.00). Heteroskedastic, heavy-tailed, outlier, and AR(1) noise layers leave A–D sizes at .03–.07: violations of the noise assumption are *not* what breaks these tests — the design is.

📝 **Results craft — Table 2 is the paper's center of gravity.** Two tables with the same H0 semantics, different data: identical behavior on synthetic, divergence on real physical constants. That contrast *is* the finding "small-sample theory does not transfer to real physical designs", and it required no new math to discover — only measurement. Note the bolding discipline: bold marks the pathological numbers, never our own wins. The final sentence assigns cause (collinearity of built columns) and points to S3's degenerate case (FPR 1.00) as the limiting signature — concrete, checkable, and it sets up the portfolio argument.

### 4.2 Power: raw and calibrated

At α = 0.05 on S2 H1 (n = 200): A .52, B .28, C .63, D .94 — but C and D pay Tables 2's false-alarm rates for that power. Calibrated to a common 5% empirical FPR, the ordering inverts at small n and re-inverts at large n:

| Table 3 — SRSD power @empirical FPR 5% | n=50 | n=200 | n=1000 |
|---|---|---|---|
| A | **.41** | **.51** | .59 |
| B | .19 | .28 | .42 |
| C | .05 | .28 | .41 |
| D | .19 | .37 | **.65** |

On S1 structural H1 (rel = 1%, n = 200, raw α): quadratic, threshold, exponential-family shapes are caught by everything (≥.99); the interesting cells are the blind spots — Table 4.

| Table 4 — S1 blind-spot cartography (n = 200, rel = 1%) | A | B | C | D |
|---|---|---|---|---|
| cubic (odd symmetry, large n shrinks quadratic projection) | .40 | 1.00 | 1.00 | 1.00 |
| interaction (lx1·lx2) | .21 | 1.00 | 1.00 | .10 |
| opposing curvature (+c, −c, exchangeable design) | **.02** | .98 | .98 | 1.00 |
| dilution (one curved variable of two) | 1.00 | .88 | .87 | 1.00 |

A's average-of-squares column is exactly blind when two variables' curvatures cancel (OPPOSING) and weak when a single odd-cubic term's quadratic projection shrinks with n; B's fitted-value powers absorb interaction locally (that is what a Taylor expansion does), which is the same mechanism that makes B nearly blind to rational laws in S3 (.13–.19) while dominating synthetic quadratic benchmarks. **No ranking in Table 3 or 4 survives transfer across task families.**

### 4.3 Family-stratified results and sensitivity

On S3 (Table 5, official domains, pooled noise): A misses the exponential family (calibrated power .11–.13 — an exponential law's curvature in log-log space is modest at official scales); calibrated trigonometric power collapses for C (.08) and D (.05) — their raw power was pure false-alarm budget. Sampling-range scaling (×0.5, ×2) shifts no family/method cell by more than 15 percentage points except rational×A (.66→.84 at ×2: larger domains make curvature easier to see — the expected direction); heteroskedastic versions of every cell reproduce the homoskedastic ones within .02. **Conclusions are insensitive to official-range perturbation and to the variance model.**

### 4.4 Do alarms predict extrapolation failure?

For every S2 fit we evaluate a same-distribution test set (2,000 samples) and a ×3-domain extrapolation probe against the closed-form truth. Base rate of far-field relative error >10%: .20. Alarm precision: A .39 (recall .63), B .55 (recall .49) — both enrich the base rate 2–2.8×, and their caught sets are largely *disjoint* (portfolio effect). The alarms are not calibrated to error magnitude: a rejection means "form evidence exists", not "extrapolation will be bad", and the two blind-spot tasks of §4.2's S2 table (i.39.11, i.43.43: near-zero in-domain curvature, 98% far-field error) are missed by every test. **No rejection is evidence of extrapolation safety.**

📝 **Results-section craft.** Order is fixed: size first (you may not talk about power until you've priced false alarms), then raw power *immediately undercut* by calibrated power in the same subsection, blind spots as their own table (they are the contribution, not an embarrassment — a benchmark that found no failure modes would be unpublishable), sensitivity as a robustness paragraph, and the extrapolation question last with the negative result stated as a negative result. Every table caption carries the seed count; every surprising number (.40 for cubic×A, 1.00 for the degenerate monomial) is explained mechanistically in the following sentence, because an unexplained number is a reviewer question, an explained number is a finding.

## 5 Analysis: the four questions

**Q1. Does the low-df curvature test beat RESET on univariate quadratic deviations?** On synthetic QUAD-like truths both saturate (1.00 at n ≥ 50). The asymmetries appear elsewhere: cubic odd terms (A .40 vs B 1.00) and interaction terms (A .21 vs B 1.00, D .10). Low degrees of freedom buy calibration, not power.

**Q2. Does averaging squared columns miss feature-specific curvature and interactions?** Yes, by construction — and the benchmark quantifies exactly when it matters: cancellation (OPPOSING .02) and dilution are detectable blind spots in S1, but *neither structure appears in the 32 real physical laws*, where A has no family below .11. Synthetic blind spots are real properties of the statistic, not of its use on real data; a benchmark needs both suites to say that.

**Q3. Does HC3 fix false alarms under heteroskedasticity?** Under synthetic heteroskedastic noise, sizes stay at .03–.07 for all methods — heteroskedasticity was never the problem. On real physical monomials HC3 itself is the failure (Table 2): the fix for the wrong disease, applied where the design breaks its assumptions.

**Q4. Do alarms predict extrapolation error?** As events, usefully (2–2.8× enrichment, complementary between A and B). As magnitudes, no — and the two worst extrapolators in S2 are missed by all four tests. The honest summary: specification tests price *form evidence*, not *field behavior*.

📝 **Analysis craft.** Each Q gets a one-paragraph answer whose first sentence is the answer. Q2 and Q4 contain the paper's epistemic fingerprint: Q2 refuses to over-read synthetic blind spots ("real properties of the statistic, not of its use"); Q4 states a negative result plainly. Benchmarks that only confirm wins get asked "what did you find that we didn't want to find?" — have the answer ready, in the analysis section, not in a rebuttal.

## 6 Limitations

Single-model-family tests only (all four tests evaluate a *power law* against alternatives; we do not benchmark model *selection* between families). Built-regressor powers fixed at (2, 3); higher-order or adaptive constructions may close B's rational blind spot at calibration cost we have not measured. All tests operate in log-log space by design; additive-error tasks are out of scope and marked not-applicable rather than log-shifted. Dimensional information was not manipulated as a factor (§5 of S3 report records the honest blank). Single-noise-draw-per-seed pairing; per-cell seed counts (500/30) give Wilson intervals within ±.02–.04 on headline rates, but family-level S3 cells (30 seeds) support directional claims only.

## 7 Reproducibility statement

All task generators draw from seeded NumPy generators with a documented seed rule embedding the grid position (any cell is independently reproducible); harness content hashes (including the tested `form_check` implementation) are recorded in every run's `env.json`; failure registrations, per-equation tables, raw parquet/CSV, and analysis scripts are released together. Zero-noise cells deliberately exercise the failure-registration path (all four tests register degeneracy there rather than emitting p-values).

## References

- Ramsey, J. B. (1969). Tests for specification errors in classical linear least-squares regression analysis. *JRSS B*, 31(2), 539–551.
- White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4), 817–838.
- MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. *Journal of Econometrics*, 29(3), 305–325.
- Long, J. S., & Ervin, L. H. (2000). Using heteroskedasticity consistent standard errors in the linear regression model. *The American Statistician*, 54(3), 217–224.
- La Cava, W., et al. (2021). Contemporary symbolic regression methods and their relative performance. *NeurIPS Datasets and Benchmarks*.
- Mundhenk, T. N., et al. (2022). Rethinking symbolic regression datasets and benchmarks for scientific discovery. *arXiv:2206.10540*.
- Udrescu, S.-M., & Tegmark, M. (2020). AI Feynman: A physics-inspired method for symbolic regression. *Science Advances*, 6(16).
- Bakarji, J., et al. (2022). Discovering governing equations from partial measurements. *Nature Computational Science* (data-driven dimensional analysis).
- Cranmer, M. (2023). Interpretable machine learning for science with PySR and SymbolicRegression.jl. *SciPy*.

---

## 📝 写作总纲（交稿前删除）

1. **主张结构**：全文只 claim 三件事——(a) 一个缺失的测量轴（falsification-side），(b) 一套协议（common-α + empirical-FPR calibration + failure registration），(c) 一张失效模式地图。每个 claim 都有两套独立任务池背书。这比"我们提出了更好的检验"安全得多——后者需要打赢所有 baseline，前者只需要测得准。
2. **数字纪律**：摘要、正文、表格的每个数字都能溯源到 `runs/full_500`、`runs/srsd_full`、`runs/feynman_families` 的文件；写论文时不允许凭记忆填数（我们拉数的那两步 terminal 调用就是论文附录"reproducibility"的素材）。
3. **负结果的放置**：外推预测的"不成立"部分写在 §4.4 正文，不藏 limitation——benchmark 论文的可信度来自它敢于让自己的方法难看。
4. **表格的视觉语法**：粗体只标病理数字或每行最优（不标自己）；每张表后一句机制解释；表注永远带种子数。
5. **篇幅分配**：Intro 1 页、Design 1.5 页、Results 2.5 页、Analysis 1 页、Limitations 0.5 页——Results 最重，因为这是 D&B 赛道；如果投 ML4PS workshop（9-12），可把 S3 压成半页。
6. **还没写但投稿前必须补的**：算法伪代码框（四检验统一伪代码）、S1–S3 的完整格数表（附录）、Wilson CI 列（正文表格加两列或用脚注）、`validation/` 目录重排后的路径修正。
