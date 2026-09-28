# Extension: Why Not Just Compare Against a Nested Rival? (gate comparison)

**Annotated extension v1 — to be merged after Q5 (as Q6), with a forward
pointer from §6 Limitations. Status: formal-run-complete (19 synthetic
scenarios + 32 SRSD equations × n ∈ {50, 200} × 200 seeds = 20,400 datasets,
zero failures; preregistered in `gate_compare.py`'s docstring before the
runs; numbers cite `runs/gate_cmp/`).**

📝 **Placement craft.** This is the benchmark's most obvious attack —
"anyone would just fit a quadratic rival and compare" — answered *with data
instead of prose*, and it is where detector A (form_check) formally cedes the
default-gate position. A section that answers its own strongest objection
with a table is worth two pages of defense elsewhere; the repositioning
paragraph below must survive into the camera-ready, because the paper's
credibility rests on the benchmark being able to defeat its own hero.

## Q6. Is the permutation gate dominated by the "obvious" alternative — model comparison against a nested quadratic?

The four tests of §3.2 are all *test statistics*. The workflow alternative
that requires no new statistics is a **model-comparison gate**: fit the
power law and a nested rival — the full quadratic expansion of the log
features (linear terms, per-feature squares, pairwise cross-terms) — and
escalate if the rival wins. We benchmark three gates under identical
conditions (n ∈ {50, 200}, rel = 3%, 200 seeds, preregistered):

- **G1 = A (form_check):** permutation curvature p < .05; size valid
  by construction (the null is resampled from the data itself).
- **G2 = CV-quad:** 5-fold cross-validated log-RMSE; escalate if
  RMSE_cv(power) − RMSE_cv(quad) > δ. Naive δ = 0 and δ\* calibrated on
  an H0 pool to FPR = 5% (per bench × n).
- **G3 = AICc-quad:** escalate if AICc(power) − AICc(quad) > m; naive
  m = 0 and calibrated m\*. k = column rank (the expansion is rank-deficient
  by construction under collinearity).

Structural-H1 detection at calibrated FPR = 5% (mean over tasks; SRSD = 19
non-monomial equations, synthetic = 11 structural scenarios):

| | synth n=50 | synth n=200 | SRSD n=50 | SRSD n=200 |
|---|---|---|---|---|
| A (form_check) | .774 | .782 | .348 | .493 |
| CV-quad | .929 | .934 | .462 | .628 |
| AICc-quad | **.930** | **.940** | **.525** | **.669** |
| OR(A, CV-quad) | .937 | .941 | .547 | .675 |

Three findings, stated in the order a reviewer will care about:

**(a) The obvious alternative is the better default detector, and the margin
is exactly A's blind zone.** The quadratic expansion carries per-feature
squares *and cross-terms*, so INTERACT (A: .24/.26) and OPPOSING (A:
.03/.02) are caught at 1.00 by both rival gates, and A's two known SRSD
extrapolation-blind equations (i.39.11, i.43.43; §4.3) are detected at
.99–1.00. The direction-set, not the null construction, was the binding
constraint. CUBIC curvature defeats *all three* gates equally (~.3): odd
curvature is near-orthogonal to even directions in every even-column design
— a property of any fixed direction set, which is the honest general form
of the blind-zone lesson.

**(b) The rival gates do not come calibrated.** Because the rival nests the
null, ties and noise break toward escalation: naive thresholds false-alarm
at .15–.21 (CV, δ = 0) and .09–.14 (AICc, m = 0). Calibration to 5%
required an H0 reference pool — knowledge a production gate does not have.
A's permutation null is size-valid out of the box (.038–.080 measured;
synth n=50 sits at .080, slightly over nominal and reported as measured).
This is A's remaining structural advantage: *no reference set needed*.

**(c) Complementarity is real but thin.** OR(A, CV) adds +.085 over CV
alone at SRSD n=50 and +.047 at n=200 (equations where curvature is weak in
log-log space but permutation power survives), and nothing on synthetic
(+.008, where both gates saturate). Noise-layer scenarios (t(3),
heteroskedastic, AR(1), outliers) leave all three gates quiet (.04–.07) —
the gates measure form, not noise assumptions, and the §3.3 semantic
separation survives its most obvious stress test.

📝 **Repositioning craft (must survive revision).** The recommended-gate
sentence for §5/§7, verbatim intent: *default AICc-quad if a reference pool
or ~10% naive FPR is acceptable; A if an out-of-the-box valid p-value is
required; OR both if recall is the binding constraint.* A is repositioned
from "the gate" to "the self-calibrating member of the portfolio" — the
same portfolio argument §4.4 already made against single-detector
superiority, now applied to A itself. The paper's contribution claims are
unchanged: the fitted-but-wrong phenomenon (§2), the HC3 trap (§4.1), and
the benchmark itself all survive A losing this head-to-head; what would
*not* survive is hiding the head-to-head.

## Scope and what this section does not claim

- **One rival family, one decision rule each.** The rival is the quadratic
  expansion (the reviewer's canonical "obvious fix"); a rational/log rival
  battery is Q5's machinery, not this section's. The gates answer one
  binary question — "is the power law adequate?" — and do not rank models;
  selection among viable families remains the discovery layer's job (Q5's
  division of labor, unchanged).
- **Calibration asymmetry is a workflow property, not a statistical
  objection.** An H0 pool is cheap to simulate for the *power-law* null
  (that is literally A's permutation null); the point is that G2/G3 must
  spend it while A spends nothing, at the price of lower power everywhere
  it overlaps the rival.
- **No extrapolation claim.** i.39.11/i.43.43 are *detected* by the rival
  gates here; whether alarm→correction reduces out-of-domain error is
  §4.5's extrapolation axis and is not re-measured in these coordinates.
- **Cost accounting is fits, not wall-clock.** A ≈ 2 full fits + 1
  vectorized pinv (999 permutations); AICc-quad = 2 fits; CV-quad = 12.
  All three are negligible against the family battery or symbolic
  regression — which is the actual status quo's cost point: the workflow
  runs *none* of them today (§2's silent-acceptance finding).

## References to add (beyond the draft's list)

- Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and
  Multimodel Inference: A Practical Information-Theoretic Approach*
  (2nd ed.). Springer. (AICc)
- Arlot, S., & Celisse, A. (2010). A survey of cross-validation procedures
  for model selection. *Statistics Surveys*, 4, 40–79.

📝 **Citation craft.** Two method citations only (AICc, CV), for the same
reason as Q5: the section introduces no new statistic, so it cites the
procedures it benchmarks, not neighbors. The nesting tie-bias finding
(needs no reference) is stated as a property of nested comparisons — if a
reviewer names prior art on nested-model false-alarm rates, the text
already concedes the mechanism.
