# Prospective Structure Test + Real-Data Case Studies (2026-09-09)

Two additions answering the reviewer suggestions: (A) a prospective
blind test of §5 Q3's claim that silence is predictable from symbolic
structure before running anything; (B) real-data case studies grounding
the 7.1% silent-tail figure outside simulation.

---

## A. Prospective structure test (§5 Q3 made falsifiable)

**Claim under test.** §5 Q3: "a practitioner can predict silence from
structure before running anything: if the residual story of a candidate
form is 'curvature that oscillates' or 'a cross-term', a four- or
five-gate panel will not save it."

**Rules fixed before scoring** (structure-only, no silence values):

1. **trig_O1** — the equation contains sin/cos/tan with a non-suppressed
   (non-tiny-coefficient) appearance → the in-domain curvature is large →
   **visible**.
2. **suppressed_trig** — trig appears only inside tiny-coefficient terms
   (relative magnitude < 1e-5) → oscillation present but invisible in-domain
   → **blind**.
3. **cross_term** — some Add contains a product of ≥ 2 distinct symbols →
   ŷ-powers, regressor powers, and the quadratic rival all see the same
   cross-terms only at second order → **blind**.
4. **small_correction** — any Float coefficient < 1e-5 multiplying a
   nonlinear correction (the relativistic sqrt(1−β²) family) → curvature
   second-order in-domain → **blind**.
5. otherwise → **plain**, visible.

**Script:** `blind_structure_test.py`. Equations parsed from the upstream
SRSD-Feynman registry (`srsd_data.load_feynman()`), not from any
silence-carrying artifact.

**Score against measured per-task silence (n = 50, noise-present):**

| | predicted blind (13) | predicted visible (5) |
|---|---|---|
| mean measured silence | **.515** | **.024** |
| min / max | .029 / .849 | .000 / **.090** |

- 11/18 exact hits at the .5 cut; **zero false negatives** — every task
  with measured silence > .5 is predicted blind.
- Point-biserial r = .665 (p = .0026).
- The rule's errors are one-directional: the 7 false positives sit in the
  paper's middle band (.16–.39), where §4.3 already locates gate-specific
  partial detection; the worst over-prediction is i.43.43 (.029, the
  near-linear-in-log family whose curvature is real but too mild at n = 50).
- The separation is the headline: predicted-blind tasks average .52
  measured silence; predicted-visible tasks never exceed .09.

**Conclusion:** §5 Q3's prediction claim survives falsification. Structure
classifies the blind spots with no misses; the uncertainty is one-sided
(structure over-predicts blindness in the partial-detection band, never
under-predicts it). For a practitioner: a cross-term or suppressed
oscillation in the candidate's residual story means the panel cannot see
the misfit; no structural red flag means the panel has never been observed
to miss above 9%.

**Output:** `runs/joint_silence/blind_test_per_task.csv`.

---

## B. Real-data case studies

**B1. Concrete curing (uci_concrete_strength, Yeh 1998; n = 1030).**
Deployed relation: strength ~ age. The 28-day strength is the design
standard, so fitting a power law on ages ≤ 28 d is exactly what practice
does; curing is a saturating process, so the power law must over-predict
at long ages. Fit: bootstrap n = 100 from the age ≤ 28 pool (300 fits);
far-field truth: **real measured strengths at ages 90–365 d** (190
specimens).

| verdict (4-gate panel, α = .05) | rate | median far-field error | P(err > 10%) |
|---|---|---|---|
| **jointly silent** | **.863** | **+27.9%** | **1.00** |
| flagged | .137 | +31.4% | 1.00 |

86.3% of fits receive a full silent certificate on real data; every
silent fit carries a far-field error above 10% (median +24.5% signed
over-prediction — the saturation direction, invisible to every gate and
to the panel jointly). This is the paper's simulated tail (7.1% > 10% at
n = 50) realized at full strength outside simulation: the silent
certificate is not conservative, it certifies a +25% over-prediction of
real 90–365-day strengths.

Note the honest asymmetry: here even flagged fits carry the same
far-field error (the alarm fires on in-domain structure, not on the
far-field failure mode). Alarm-informativeness on real data is weaker
than the simulation's 7.6–157× enrichment; silence and alarm are both
uninformative about this failure mode — consistent with the paper's
"silence means nothing" reading, and worth stating.

**B2. Mercury vapor pressure (rdataset_pressure, Weast handbook; n = 19).**
Power law fitted in the low-T half, validated against the real high-T
measurements (T = 220–360 °C, Clausius–Clapeyron curvature): the panel
detects (silence 0/199 evaluable bootstrap fits; median far-field error
77%). A negative control showing the panel does fire on large in-domain
curvature — the blind-spot analysis is about direction, not about gates
being globally asleep.

**B3. USGS rating curve (11447650; n = 34,928).** Power law on the
bottom 90% of stage, validated on real top-decile flows: detected
(0 silent; median far-field error 116% — the datum-offset failure mode
is panel-visible). Second negative control.

**Provenance:** `real_data/concrete_case_study.py` (+ `_rows.csv`),
`real_data/pressure_case_study.py`, `real_data/usgs_case_study.py`; all
fits via `common.run_all_methods` (identical gate paths as the
simulation), far-field scored against held-out real measurements.

---

## Suggested paper integration

1. **§5 Q3**: append the prospective-test paragraph (blind-predicted mean
   silence .515 vs visible .024; FN = 0; r = .665, p = .003).
2. **New §4.6 or §4.7 "Real data"**: the concrete case as the lead
   (86.3% silence, all silent fits > 10% far-field error), the two
   negative controls in one paragraph each.
3. **§6 Limitations**: note the real-data alarm-informativeness asymmetry
   (flagged fits carry the same far-field error as silent ones in B1 —
   alarm informativeness measured in simulation does not automatically
   transfer to every real failure mode).
