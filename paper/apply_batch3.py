# -*- coding: utf-8 -*-
"""apply_batch3.py — §4 结果区修订(HC3 任务级/电池集合层/OR 预算/S6 口径)。"""
p = 'draft_v3_complete.md'
md = open(p, encoding='utf-8').read()
subs = [
 ("""**Table 2 — S2 true physical monomials (13 tasks), false-positive rate at α = 0.05 (26,000 fits per cell).**""",
  """**Table 2 — S2 physics-derived monomial tasks (13 SRSD-Feynman equations; synthetic draws from real equations, not measurements), false-positive rate at α = 0.05 (26,000 fits per cell).**"""),
 ("""On synthetic Gaussian data, all four tests are calibrated (Table 1). On
real physical laws the picture splits (Table 2): A and B remain
calibrated (.049–.053, flat in n); both HC3 variants break identically at
every sample size — C at .417–.422, D at .563–.565, intervals that do not
move from n = 50 to n = 1000. This is not a small-sample artifact that
more data cures; it is a structural miscalibration of the robustness
correction itself on near-collinear built columns (Section 4.4 rules out
noise as the cause; the data in Table 2 are clean by construction — what
differs is the *design*).""",
  """On synthetic Gaussian data, all four tests are calibrated (Table 1). On
the physics-derived designs the picture splits (Table 2): A and B remain
calibrated (.049–.053, flat in n); both HC3 variants break identically at
every sample size — C at .417–.422, D at .563–.565, intervals that do not
move from n = 50 to n = 1000. This is not a small-sample artifact that
more data cures; it is a structural miscalibration of the robustness
correction itself on near-collinear built columns (Section 4.4 rules out
noise as the cause; the data in Table 2 are clean by construction — what
differs is the *design*).

The pooled rates also conceal the failure's granularity, which matters
for any claim of generality over "physical laws": at n = 200 the HC3
failure is equation-specific, not uniform. Five of the 13 tasks
(i.34.27, i.38.12, ii.27.18, iii.12.43, iii.15.14) reject at .998–1.00
for C while the other 8 sit at .048–.060, essentially nominal. The
Wilson intervals of Table 2 treat all 26,000 fits as independent draws;
because seeds repeat within a fixed equation, the honest uncertainty
unit for generalization over equations is the task. Task-level bootstrap
(resampling the 13 equations; Appendix C) gives C a 95% interval of
[.197, .707], and leave-one-equation-out moves its pooled rate over
[.368, .448], against [.053, .054] for A. Two qualifications follow.
First, "42–56% on physical laws" is a pooled summary of a bimodal
phenomenon — near-nominal calibration on most equations, total failure
on a specific subset — and the open question is what makes that subset
special (Section 5, Q3, shows the answer tracks design conditioning).
Second, the claim we defend is the conditional one: HC3-Wald variants
are severely miscalibrated under the evaluated physics-derived designs;
we do not claim miscalibration under every design a user could
construct."""),
 ("""With the permutation null, the battery's own rational member rejects its
own true family a third of the time — nothing about rational laws is at
fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The
wild bootstrap repairs exactly that and nothing else.""",
  """With the permutation null, the battery's own rational member rejects its
own true family a third of the time — nothing about rational laws is at
fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The
wild bootstrap repairs exactly that and nothing else.

**Set-level errors and multiplicity.** The battery's decision is the
viable *set*, so the errors that matter are set-level, and each member
runs at α = .05 without correction — a union decision whose two error
types move oppositely under multiplicity. On H0-side cells the reported
containment is one minus the exclusion rate of the true family's own
member (.0496–.066 across the four truth families; the viable set is
never empty), and per-member α does not inflate that rate because
exclusion depends on one member only. Inflation appears on the other
side: some wrong member survives alongside the truth in 23.3% of
power-truth cells (near-neighbour overlap: exponential/logarithmic
members retain power data) and 92% of logarithmic-truth cells (rational
retains log data) — reported as near-neighbour overlap rather than
suppressed in the released tables. On the structural-H1 pool (72,000
cells whose truth is no battery member by design) some family is
retained in 100% of cells: routing there means *which wrong family is
modal* (Table 10), and a genuinely empty set never occurs — the battery
is a candidate filtering stage whose surviving set the discovery layer
must still rank, not a final selection. Downstream selection error after
the battery is out of scope here and unmeasured."""),
 ("""have. A's permutation null is size-valid out of the box (.038–.080
measured; synth n = 50 sits at .080 and is reported as measured). This is
A's remaining structural advantage: no reference set needed. **(c)
Complementarity is real but thin.** OR(A, E) adds +.085 over E alone at
SRSD n = 50 and +.047 at n = 200, and nothing on synthetic (+.008, both
gates saturate). Noise-layer scenarios leave all three gates quiet
(.04–.07): the gates measure form, not noise assumptions, and the §4.4
separation survives its most obvious stress test.""",
  """have. A's residual-permutation null is an approximation (§3.2) but was
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
E's threshold tightened on the H0 pool) collapses the union: at n = 50
A's own size consumes the budget and the union *is* A (TPR .348 SRSD /
.774 synthetic); at n = 200 the recalibrated union reaches .664/.940
against F's .669/.940 — parity with F, never a gain. The corrected
reading: at any common error budget F is the frontier; the union is
defensible only when a workflow explicitly budgets for its larger naive
size. Noise-layer scenarios leave all three gates quiet (.04–.07): the
gates measure form, not noise assumptions, and the §4.4 separation
survives its most obvious stress test."""),
 ("""Where signal is strong the gates agree — three relations rejected
unanimously every year, one true scaling relation cleared unanimously —
and A's permutation p is stable across seeds (max spread ±0.05, no
verdict flips). Divergence concentrates exactly where the synthetic map
puts weak curvature: the life-expectancy cells flip year to year (signal
at the detection boundary), and the two-feature cell shows B rejecting
where A stays quiet — the dilution mechanism (§4.2, Table 4) operating
with the sign the synthetic suite predicts for A's mean-square column.""",
  """Where signal is strong the gates agree — three relations rejected
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
mechanism confirmation."""),
]
missed = []
for old, new in subs:
    if old not in md:
        missed.append(old[:60])
    else:
        md = md.replace(old, new)
open(p, 'w', encoding='utf-8').write(md)
print('missed:', missed or 'NONE')
