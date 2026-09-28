# -*- coding: utf-8 -*-
"""apply_batch4.py — 分析/结论区修订: HC3 撤回与重定性、OR 预算、真实数据口径、
任务级不确定性入附录 C、参考文献修复、变更记录。"""
p = 'draft_v3_complete.md'
md = open(p, encoding='utf-8').read()
subs = [
 # ---------- 摘要 Second(在 batch1 基础上再更新为独立复核结论) ----------
 ("""**Second**, heteroskedasticity-robust (HC3)
RESET variants — the default "safe" choice in mainstream packages — are
severely miscalibrated on our physics-derived monomial tasks: pooled
false-positive rates of 42–56%, flat in n, driven by 5 of 13 equations
whose built designs are near-collinear (per-task rejection ≈1.00) while
the other 8 remain near nominal (.048–.060).""",
  """**Second**, the benchmark caught — and its revision retracted — its own
most striking claim: our harness's HC3-Wald path rejects at rate 1.00 on
5 of 13 physics-derived tasks whose built designs are near-collinear
(pooled 42–56% false alarms, flat in n), but an independent SVD-based
implementation (statsmodels) is calibrated (.067) on the same data while
agreeing with ours to three digits elsewhere: the failure is a silent
numerical pathology of pinv-squared-conditioning inference paths, not of
the robustness correction — a deployment hazard discovery code should
test for."""),
 # ---------- §4.1 尾段: 更新为复核结论 ----------
 ("""Second, the claim we defend is the conditional one: HC3-Wald variants
are severely miscalibrated under the evaluated physics-derived designs;
we do not claim miscalibration under every design a user could
construct.""",
  """Second, an independent-implementation check run in revision (Section 5,
Q3) reattributes the failure: a standard SVD-based HC3-Wald
implementation (statsmodels) is calibrated on the very tasks where our
pinv-based path rejects at 1.00 (.067 over 60 seeds each), and the two
implementations agree to ≥3 significant digits on well-conditioned
tasks. Tables 2, 3, and 6 therefore report properties of the harness's
HC3 path — a real, silent, deployment-relevant hazard, but no longer
evidence against the HC3 correction itself."""),
 # ---------- Q2 观察措辞微调 ----------
 ("""routed (S4), and observed (S6). A synthetic blind spot stops being
a laboratory artifact the moment it predicts which real relation silences
which gate.""",
  """routed (S4), and signature-consistent on observational data (S6). A
synthetic blind spot stops being a laboratory artifact the moment it
anticipates which real relation goes quiet — though on S6 that
correspondence remains pattern-level, not mechanism-proof."""),
 # ---------- Q3 整段重写 ----------
 ("""**Q3. Does HC3 fix false alarms under heteroskedasticity?**
Heteroskedasticity was never the disease here (Table 6). What breaks on
real physical monomials is HC3 itself (Table 2): the correction presumes
a design whose leverage structure is informative; on near-collinear built
columns it multiplies noise into a statistic whose reference distribution
is wrong at every n measured. The robust estimator is the failure mode.""",
  """**Q3. Does HC3 fix false alarms under heteroskedasticity — and was the
observed failure HC3's or ours?** Heteroskedasticity was never the
disease (Table 6). Table 2 shows *our* HC3-Wald path breaking on the
near-collinear built designs of 5 of 13 physics-derived tasks; the
original draft attributed this to the correction itself ("the robust
estimator is the failure mode"). The independent-software check required
in revision retracts that attribution. On the two worst tasks, 60 seeds
each at n = 200: our pinv(X′X)-based path rejects at 1.00 while a
standard SVD-based implementation (statsmodels HC3-Wald, same F
reference) rejects at .067 — essentially nominal. On well-conditioned
tasks the two implementations agree to ≥3 significant digits (e.g.,
p = .0906 vs .0906), confirming the definition, not the data, is
implemented identically and the divergence is numerical. Our path's
rejection tracks design conditioning (point-biserial r = .78 with
log₁₀ condition number; top-decile conditioning rejects at 1.00,
bottom-decile at .05). Corrected conclusion: a silent numerical failure
of pinv-based inference on squared-condition-number designs — cond(X′X)
squares cond(X), and near-collinear power-augmented designs sit exactly
where that squared conditioning meets tolerance cutoffs. The practical
hazard is real (research code commonly pinv-solves; the false alarms are
total and silent), but it is an implementation robustness lesson, not
evidence against HC3; whether any HC3 miscalibration survives a
re-implemented sweep is open, and C/D rows in Tables 2, 3, and 6 should
be read as path-specific until then."""),
 # ---------- Q6 OR/校准迁移段改写 ----------
 ("""The purchase price is a threshold: naive E/F false-alarm at .09–.21 until
an H0 reference pool calibrates them, while A is size-valid by
construction. The portfolio reading follows: a pipeline that can afford a
reference pool (or a ~10% naive false-alarm budget) should default to the
rival gate; one that cannot, or that must report a valid p-value per
candidate, should default to A; recall-critical deployments should run
both — OR(A, E) is the best S5 cell on three of four benches.""",
  """The purchase price is a threshold: naive E/F false-alarm at .09–.21 until
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
should budget for F rather than pay for an uncalibrated union."""),
 # ---------- Q7 收尾句已有 consistency 措辞(§4.8), 这里改 Q7 段 ----------
 ("""**Q7. Does any of this transfer to real data?** Where the synthetic map
predicts strong signal, four gates built for physical laws agree on
macroeconomic relations they were never designed for (Table 13) — cheap
evidence that the map is not a laboratory artifact. Where the map
predicts weakness, live divergence follows the predicted mechanism
(dilution, Table 14). What S6 cannot supply — rates, calibration,
extrapolation truth — it does not attempt; real data here is a
mechanism-validity probe, not a leaderboard.""",
  """**Q7. Does any of this transfer to real data?** Where the synthetic map
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
rankings."""),
 # ---------- §6 Limitations 更新 ----------
 ("""**Statistical precision.** S1/S2/S4 rates rest on 500
seeds; S3 on 30 (directional); S5 on 200; S6 on 5 seeds per cell with
p-stability, not rates, as the reported quantity.""",
  """**Statistical precision and clustering.** S1/S2/S4 rates rest on 500
seeds; S3 on 30 (directional); S5 on 200; S6's 5 seeds are Monte Carlo
replicates of one estimator, not samples (median p spread .009, max
.094). Seeds repeat within equations, so Table 2's Wilson intervals
understate cross-equation uncertainty; task-level bootstrap intervals
and leave-one-equation-out ranges (Appendix C) are the honest units for
generalization over equations, and the HC3 failure is concentrated in 5
of 13 equations. No multiple-comparison correction is applied across the
1,480 cells; rates are descriptive, and family-level set decisions are
analyzed as sets (§4.6) rather than as four independent tests.
**HC3 single-replication.** The independent statsmodels check covers the
fitted-HC3 variant on two failing tasks (60 seeds) and spot checks
elsewhere; a full re-implementation sweep across C/D, all suites, and
the calibrated-power rows is outstanding. **Permutation validity.** The
residual-permutation null is an approximation whose finite-sample
behavior under leverage, dependence, or stronger heteroskedasticity is
unmeasured here; leverage-adjusted, wild, block, cluster, and
fixed-design null variants are the required stress set. **OR budget.**
Union gates are reported at both naive and joint-5% calibration; any
deployment comparison must fix the budget first (§4.7(d))."""),
 # ---------- §7 复现: 增补修订工件 ----------
 ("""frozen in `wb_run.py` / `wb_blindspot.py` before any gate ran. The
engineering companion bench of Appendix E is released the same way
(`real_data/`): raw snapshots cached with SHA-256 hashes, a manifest of
preregistered relations with per-relation notes, and replayable
fetch/prepare/smoke scripts; every declined evaluation carries a typed
reason.""",
  """frozen in `wb_run.py` / `wb_blindspot.py` before any gate ran. The
revision's re-analyses are scripted and released with their raw inputs:
`paper/revision_analysis.py` regenerates every number in §4.1's
task-level paragraph, §4.6's set-level errors, §4.7's OR-budget and
transfer analyses, and §4.8's Monte-Carlo spreads from the released
per-cell parquet/CSV files; the independent HC3 check (Q3) uses
statsmodels 0.14 on identical data frames. The engineering companion
bench of Appendix E is released the same way (`real_data/`): raw
snapshots cached with SHA-256 hashes and access dates (2026-09-06), a
manifest of preregistered relations with per-relation notes and
per-source license fields, and replayable fetch/prepare/smoke scripts;
every declined evaluation carries a typed reason. The benchmark artifact
(repository archive with all generators, per-cell results, and
content-hashed configurations) is deposited with the publication; a
permanent DOI is assigned at release."""),
 # ---------- §8 结论: HC3 与 OR 句 ----------
 ("""Robust-covariance
RESET variants are structurally miscalibrated on real physical laws
(42–56% false alarms, flat in n).""",
  """Its most striking
claim did not survive its own revision: the 42–56% HC3 false alarms were
a silent numerical failure of a pinv-based inference path on
near-collinear built designs — an SVD-based implementation is calibrated
exactly where it fails — a retraction that is itself the benchmark
working as designed."""),
 ("""the nested quadratic
rival adjudicated by AICc is the strongest default gate on both benches —
paying for it with a naive 9–14% false-alarm rate until a reference pool
calibrates it — while the permutation gate is the only one
size-valid out of the box, and their union is the best recall on real
equations.""",
  """the nested quadratic
rival adjudicated by AICc is the strongest default gate on both benches —
paying for it with a naive 9–14% false-alarm rate until a reference pool
calibrates it (a calibration that transfers across equations) — the
permutation gate was empirically size-valid on every pool measured
without one, and the OR union's apparent recall premium is a joint-error
budget artifact that vanishes at a common budget."""),
]
missed = []
for old, new in subs:
    if old not in md:
        missed.append(old[:70])
    else:
        md = md.replace(old, new)
open(p, 'w', encoding='utf-8').write(md)
print('missed:', missed or 'NONE')
