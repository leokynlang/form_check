# -*- coding: utf-8 -*-
"""apply_batch5.py — 第二轮修订: D 独立验证入文、.067 精度措辞、H1 配对
bootstrap、数值机制精确化、Table 3 限定语、附录 C/D 增补、§6/§7 更新。"""
p = 'draft_v3_complete.md'
md = open(p, encoding='utf-8').read()
subs = [
 # ---------- 摘要 Second: .067 → 名义相容 ----------
 ("""while an independent SVD-based
implementation (statsmodels) is calibrated (.067) on the same data while
agreeing with ours to three digits elsewhere:""",
  """while an independent SVD-based
implementation (statsmodels) is nominal-compatible (.067, Wilson
[.034,.126]) on the same data and agrees with ours to three digits
elsewhere:"""),
 # ---------- §4.1: .067 精度 ----------
 ("""a standard SVD-based HC3-Wald
implementation (statsmodels) is calibrated on the very tasks where our
pinv-based path rejects at 1.00 (.067 over 60 seeds each), and the two
implementations agree to ≥3 significant digits on well-conditioned
tasks.""",
  """a standard SVD-based HC3-Wald
implementation (statsmodels) is nominal-compatible on the very tasks
where our pinv-based path rejects at 1.00 (8/120 rejections = .067,
Wilson [.034,.126] over 60 seeds × 2 tasks — sufficient to reject the
old rate of 1.00, not precision evidence of exact 5% calibration), and
the two implementations agree to ≥3 significant digits on
well-conditioned tasks."""),
 # ---------- Q3: 整段扩展(D 覆盖+精确调用+实测 cond+无滚动检查声明) ----------
 ("""**Q3. Does HC3 fix false alarms under heteroskedasticity — and was the
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
be read as path-specific until then.""",
  """**Q3. Does HC3 fix false alarms under heteroskedasticity — and was the
observed failure HC3's or ours?** Heteroskedasticity was never the
disease (Table 6). Table 2 shows *our* HC3-Wald path breaking on
near-collinear built designs; the original draft attributed this to the
correction itself ("the robust estimator is the failure mode"). The
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
lesson, not evidence against HC3; whether any residual HC3
miscalibration survives elsewhere is open, and C/D rows in Tables 2, 3,
and 6 should be read as path-specific outside the verified cells."""),
 # ---------- §4.2: Table 3 后段改写 + 配对 bootstrap ----------
 ("""At deployment-relevant sample sizes (n ≤ 200), A delivers roughly twice
B's calibrated power; at n = 1000, D overtakes (.65). Even so, no method
exceeds .65 — most of S2's H1 mass is hard for every test, which Section
4.3 decomposes by family.""",
  """At deployment-relevant sample sizes (n ≤ 200), A delivers roughly twice
B's calibrated power; at n = 1000, D shows the highest pooled power
(.65). Pooled values, however, are means over only 19 H1 equations, and
task-level paired bootstrap (resampling equations; 10,000 replicates at
the Table 3 thresholds) shows the D-at-n=1000 advantage is not
established cross-equationally: Δ(D−A) = −.224 [−.418, −.016] at n = 50
(A significantly better), −.147 [−.426, +.133] at n = 200, +.060
[−.218, +.333] at n = 1000 — every interval either favors A or crosses
zero, and D beats A on 10–47% of tasks. Δ(D−B) = −.008/+.086/+.238
([-​.213,+.185]/[−.208,+.372]/[−.076,+.534]: all crossing zero);
Δ(D−C) = +.136 [+.033,+.250] at n = 50 and +.252 [+.025,+.478] at
n = 1000 (C being the retracted path, this compares within the
path-specific family). The one paired difference that excludes zero in
favor of a rival-family gate is Δ(F−A) = +.177 [+.054,+.323] at n = 50
and +.176 [+.056,+.320] at n = 200 (F wins 63–68% of tasks; F cannot be
paired with D without integrating D into the S5 harness — outstanding).
The defensible reading: D has the highest *average* calibrated power in
the current task pool at n = 1000, not established cross-equation
superiority; τ thresholds describe decision boundaries and are no
substitute for these intervals. No method exceeds .65 — most of S2's H1
mass is hard for every test, which Section 4.3 decomposes by family."""),
 # ---------- §6: HC3 单次复核段更新 ----------
 ("""**HC3 single-replication.** The independent statsmodels check covers the
fitted-HC3 variant on two failing tasks (60 seeds) and spot checks
elsewhere; a full re-implementation sweep across C/D, all suites, and
the calibrated-power rows is outstanding.""",
  """**HC3 independent verification.** After two revision rounds: fitted-HC3
(C) is verified on 2 failing tasks × 60 seeds; exog-HC3 (D) is verified
across all 13 H0 equations with stratified sampling (176 fits) — in
both cases our SVD path matches statsmodels to ≤ .0023 in p and ≤ 3×10⁻⁴
in covariance, while the X′X path fails only at measured
cond(X′X) ≥ 4×10¹⁷. Outstanding: the same sweep across the remaining
suites, noise layers, and the calibrated-power rows."""),
 # ---------- §7: 追加第二轮脚本 ----------
 ("""`paper/revision_analysis.py` regenerates every number in §4.1's
task-level paragraph, §4.6's set-level errors, §4.7's OR-budget and
transfer analyses, and §4.8's Monte-Carlo spreads from the released
per-cell parquet/CSV files; the independent HC3 check (Q3) uses
statsmodels 0.14 on identical data frames.""",
  """`paper/revision_analysis.py` regenerates every number in §4.1's
task-level paragraph, §4.6's set-level errors, §4.7's OR-budget and
transfer analyses, and §4.8's Monte-Carlo spreads from the released
per-cell parquet/CSV files; `paper/verify_round2.py` regenerates the
second-round D-verification, Wilson and paired-bootstrap numbers; the
independent HC3 check (Q3) uses statsmodels 0.14 on identical data
frames."""),
 # ---------- 附录 D: 第 8 条 ----------
 ("""(h) Related work extended to the wider
   specification-testing literature; references [15]–[22] split with
   access dates and per-source licenses; [8] description corrected;
   Table numbering restored to consecutive (v3 left a gap at 5).""",
  """(h) Related work extended to the wider
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
   Appendix C). (d) Missing glyphs (∈, ∅, ∉) fixed via symbol-font
   fallback. (e) S2 counts unchanged (32 tasks: 13 H0 + 19 H1; no task
   was removed; 1,480 cells and 508,825 datasets stand as audited)."""),
]
missed = []
for old, new in subs:
    if old not in md:
        missed.append(old[:60])
    else:
        md = md.replace(old, new)
open(p, 'w', encoding='utf-8').write(md)
print('missed:', missed or 'NONE')
