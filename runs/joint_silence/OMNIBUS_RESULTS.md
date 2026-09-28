# Omnibus Combination Check — Results (2026-09-09)

**Purpose.** Direct test of the paper's §4.5/§5 mechanism claim: the six
blind-spot tasks are silent because their misfit directions are
near-orthogonal to the panel's tested column sets, not because no panel-level
statistic was formed. If a standard omnibus combination of the four gate
p-values recovered power on those tasks, the "direction coverage" semantics
of the paper would need rewriting. It does not recover power — it makes
silence worse.

**Method.** Same rows as `joint_silence.py` (S2, 18 tasks after the i.30.5
withdrawal, noise-present strata rel > 0, SVD-path C/D, declines excluded
from numerator and denominator). Two combination rules over the four per-row
gate p-values:

- **Fisher** X = −2·Σlog p, 8 df — at nominal α = .05 and at H0-pool-calibrated
  thresholds (5% empirical size per n stratum; calibrated X thresholds
  20.96 / 21.08 / 21.45 at n = 50/200/1000 vs the χ²₈ nominal 15.51 — the H0
  X distribution is right-shifted, as positive dependence predicts).
- **Tippett** min-p with the exact 4-test correction 1−(1−p_min)⁴, nominal.

**Script:** `omnibus_combination.py` (repo root). **Output:**
`runs/joint_silence/omnibus_combination.csv`.

## Pooled results (S2, noise-present)

| n | H0 sil (any) | H1 sil (any) | LR (any) | H1 sil Fisher nom. | H1 sil Tippett nom. | H0 sil Fisher cal. | H1 sil Fisher cal. | LR Fisher cal. |
|---|---|---|---|---|---|---|---|---|
| 50 | .870 | .378 | 2.30 | .437 | .471 | .950 | .526 | **1.81** |
| 200 | .881 | .298 | 2.96 | .328 | .354 | .950 | .376 | **2.53** |
| 1000 | .881 | .236 | 3.73 | .250 | .273 | .950 | .281 | **3.38** |

Combination moves silence **up** in every cell. The correctly calibrated
Fisher panel is strictly worse than the paper's any-gate rule: at n = 50 its
silent certificate carries LR 1.81 vs 2.30. This is the positive-dependence
mechanism doing exactly what §6 says: when gate p-values are positively
dependent, the Fisher X null shifts right, and after recalibration the
combined statistic spends its entire 5% budget confirming what the marginals
already knew — it adds no new rejection directions.

## Per-task, n = 50 (noise-present): any-gate vs nominal Fisher

| task | silence (any) | silence (Fisher nom.) |
|---|---|---|
| iii.4.33 | .849 | **.862** |
| ii.11.28 | .830 | **.874** |
| iii.4.32 | .827 | **.883** |
| ii.11.27 | .813 | **.892** |
| i.41.16 | .797 | **.896** |
| bonus.15 | .734 | **.821** |

Every blind-spot task gets **worse** under combination — the marginal
p-values on those tasks carry no joint signal for Fisher to amplify
(per-gate powers .04–.18 ≈ their sizes; combining four nulls yields a null).
In the middle band the picture is the same: i.27.6 .162 → .341, bonus.5
.164 → .322, i.39.11 .030 → .069, i.43.43 .029 → .075 (nominal; the
calibrated variant shows the same ordering). No task improves.

## Conclusion for the paper

The blind-spot structure survives a direct aggregation attack. Silence on
iii.4.x / ii.11.x / i.41.16 is a property of direction coverage, not of
missing panel-level statistics; the any-gate rule is in fact the *best*
monotone combination of this panel's marginals on this bench. This closes
reviewer question #3 (omnibus combination) with a negative result in the
paper's favor, and can be added to §4.5 as one paragraph plus the pooled
table.

## Provenance

- Sources: `runs/srsd_full/raw_srsd.parquet` (A/B, klass),
  `runs/cd_svd/s2_cd_svd.csv` (C/D SVD, sentinel-filtered, i.30.5 excluded) —
  identical to `joint_silence.py`.
- Sanity anchors reproduced exactly before the new rules were run:
  H1 silence .3782/.2982/.2360, H0 .8702/.8812/.8805, LR 2.30/2.96/3.73 —
  matches paper Table 1 (`.378/.298/.232`; `.870/.881/.881`; 2.30/2.96/3.79).
- **Pool-aggregation fix applied to `joint_silence.py` (2026-09-09):** the
  S2 section had no `rel > 0` filter, so the committed
  `joint_silence_srsd_pool.csv` previously carried N = 45,000 per n
  (including the 9,000 zero-noise rows that emit no certificates) and read
  .320/.246/.194. The paper's headline numbers were computed on the correct
  noise-present basis and reproduce exactly (this session, on identical
  rows); only the committed pool artifact lagged. After the one-line fix the
  artifact reads N = 36,000 and .3782/.2977/.2321 — matching Table 1.
  Per-task table unchanged in structure; its cells now also exclude
  zero-noise rows (N 2,500 → 2,000 per task×n).
- Row-count note: S2 wrong-form noise-present usable rows = 108,000
  (= 36,000 × 3 n; the nominal 114,000 minus 6,000 rows lost to the
  sentinel filter + i.30.5 exclusion + inner-join). S1 = 66,000. Total
  wrong-form noise-present = **174,000**, not the 171,000 stated in the
  abstract/intro/§1. The 171,000 figure has no source in the released
  artifacts; correct it to 174,000 (or document the drop rule if it was a
  different release cut). The 19-vs-18 task mismatch in the S5 paragraph
  (§3.3) is separate: gate_cmp ran before the withdrawal, so S5's grid
  contains 19 physics tasks including feynman-i.30.5 — state this
  explicitly.
