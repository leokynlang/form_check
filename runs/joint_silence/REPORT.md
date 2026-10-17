# Joint Silence: P(beautiful certificate | wrong form)

> ## ⚠️ 口径警告（2026-10-17 批注，引用本文件前必读）
>
> **本文件 Headline numbers 一节的数字与已保存原始产物不一致，不得直接引用。**
>
> 1. **SRSD 联合沉默率**：本文件记为 **32.0% / 24.6% / 19.4%**，并称覆盖
>    "19 non-monomial equations (45,000 datasets per n)"。而已保存的
>    `joint_silence_srsd_per_task.csv` 实际为 **18 个任务、每任务 2,000 份
>    数据集（每档 36,000）**，按样本量加权所得为 **37.8% / 29.8% / 23.2%**
>    （校准后 38.6% / 30.2% / 23.6%）；`summary.json` 的 `srsd_pool` 与后者
>    完全一致。**两者不可互替。**
> 2. **H0 基线**：本文件第 54–61 行用于计算后验的 "on true power laws (H0),
>    joint silence is 70%"，与 `toolbox_silence_h0_reference.csv` 记录的
>    **90.5%**（srsd, n=50, calibrated）不一致。该文件在 n=200 为 91.4%、
>    synth n=50 为 89.3%。
> 3. **因此本文件的 "P(wrong form | all-clear) ≈ 32/(32+70) ≈ 31%" 是用
>    两个不同口径的数字相除得到的。** 若改用同口径重算：完整工具箱的
>    H1 沉默率 42.2%（calibrated）与 H0 沉默率 90.5% 均来自可核实产物，
>    得 `42.2/(42.2+90.5) ≈ 31.8%` —— 结论方向偶然一致，但依据必须替换。
> 4. **条件代价亦不符**：第 101–105 行记 n=50 的 `P(err > 10%) = 8.0%`，
>    而 `summary.json` 的 `conditional_extrapolation` 为 **7.1%**（alpha）
>    与 **7.3%**（calibrated）。p90 亦为 8.6%–8.7% 而非 9.0%。量级接近，
>    但引用时须以 JSON 原始值为准。
> 5. **0.668% 与 6.874% 不在本文件**：本文件未出现这两个比例。其原始出处为
>    `paper/draft_v3_complete.md` 第 43–50 行（"0.67% of screened wrong forms
>    (481/71,992), but 6.9% of dangerous ones"）；`350/12,369 = 2.83%` 的出处为
>    `paper/review3_results.md` 第 38 行（10% 危险阈值）。引用这些比例时**不要
>    标注为本文件来源**。
> 6. **正确处理**：本文件保留作历史记录，不删除、不改动下方原文。正式引用
>    一律采用原始计数重算值，并以 `joint_silence_srsd_per_task.csv`、
>    `toolbox_silence.csv`、`toolbox_silence_h0_reference.csv`、
>    `summary.json` 为准。32.0% 与 45,000 份数据集的记载出处尚未查明，
>    可能是另一次运行或早期版本，须在复现阶段定位。
>
> 该批注已同步写入 D:/dimreflex 计划书 v1.4 第五章第 1 节及其修订记录。
> **在重新绑定并重跑之前，本文件的汇总数字不得进入论文表格或正文。**

**New analysis, 2026-09-09.** Question posed by the author's reframing:
not "which gate is better" but *when the fitted functional form is wrong,
how often does every falsification check stay quiet — so the pipeline
certifies a misfit model as if nothing were amiss?*

Scripts: `joint_silence.py` (four spec-test gates, S1+S2 data),
`joint_silence_toolbox.py` (adds the model-comparison gates, gate_cmp data).
Outputs: `runs/joint_silence/`. All numbers reproducible; preregistered
event definition: **joint silence = no gate rejects at threshold AND no
gate declined evaluation** (a declined fit is an honest refusal, not a
certificate; counting it as silence would flatter the gates).

## Definitions

- **Wrong-form pool**: structural H1 cells only (synthetic: 11 scenarios ×
  rel ∈ {0.5%,1%,3%,5%}; SRSD: 19 non-monomial equations). Noise-layer H1s
  are wrong-error-model, not wrong-form, and are excluded (§3.3 semantics).
- **Thresholds**: α = .05 nominal, and H0-pool-calibrated (5th-percentile
  p per gate) for a matched 5% false-positive budget.
- Gates: A = form_check permutation; B = RESET fitted-F; C/D = RESET
  HC3 variants on the SVD paths (v5 configuration); MC = CV-quad and
  AICc-quad (calibrated d\*/m\* from gate_cmp H0 pools).

## Headline numbers

### 1. Synthetic bench (S1): the toolbox almost never stays silent

| n | joint silence (4 gates, calibrated) |
|---|---|
| 50 | 0.15–0.22% |
| 200 | 0% |
| 1000 | 0% |

tpr_any ≈ .998 at n=50 and 1.000 above. **The residual silence concentrates
exactly where the paper already locates it**: OMITTED_M at rel ≤ 1% noise
(silence 1.6–2.4% — an omitted variable entering quadratically through a
correlated channel is genuinely hard at low noise) and CUBIC_LIN at rel=5%
(0.4%). At matched budgets the union of the four gates nearly closes on the
synthetic gallery: a beautiful-but-wrong certificate is a ~1-in-500 event at
worst, and does not occur at n ≥ 200.

### 2. SRSD (S2, SVD paths): real equations silence the whole panel far more often

Pool over 19 non-monomial equations (45,000 datasets per n):

| n | joint silence (α) | joint silence (calibrated) |
|---|---|---|
| 50 | **32.0%** | 32.6% |
| 200 | 24.6% | 25.0% |
| 1000 | 19.4% | 19.7% |

**Baseline for reading these: on true power laws (H0), joint silence is
70%** (calibrated; per-gate FPR 5% each, moderately overlapping → ~30%
false-alarm rate for the panel). So among wrong-form fits, the "all-clear"
rate is 32% vs 70% on correct forms — the panel *is* informative, but the
likelihood ratio is only ~2.4×: a silent certificate is weak evidence.
Equivalently: **P(wrong form | all-clear) ≈ 32/(32+70) ≈ 31% under the
bench's 50/50 class mix** — the posterior that a certified fit is
nevertheless wrong is barely reduced from its base rate.

**Every fourth wrong function fit to real physics data at n=50 receives a
clean certificate from all four spec tests simultaneously**, and the rate
decays slowly — still ~20% at n=1000. Per-task worst cases at n=50:
i.41.16 (78.7% silence), iii.4.33 (68.9%), ii.11.28 (66.3%), iii.4.32
(66.1%). This is the quantitative form of "we missed the real function all
the time": the panel's marginal powers (A .44, B .20, C .15, D .51 at n=50)
combine into a *coherent* wrong-form acceptance a third of the time.

### 3. Full toolbox (spec tests ∪ model comparison, gate_cmp coordinates)

Adding CV-quad and AICc-quad (the S5 gates that beat A head-to-head):

| bench | n | silence (calibrated) |
|---|---|---|
| synthetic | 50 | 6.3% |
| synthetic | 200 | 5.9% |
| SRSD | 50 | **42.2%** |
| SRSD | 200 | 30.5% |

Two readings, both honest and both uncomfortable:

- On **synthetic** data the union drives wrong-form silence to ~6% — but the
  H0 reference shows the union's calibrated silence on *true* power laws is
  89–92% (vs nominal 95%): the gates overlap, and the remaining 8% of true
  forms are flagged. On the synthetic gallery the toolbox is nearly closed.
- On **real-equation** data the number *rises* relative to the four-gate
  panel (32→42% at n=50). The gate_cmp bench is n ∈ {50,200} × rel=3% only
  (a harder noise regime than the S2 pool average), so the levels are not
  directly comparable across benches — but within the same coordinates,
  AICc-quad's higher marginal TPR does **not** seal the joint leaks: tasks
  where curvature is invisible to even-direction bases (ii.11.27: 88.5%
  silence at n=50; CUBIC_LIN: 53–66%) stay silent to every fixed gate.

### 4. Conditional cost: what does a silent certificate cost downstream?

Among silent-certified wrong-form fits on SRSD, the far-domain (one decade
out) extrapolation error distribution:

| n | median | p90 | P(err > 10%) |
|---|---|---|---|
| 50 | 0.8% | 9.0% | 8.0% |
| 200 | 0.3% | 2.0% | 1.0% |
| 1000 | 0.1% | 0.7% | 0.3% |

The median silent fit is *fine* — that is exactly why it is silent. The
risk is in the tail: ~8% of silent certificates at n=50 hide >10%
far-domain errors. Per-task "dangerous silence" (median far-error of the
silent fits > 2%): **i.27.6 (73.5% median far-error among its silent
fits, 16% of fits silent), i.37.4 (50.3%), i.39.11/i.43.43 (98.7% — the
known extrapolation-blind tasks; silence 2.7%), bonus.15/bonus.5
(8–13%)**. i.27.6 is the sharpest new exhibit: one fit in six is
silently certified and its median far-domain error is 73%.

## Framing (for the paper)

The benchmark's existing claims survive unchanged; this analysis adds the
*joint* statistic none of the per-gate tables compute:

> On real-equation tasks, the probability that an incorrect functional form
> receives a certificate — no specification test objects, no model
> comparison objects, no evaluation declines — is 32% at n=50 and still
> 19% at n=1000. The silence is not uniform: it concentrates in tasks
> whose misfit direction is near-orthogonal to every fixed direction set
> tested (odd-order curvature, cross-terms absent from the rival family),
> and 8% of silent certificates at n=50 carry >10% far-domain
> extrapolation error. The practical reading: passage through the current
> falsification toolbox is evidence of consistency with the tested
> alternatives, not evidence of correctness.

This *strengthens* rather than competes with the repositioning already in
§5/§7: if the whole toolbox certifies wrong forms a third of the time on
real equations, no single-gate ranking (the S5 question) is the first-order
problem — the direction-set coverage of the panel as a whole is.

## Disclosures

- rel = 0 cells excluded (no certificates exist: fits are exact,
  fail-registered).
- feynman-i.30.5 withdrawn (v5 App B); SVD rerun replay-side runner_error
  rows excluded via sentinel_ok.
- HC3 gates C/D are the SVD paths (v5 re-run); pinv-path numbers are
  path-specific and were not used here.
- SRSD "conditional cost" medians exclude data_fail rows (no extrapolation
  recorded); counts reported.
- gate_cmp bench covers rel=3% only; silence levels are not comparable
  across benches with different noise mixes, only within coordinates.
- Observational benches (World Bank, engineering) excluded: no ground
  truth for "the form is wrong".
