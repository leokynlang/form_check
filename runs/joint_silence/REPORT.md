# Joint Silence: P(beautiful certificate | wrong form)

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
