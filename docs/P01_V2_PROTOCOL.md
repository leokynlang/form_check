# P0-1 v2 frozen protocol and execution

Protocol: `p01-independent-recalibration-v2`.
Protocol file: `experiments/release_validation/protocol_v2.json`.
SHA-256: `91a2d14094869ce1c53dcafc4600b31ffc5615460faef30083086bc5ff59b525`.

The protocol and source snapshots were frozen before new data generation.
This is independent-sample validation on the same known equations; fixed-domain
and unseen-structural-group confirmation remain separate P0-2 work.
V1 outputs remain byte-for-byte preserved at `runs/p01_pilot_v1/`.
V2 writes only to `runs/p01_v2/`.

## Statistical change from v1

V1 used a finite empirical reference CDF. All raw values above the reference
maximum could tie at a bounded score; a low-budget panel threshold could reach
that ceiling and disable every valid alarm.

V2 decodes the existing A-D raw score `-p` and computes
`-log10(max(p, float64.tiny))`. E/F/G retain their higher-is-more-anomalous raw
comparison scores. For each method, each independent screened-H0 reference set
defines its median and the first positive finite scale in this frozen order:

1. IQR / 1.3489795003921634, with linear empirical quantiles.
2. MAD * 1.482602218505602.
3. Population standard deviation (`ddof=0`).

The member score is `(transformed raw - median) / scale`; a panel uses its maximum.
There is no CDF, clipping, tail extrapolation fit, or random tie jitter. An empty
or constant reference is unusable and produces a typed policy decline. Invalid
p-values or nonfinite standardized scores likewise decline.

This removes the artificial reference-CDF ceiling. It does not remove intrinsic
ties: A still uses 999 permutations, and machine-underflow p-values are clamped
at the smallest positive normal float. Both limitations are counted explicitly.
The scores are empirical ranking statistics, not claims of standard-normal null
distributions or new analytic p-values.

## Data design

- Existing fixed task registry: 13 H0, 19 H1.
- n = 50 and 200; log-space noise SD = 0.03.
- Twenty independent reference/calibration pairs per stratum.
- Each pair: 32 reference and 64 calibration samples per H0 task.
- Shared independent test bank: 128 samples per task per stratum.
- Total: 58,112 generated datasets; 8,192 are independent test datasets.
- H1 is generated only for test, never for score normalization or thresholds.
- Seed entropy includes v2 identifier, new master seed, phase, calibration round,
  task identity, n, noise and replicate. Test has round `-1`.
- All methods share exactly the same data within each record.

The test bank is intentionally reused across calibration policies to isolate
calibration sensitivity. Twenty decisions on one test sample are not twenty
independent test samples.

Methods and candidate fits are unchanged from the v1 adapter. G remains the fixed
RBF residual-CV rival, with no outcome-based bandwidth or regularization tuning.

## Policy calibration

Policies: D, F, A-D panel, A/E/F panel and G.
Targets: 1%, 5% (primary), 10%.
Screen: finite training **log-space** R2 >= 0.99.
Budget event: rejection OR any required member decline, conditional on screened H0.
Choose the smallest observed finite calibration threshold whose strict
exceedances plus declines stay within `floor(budget*N_screened)`.
If declines alone exceed the budget, the policy is infeasible.

There is no H1 tuning, test-threshold retuning, optional stopping or selection of
a best calibration round. Every actual test H0 withholding rate is reported,
including the fraction of calibration rounds exceeding the nominal target.
Empirical calibration compliance is not a population error-control guarantee.

## Repeated-calibration uncertainty

The point estimate is the equal mean of round-specific ratios. Pooling counts
can change the estimand when denominators vary, so it is not substituted for
the registered mean. Primary estimates require all20 policy rounds feasible
and each ratio defined; any feasible-only sensitivity is labelled separately.

Four complementary descriptions are saved:

1. Round sample SD, min/max and empirical 2.5%-97.5% percentile **range**:
   sensitivity to new reference/calibration samples conditional on the fixed test
   bank. This range is not a confidence interval; with20 rounds its tails are
   estimated imprecisely.
2. Conditional task bootstrap: keep the20 fitted policies fixed and resample
   test tasks within H0/H1, retaining13/19 task counts.
3. Round-only bootstrap: resample20 calibration rounds with original test-task
   weights fixed; interval for the empirical round-mean estimator.
4. Crossed bootstrap: in each draw use one shared round-index vector and one
   shared class-stratified task-weight vector across all rounds and policies.
   Compute within-round ratios and then their equal mean.

All bootstrap procedures use2,000 draws and paired samples for method effects.
Task weights are shared across rounds; independently resampling tasks per round
would incorrectly wash out persistent between-equation variation.

These intervals condition on the identities of the13 H0 calibration equations,
the empirical20-round bank and realized shared test samples. They do not cover
uncertainty about selecting a different set of calibration equations, refit H0
task-level policies under new task populations, or include record-level test
resampling. They do not establish unseen-equation generalization.

## Outcome and preserved scope

The legacy v1 endpoint remains unchanged for this P0-1 comparison: first feature
at three times the observed training maximum, others at training log medians.
Danger means relative error >10% there. Unknown outcomes remain in operational
accounting and are excluded only from explicitly labelled risk denominators.
P0-2 will replace this with independently fixed external domains.

## Reproduction and audit

```powershell
python -m pytest tests/test_release_validation.py tests/test_release_validation_v2.py tests/test_release_uncertainty_v2.py -q -o cache_dir=runs/p01_test_cache
python -m experiments.release_validation.runner_v2 --out runs/p01_v2 --workers 8
python -m experiments.release_validation.analyze_v2 --run runs/p01_v2
python validation/audit_release_policy_v2.py --run runs/p01_v2
```

Windows spawned workers use one BLAS thread each. Worker count affects execution
only, not job identities or samples. Complete cell checkpoints resume with the
same frozen fingerprint; source or protocol drift rejects resumption. A complete
run is verified and not regenerated by the generator. Code changes require a new
protocol/run; v1 is never an output destination.

The independent auditor reconstructs transforms, thresholds, decisions, counts,
round statistics and all three bootstrap intervals without importing the new
calibration or analysis implementations. It verifies every v1 file hash before
and after v2 and checks v2 seeds and dataset hashes against v1. For older legacy
caches lacking input dataset hashes, it verifies frozen source inputs and seed
separation rather than claiming exhaustive dataset-content reconstruction.
