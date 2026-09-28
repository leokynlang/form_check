# Real engineering holdout validation — 2026-09-05

This is an exploratory, outcome-independent analysis plan, saved before computing
any model/test outcomes in this run. Dataset descriptions and source data previews
have already been read. This is not a registered confirmatory study. No dataset
will be dropped because results are unfavorable. No product code will be changed.

## Data and rationale

- NIST Chwirut1: 214 observed ultrasonic calibration measurements; response versus
  metal distance. A smooth departure from a power law may be detectable despite
  high training R-squared. Group by exact distance; unidentified experimental-run
  dependence cannot be removed because run identifiers are unavailable.
- NIST Nelson: 128 observed dielectric breakdown measurements versus time (weeks)
  and temperature (Celsius). Convert temperature to Kelvin before taking logs.
  Group by exact time-temperature combination. This tests new combinations, not
  independent specimens or laboratories; specimen identifiers are unavailable.
- Official NIST nonlinear model forms are reference comparators, NOT known true
  laws. Refit all coefficients using training observations only. Never use the
  certified full-data coefficients for prediction or initialization.

## Fixed evaluation

- Seeds 20260905 through 20260924 (20 resplits of the same two datasets, not 20
  independent datasets). Random-row holdout 25%; random-group holdout 25% groups.
- Chwirut1 domain holdout: highest 25% of unique distances. Also a narrow-domain
  stress test: middle 50% of unique distances for training, both tails for testing.
- Nelson domain holdouts: highest 25% of unique times; separately highest
  temperature level. Ties kept together. No target values determine splits.
- On training data only: unchanged native form_check (200 permutations), existing
  common.py empirical-p wrapper (999 permutations, alpha .05), three RESET
  variants. Native 3-sigma alarm and wrapper .05 decisions reported separately.
- Sensitivity check: separate experimental HC3-adjusted wild bootstrap for the
  same quadratic direction, not credited to the product form_check.
- Models: power law; log-space polynomial total degrees 2 and 3; degree 9 for
  Chwirut1 only as a deliberately flexible stress comparator; unchanged
  discover_auto (default families and backend); training-refitted NIST form.
- A training-only grouped 4-fold selector among polynomial degrees 1,2,3
  (plus 9 on Chwirut1). Minimum mean validation log-MSE, degree breaks ties.
  Compare with built-in discover_auto, whose internal CV is row-wise.
- Primary metric log-RMSE; also training/test log-R2 and median absolute relative
  error. No clipping of divergent predictions; nonfinite predictions and
  rank-deficient polynomial fits are recorded as failures.
- Record all splits, per-observation predictions, model failures, input/code
  SHA256s and environment. Summaries use median and IQR; no false-positive/power
  claims on real data with unknown conditional means. Resplit frequencies are
  stability summaries, not independent-trial significance tests.

## Interpretation

Evidence of usefulness: training-only warning plus lower held-out error after
changing model. Evidence of model-selection optimism: low training/internal-CV
error coupled with materially worse new-group or new-range error. An alarm alone
does not prove overfitting. Failure to alarm is not model certification. Classical
and permutation p-values can be affected by heteroscedasticity or dependence.
Only the form-check/model-selection portions are evaluated; this is not a full
automatic dimensional-analysis validation or proof of superiority to RESET.
