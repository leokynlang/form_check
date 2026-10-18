# R3 minimum acceptance route: P0-1 / P0-2 / P1-1

Status: 2026-10-18. R2 manuscript, old runs and original methods remain unchanged.
P0-1 v2 and both frozen P0-2 stages are complete; see the execution updates below.
The initial design and stage-specific pilot limitations remain recorded for
provenance. P1-1 has not been executed, and the new physical-coverage protocol is
a source/rule registration with execution obligations still pending.

## 1. Repository audit and reuse

| Asset | Reuse | Required addition / limitation |
|---|---|---|
| `runs/srsd_full/raw_srsd.parquet` | 240,000 rows / 32 equations; original A/B, seeds, legacy endpoint | Development/replay only; no E/F on these data; original C/D must not be used |
| `runs/srsd_full/quality_recompute.parquet` | 240,000 matching rows; training **log-space** R2 and diagnostics | Assert one-to-one keys before joining; preserve failed rows |
| `runs/cd_svd/s2_cd_svd.csv` | 240,000 rows; corrected C/D and replay sentinel | Use sentinel provenance; default `common.run_all_methods` still uses historical HC3 path |
| `runs/gate_cmp/results.csv` | 20,400 rows, 51 synth/SRSD tasks; A/E/F | Different seeds/grid from S2, no direct row join with S2; no cached C/D or fit screen |
| `srsd_data.py` and cached upstream SRSD classes | Equation registry, declared sampling, log-noise generator | Pin/hash upstream source; record data hashes; validate sampling semantics before confirmation |
| `common.form_check_with_p`, classical RESET | A statistic sentinel, B, typed failures | Explicitly invoke these; never silently fall back from corrected C/D |
| `paper/revision_loc/rerun_cd_svd.py::reset_test_svd` | Registered stable C/D implementation | Wrap as one adapter; numerical reference tests on new input data |
| `gate_compare.gate_model_compare` | Common full-quadratic E/F and deterministic five folds | Run on exactly the A-D datasets; finite-statistic checks per gate |
| `scenarios.py` | Mechanism development and smoke tests | Some truths depend on sample mean/std: cannot directly define a fixed external truth |
| R2 `audit_*.py`, threshold decomposition | Denominator conventions, task bootstrap, historical anchors | New independent score calibration, policy budget and evaluation; old audits are not confirmation |
| `runs/real_engineering_holdout` | Existing Chwirut1/Nelson held-out prediction workflow | Useful implementation reference; not known-truth specification validation |
| `real_data/` | Compatibility and source manifests | No truth labels: cannot estimate false alarms or dangerous release |

Do not increase old seed counts as a substitute for new independent equations.
There is no applicable AGENTS.md under this repository at inspection time.

## 2. P0-1 tasks (first implementation)

| ID | Executable work | Acceptance / output | Depends on |
|---|---|---|---|
| 01a | Inventory input sizes, schemas, hashes, key uniqueness and seed grids | `runs/p01_inventory/inventory.json`; all joins verified; no data mutation | None |
| 01b | Shared A-G adapter: A/B reuse, C/D stable, E/F reuse, G richer RBF residual-CV rival | Typed per-method failures and timing; all methods share dataset hash | 01a |
| 01c | Independent reference / calibration / test seed namespaces | Dataset and seed disjointness assertions; fixed protocol digest | 01b |
| 01d | Calibrate complete policies at 1/5/10% H0 withholding budgets after R2 screen | Include declines in budget; infeasible budgets recorded; no H1 used | 01c |
| 01e | Run fixed fresh-sample pilot on existing 32 SRSD tasks | 6,592 datasets for v1; known tasks, independent samples; H1 only in test | 01d |
| 01f | Evaluation and paired uncertainty | Policies, decisions, per-task counts, paired differences, trade-off plots; primary 5% | 01e |
| 01g | Independent held-out-equation confirmation | New grouping/split manifest, frozen method config and thresholds; actual held-out H0 size | P0-2 registry |

G is a **fixed richer rival**, not a proven omnibus test: standardize log features
using each CV training fold only; at most 32 training-fold RBF landmarks;
bandwidth squared = feature count; ridge = 0.01 per training observation;
unpenalized affine terms; score = CV RMSE(power) - CV RMSE(RBF rival).
No pilot-driven tuning. Primary panels remain AD and AEF, with G standalone.

Calibration uses three non-overlapping streams: reference H0 defines member CDFs,
calibration H0 defines the **joint** threshold, test H0/H1 estimates outcomes.
All phases retain the same registered fit screen. A missing member withholds
release; it is not an alarm, a success, or a row to drop. No threshold is selected
by maximizing H1 detection. Strict score exceedance makes ties conservative.
If declines alone exceed the budget, the policy is marked infeasible and is not
presented as satisfying that budget.

Report the following with explicit counts: H0 screened withholding, separate
rejection/decline, H0 overall retention, H1 detection yield, all/screened dangerous
release, screened-safe retention, total release, unavailable outcomes, and
dangerous-H1 share among all releases **only on finite-outcome records**.
The latter is a descriptive mixture-dependent risk, not a deployment posterior;
dangerous true-form H0 outcomes are counted separately as well.

Pilot inference resamples equations jointly for all methods, freezes thresholds,
and reports the number of denominator-contributing tasks. It does not estimate
uncertainty from choosing the calibration policy. Confirmation must add an outer
calibration-task resample / multiple independent calibration draws.

Suggested execution:

```powershell
python -m experiments.release_validation.inventory --out runs/p01_inventory
python -m pytest tests/test_release_validation.py -q
python -m experiments.release_validation.runner --out runs/p01_pilot_reproduce
python -m experiments.release_validation.analyze --run runs/p01_pilot_reproduce
```

Use a new directory for a new protocol. Resume accepts only identical protocol
and source hashes, writes atomic checkpoints, and never overwrites R2 outputs.
The pilot is diagnostic evidence, not a preregistered test on untouched tasks:
all 32 equations have already appeared in the research process.

The first completed run is `runs/p01_pilot_v1/`. Its exact generation sources are
preserved under `source_snapshot/` with a verified `source_snapshot.json` map.
Subsequent import cleanup and a retry for transient Windows metadata replacement
locks did not alter statistics or redraw samples. Source hashes intentionally
prevent resuming a frozen run with modified code; reproduce into a new directory.
Install `requirements-validation.txt` for the independent statsmodels checks.
Run `python validation/audit_release_policy.py --replay` for historical parity,
and `python validation/audit_release_policy.py --run runs/p01_pilot_v1` for saved
policy accounting checks (the current auditor's fixed row count targets v1).

### Current execution status

- 01a complete: S2 three-way keys unique and complete; no S2/S5 dataset-seed overlap.
- 01b-01f pilot complete: 6,592 datasets, 4,096 independent test samples, 30
  policy/budget/n cells, 61,440 saved policy decisions. Nine invariant/parity
  tests and an independent accounting audit pass. Historical replay checks 64
  records, with zero mismatches on sentinel-valid records.
- Important pilot finding: empirical CDF scores saturate in the extreme tail.
  At n=50 / 1% target, AEF's threshold hits the score ceiling and no valid
  observation can alarm. This is a v1 calibration design limitation, not evidence
  that AEF inherently cannot detect departures.
- Before 01g: replace the saturating score map with a frozen continuous
  H0-reference-only scaling or another explicitly validated panel combination;
  validate at the required tail resolution; use a new protocol ID and fresh
  samples. Never overwrite or optimize the v1 results.
- Also pending: multiple calibration draws / calibration uncertainty, n=1000
  and the full noise grid, the P0-2 structural-group task split and fixed regions.
- No confirmatory method ranking or new-equation generalization claim is made.

### V2 follow-up completed

The registered v2 follow-up is complete; see `docs/P01_V2_PROTOCOL.md` and
`runs/p01_v2/RUN_NOTES.md`. The v1 notes above remain historical descriptions.
V2 replaced the empirical CDF with H0-only robust affine scaling, generated
58,112 wholly fresh records, and repeated reference/calibration20 times with an
independent8,192-sample test bank. All600 policies were feasible. The independent
auditor verified2,457,600 decisions,390 policy metric summaries and90 paired
metric summaries, including all three bootstrap interval procedures. All150
v1 files remain byte-for-byte unchanged.

The n=50 /1% AEF ceiling pathology is resolved:20/20 rounds have nonzero H1
alarms. Calibration uncertainty is now explicitly measured, but task uncertainty
still dominates; the5% AEF dangerous releases remain concentrated in bonus.15.
At V2 completion, held-out structural groups, fixed domains, full noise/n coverage
and uncertainty from changing H0 calibration task identities remained pending.
The later P0-2 study addresses the first three within its registered task library;
calibration-design identities remain fixed, and the broader physical mechanism
coverage limitation remains unresolved.

## 3. P0-2 registered design and completion update

| ID | Work | Reuse / new | Completion criterion |
|---|---|---|---|
| 02a | Audit all 120 source formulas and exclusions; record equivalence/structure groups | Reuse SRSD; add explicit group registry | No near-equivalent equations leak across development/calibration/test |
| 02b | Build ~40-60 H1 equations and ~20-30 H0 designs across six mechanisms | New registry; distinguish genuinely new formulas from parameter variants | Domain, truth, noise, mechanism, physical/synthetic status fixed before outcomes |
| 02c | Freeze training and external domains before sampling | New per-task domain manifests | Same positive, physically valid evaluation coordinates for every n/seed |
| 02d | External grids: every feature direction, joint directions, interaction directions | Reuse equation callables, new grid builder | Grid hashes; explicit truth-near-zero and domain-invalid handling |
| 02e | Region metrics and sensitivity | New metrics, retain legacy endpoint | Primary weighted mean relative error; p95/max; 2/5/10/20% thresholds; missing outcomes explicit |
| 02f | Independent confirmation at n=50/200/1000, rel=.005/.01/.03/.05 | Reuse P0-1 runner and policies | Freeze sample count before confirmation; suggested start 100 seeds/cell |

H0/H1 new tasks are allocated into development/calibration/confirmation at the
**structural-group** level. Numeric parameter variants are not counted as
independent mechanisms. New task quantities are a design target, not a claimed
power calculation. Set the final size using development variability and compute
cost, then freeze it without inspecting confirmation outcomes.
Full Cartesian grids are not mandatory; task-specific weighted domains are.
Any invalid grid must be corrected before model outcomes are generated, and the
correction/version must remain in the manifest.

### Completed frozen execution

The finalized 110-task registry, eight H1 confirmation structure groups, four
H0 confirmation design groups and exact fixed regions are documented in
[P02_EXECUTION_PLAN.md](D:/form_check/docs/P02_EXECUTION_PLAN.md). Both stages
complete all twelve n/noise strata with 246,240 independent datasets and 108,000
test records. Each stratum has ten independent reference/calibration rounds.
[Complete results and audit links](D:/form_check/runs/p02_fixed_domain/complete_grid/RUN_NOTES.md)
retain strata separately; v1, v2, P0-2 stage1 and R2 remain unchanged.

This is structural confirmation on a combined source/control library, not broad
physical validation. All primary screened-dangerous confirmation cases come
from designed mechanisms. The physical-source conditional primary denominator
is zero in every stratum; its risk is undefined. The
[physical-coverage expansion registration](D:/form_check/docs/PHYSICAL_COVERAGE_PROTOCOL.md)
therefore uses an outcome-independent whole-source census, uniform admission and
domain rules, and conservative exposure groups. It has generated no new study
samples and does not yet supply multiple uncontested new physical mechanisms.

## 4. P1-1 tasks (mechanism prediction and repair)

| ID | Work | Reuse / new | Completion criterion |
|---|---|---|---|
| 11a | Explain bonus.5 / bonus.15 residual geometry and threshold proximity | Reuse cached decomposition and truth functions | Distinguish extrapolation amplification from cutoff-borderline behavior |
| 11b | Freeze two mechanism hypotheses: opposing curvature, interaction omission | Reuse synthetic prototypes; new population-defined truth generators | Predictions on new groups made before diagnostic outcomes |
| 11c | Compute coverage/projection measurements from independent design grids | New geometry module | No danger labels used to choose direction libraries or thresholds |
| 11d | Evaluate repair: separate squares / interaction augmentation | Reuse E/F, D; new ablation adapter only if needed | Compare with D/F and original panel at the same P0-1 joint budgets |
| 11e | Confirm on held-out mechanisms/parameter groups | New confirmation records | Report prediction failures as well as successes; paired task effects |

No need to invent a new method if F already repairs the missing directions.
A useful result can be that a simpler calibrated rival matches a larger panel.

## 5. Decision gates and manuscript changes

1. Gate A: pilot implementation validated; budgets actually enforced, no leakage.
2. Gate B: fixed-domain and held-out-equation results complete. If H0 budget
   fails to transfer, report the failure; do not silently recalibrate on test.
3. Gate C: mechanism prediction and repair evaluated under equal costs.
4. Rewrite claims around observed trade-offs and boundaries, not desired wins.

Immediate R2 edits are queued, not mixed into this development run: correct
Appendix H's claim of a failure counterpart; define training log-space R2;
scope 7.5% to its task mixture; add SRBench++ and richer specification testing
related work. After confirmation, move completed validation out of future work,
make policy trade-offs and fixed-domain results main figures, and keep old
endpoint/observational compatibility/audit history as supporting evidence.

Acceptance-relevant evidence is independently validated trade-offs plus a
cross-task mechanism/result. A positive panel win is not required. A result
confined to two already-known tasks must remain a case study. Never sample
until significance or until a dangerous release appears.
