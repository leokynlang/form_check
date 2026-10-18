# P0-1 v2: fresh data and repeated independent calibration

Known 32 equations; fully fresh data and20 independent reference/calibration pairs; shared independent test bank; legacy endpoint outcome.

20 independent H0 reference/calibration pairs; mean of round-specific ratios.
All targets include declines. Independent test error budgets are measured, not guaranteed.
The affine score has no reference-CDF ceiling; intrinsic p-value floors/ties remain.

| n | target | policy | feasible rounds | mean test H0 withheld | round sensitivity range | mean screened-dangerous H1 release | crossed interval |
|---|---|---|---|---|---|---|---|
| 50 | 1.00% | D | 20/20 | 0.43% | [0.19%, 1.28%] | 43.74% | [0.00%, 66.90%] |
| 50 | 1.00% | F | 20/20 | 0.99% | [0.39%, 1.69%] | 7.98% | [0.00%, 48.15%] |
| 50 | 1.00% | AD_panel | 20/20 | 0.58% | [0.23%, 1.05%] | 38.74% | [0.00%, 60.46%] |
| 50 | 1.00% | AEF_panel | 20/20 | 0.76% | [0.33%, 1.27%] | 7.98% | [0.00%, 48.15%] |
| 50 | 1.00% | G_rbf_cv | 20/20 | 1.03% | [0.45%, 1.99%] | 8.47% | [0.00%, 47.97%] |
| 50 | 5.00% | D | 20/20 | 4.12% | [3.08%, 5.89%] | 34.72% | [0.00%, 56.42%] |
| 50 | 5.00% | F | 20/20 | 4.68% | [3.54%, 6.18%] | 7.98% | [0.00%, 48.15%] |
| 50 | 5.00% | AD_panel | 20/20 | 4.40% | [3.67%, 5.35%] | 21.35% | [0.00%, 48.15%] |
| 50 | 5.00% | AEF_panel | 20/20 | 4.67% | [3.53%, 5.99%] | 7.98% | [0.00%, 48.15%] |
| 50 | 5.00% | G_rbf_cv | 20/20 | 5.46% | [4.32%, 7.33%] | 7.98% | [0.00%, 44.60%] |
| 50 | 10.00% | D | 20/20 | 9.27% | [7.57%, 10.75%] | 28.80% | [0.00%, 48.61%] |
| 50 | 10.00% | F | 20/20 | 10.02% | [7.74%, 12.72%] | 7.98% | [0.00%, 48.15%] |
| 50 | 10.00% | AD_panel | 20/20 | 9.57% | [7.79%, 11.04%] | 15.95% | [0.00%, 47.78%] |
| 50 | 10.00% | AEF_panel | 20/20 | 10.09% | [7.74%, 12.72%] | 7.36% | [0.00%, 44.44%] |
| 50 | 10.00% | G_rbf_cv | 20/20 | 10.22% | [9.20%, 12.08%] | 7.36% | [0.00%, 40.98%] |
| 200 | 1.00% | D | 20/20 | 1.20% | [0.85%, 1.77%] | 8.30% | [0.00%, 88.12%] |
| 200 | 1.00% | F | 20/20 | 0.85% | [0.49%, 1.47%] | 5.56% | [0.00%, 100.00%] |
| 200 | 1.00% | AD_panel | 20/20 | 1.47% | [0.87%, 2.42%] | 5.62% | [0.00%, 96.25%] |
| 200 | 1.00% | AEF_panel | 20/20 | 1.07% | [0.71%, 1.50%] | 5.56% | [0.00%, 100.00%] |
| 200 | 1.00% | G_rbf_cv | 20/20 | 1.40% | [0.81%, 2.38%] | 5.45% | [0.00%, 96.88%] |
| 200 | 5.00% | D | 20/20 | 5.06% | [3.67%, 6.22%] | 5.56% | [0.00%, 87.50%] |
| 200 | 5.00% | F | 20/20 | 4.79% | [3.54%, 5.93%] | 5.56% | [0.00%, 100.00%] |
| 200 | 5.00% | AD_panel | 20/20 | 4.69% | [3.38%, 6.23%] | 4.41% | [0.00%, 77.81%] |
| 200 | 5.00% | AEF_panel | 20/20 | 4.22% | [3.28%, 5.89%] | 5.42% | [0.00%, 95.62%] |
| 200 | 5.00% | G_rbf_cv | 20/20 | 5.36% | [4.33%, 6.94%] | 4.83% | [0.00%, 86.25%] |
| 200 | 10.00% | D | 20/20 | 9.79% | [8.26%, 11.69%] | 4.10% | [0.00%, 62.50%] |
| 200 | 10.00% | F | 20/20 | 10.59% | [9.28%, 12.44%] | 5.56% | [0.00%, 100.00%] |
| 200 | 10.00% | AD_panel | 20/20 | 10.38% | [8.01%, 12.89%] | 4.13% | [0.00%, 73.75%] |
| 200 | 10.00% | AEF_panel | 20/20 | 9.47% | [7.88%, 11.31%] | 3.54% | [0.00%, 63.12%] |
| 200 | 10.00% | G_rbf_cv | 20/20 | 10.45% | [8.64%, 12.15%] | 4.17% | [0.00%, 75.00%] |

## Uncertainty interpretation

- Calibration SD/min/max/percentile range measures sensitivity to independently redrawn reference/calibration samples on the same13 H0 equations; the range is not a CI.
- Conditional task intervals keep all20 fitted policies fixed and resample test equations stratified13 H0/19 H1.
- Round-only intervals resample the20 empirical calibration rounds with test tasks fixed.
- Crossed intervals jointly resample rounds and class-stratified test-task weights; the same vectors are paired across policies.
- All intervals condition on the original H0 calibration task identities and shared realized test bank. No new-equation guarantee or complete record-level Monte Carlo interval is claimed.
- Primary aggregates require every round feasible and its ratio defined. Feasible-only estimates are explicitly secondary.
- Legacy endpoint and known tasks remain; fixed domains and unseen structural-group confirmation are pending.
- V1 remains byte-for-byte unchanged; v2 uses no v1 sample seeds or dataset hashes.

## Calibration implementation diagnostics

AEF at n=50, 1% target has nonzero test alarms in 20/20 rounds.
This removes the v1 CDF-ceiling pathology, without claiming that all remaining p-value ties disappear.
