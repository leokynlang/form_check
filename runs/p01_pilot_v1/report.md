# P0-1 independent-sample development pilot

Fresh mutually independent reference/calibration/test samples on the 32 already studied SRSD tasks. NOT unseen-equation confirmation. H1 is generated ONLY for test.

Reference and calibration use H0 only. Every method sees identical data.
Budget = rejection OR decline among screened H0; calibration and test are separate.
Test size is observed, not guaranteed equal to its target. Legacy endpoint only.
Task-bootstrap intervals freeze calibration; they omit calibration uncertainty.

| n | noise | target | policy | calibration withheld | test H0 withheld | H0 retained overall | dangerous H1 released / screened dangerous | unavailable outcomes |
|---|---|---|---|---|---|---|---|---|
| 50 | 0.03 | 1.00% | D | 0.26% | 0.13% | 92.31% | 43/90 (47.78%) | 64 |
| 50 | 0.03 | 1.00% | F | 0.39% | 1.04% | 91.47% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 1.00% | AD_panel | 0.65% | 0.91% | 91.59% | 26/90 (28.89%) | 64 |
| 50 | 0.03 | 1.00% | AEF_panel | 0.00% | 0.00% | 92.43% | 90/90 (100.00%) | 64 |
| 50 | 0.03 | 1.00% | G_rbf_cv | 0.78% | 1.04% | 91.47% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 5.00% | D | 4.95% | 5.46% | 87.38% | 22/90 (24.44%) | 64 |
| 50 | 0.03 | 5.00% | F | 4.43% | 2.99% | 89.66% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 5.00% | AD_panel | 4.82% | 5.46% | 87.38% | 14/90 (15.56%) | 64 |
| 50 | 0.03 | 5.00% | AEF_panel | 4.69% | 4.16% | 88.58% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 5.00% | G_rbf_cv | 4.95% | 4.68% | 88.10% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 10.00% | D | 9.90% | 9.10% | 84.01% | 21/90 (23.33%) | 64 |
| 50 | 0.03 | 10.00% | F | 9.64% | 8.32% | 84.74% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 10.00% | AD_panel | 9.64% | 8.32% | 84.74% | 9/90 (10.00%) | 64 |
| 50 | 0.03 | 10.00% | AEF_panel | 9.90% | 8.84% | 84.25% | 3/90 (3.33%) | 64 |
| 50 | 0.03 | 10.00% | G_rbf_cv | 9.64% | 10.01% | 83.17% | 3/90 (3.33%) | 64 |
| 200 | 0.03 | 1.00% | D | 0.78% | 1.17% | 91.23% | 6/77 (7.79%) | 64 |
| 200 | 0.03 | 1.00% | F | 0.65% | 0.65% | 91.71% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 1.00% | AD_panel | 0.39% | 0.78% | 91.59% | 4/77 (5.19%) | 64 |
| 200 | 0.03 | 1.00% | AEF_panel | 0.39% | 0.52% | 91.83% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 1.00% | G_rbf_cv | 0.00% | 0.13% | 92.19% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 5.00% | D | 4.69% | 6.38% | 86.42% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 5.00% | F | 4.95% | 3.91% | 88.70% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 5.00% | AD_panel | 4.95% | 5.08% | 87.62% | 4/77 (5.19%) | 64 |
| 200 | 0.03 | 5.00% | AEF_panel | 4.82% | 3.78% | 88.82% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 5.00% | G_rbf_cv | 4.95% | 3.78% | 88.82% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 10.00% | D | 9.90% | 8.85% | 84.13% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 10.00% | F | 9.51% | 8.46% | 84.50% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 10.00% | AD_panel | 9.90% | 11.46% | 81.73% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 10.00% | AEF_panel | 9.51% | 9.24% | 83.77% | 3/77 (3.90%) | 64 |
| 200 | 0.03 | 10.00% | G_rbf_cv | 9.90% | 8.98% | 84.01% | 2/77 (2.60%) | 64 |

## Interpretation limits

- The equations were previously studied. Fresh seeds give sample independence, not untouched-equation confirmation.
- A budget-compliant calibration fit can exceed its budget on independent test data; report that gap without retuning.
- Endpoint-unavailable records remain in operational metrics, with risk unlabelled.
- Fixed external domains, structural-group holdouts and mechanism confirmation remain pending.
- Finite reference CDFs can saturate: at n=50 and 1% target the AEF threshold reaches the score ceiling, disabling all valid alarms. This diagnoses the v1 score mapping, not inherent panel power. Do not tune v1 from its test outcomes; revise the protocol and use new samples.
- No winner is declared from this development pilot. Per-task counts and paired intervals are in summary.json.
