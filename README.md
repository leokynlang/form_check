# form_check

A research benchmark for specification gates, diagnostic blind spots, and
extrapolation-related joint release in physical model assessment.

The study compares residual-permutation curvature (`form_check`), RESET variants,
HC3 inference, and nested-quadratic model comparisons on prespecified power-law
candidates. It records false alarms, detection, failed/declined evaluations, and
the models released by a fit screen plus diagnostic panel.

## Current manuscript and evidence

- [Research manuscript R3 (PDF)](paper/research_R3/Silent_Acceptance_research_R3_2026-10-18.pdf)
- [Editable manuscript R3](paper/research_R3/manuscript.md)
- [Revision evidence](paper/research_R3/evidence/)
- [Revision notes](paper/research_R3/REVISION_NOTES.txt)
- [Preserved manuscript R2](paper/research_R2/Silent_Acceptance_research_R2_2026-10-18.pdf)
- [Public artifact and licensing corrections](docs/PUBLICATION_NOTES.md)
- [Data snapshot](https://github.com/leokynlang/form_check/releases/tag/data-2026-10-18)

The current manuscript is a working research paper. R3 preserves the historical
R2 evidence (18 eligible H1 equations for S2 power and 17 with finite legacy
extrapolation outputs) and adds independently sampled, repeatedly calibrated
follow-ups. It keeps legacy endpoint outcomes separate from the newer frozen
external-region outcomes. Older reports and scripts are retained as history;
their counts, calibration claims, and goals may differ. Legacy filenames containing
`no_causal_link` do not establish statistical independence or a causal conclusion.

The PySR outputs are exploratory snapshots with a separate protocol. They do not
extend the main benchmark into a validated test of arbitrary search-selected
expressions. Consult their reports and the manuscript's scope.

The next validation phase is tracked in
[the P0-1 / P0-2 / P1-1 execution plan](docs/ACCEPTANCE_EXPERIMENT_PLAN.md).
Its first independent-sample development pilot is in
[`runs/p01_pilot_v1/report.md`](runs/p01_pilot_v1/report.md), with a frozen
protocol, source snapshots, per-task accounting and an independent audit.
It is separate from the R2 results and does not establish unseen-equation
generalization. Additional numerical-reference tests use
`requirements-validation.txt`.

The completed [v2 follow-up](runs/p01_v2/RUN_NOTES.md) removes the reference-CDF
tail ceiling and repeats independent reference/calibration sampling 20 times,
using 58,112 fresh datasets and a shared 8,192-sample independent test bank.
[The frozen protocol](docs/P01_V2_PROTOCOL.md) and
[independent audit](runs/p01_v2/audit.json) accompany the results; v1 is preserved.
The uncertainty analysis distinguishes calibration-round sensitivity, conditional
test-task bootstrap and crossed round/task bootstrap. These results remain
conditional on the known equation pool and legacy endpoint outcome.

The [P0-2 execution plan](docs/P02_EXECUTION_PLAN.md) adds frozen valid external
regions and whole structural-group allocation. The completed
[full-grid results](runs/p02_fixed_domain/complete_grid/RUN_NOTES.md) cover
n=50/200/1000 and log-noise SD=0.005/0.01/0.03/0.05, with 246,240 fresh datasets,
108,000 independent test records, ten calibration rounds per stratum and paired
group/round uncertainty. The [first stage](runs/p02_fixed_domain/stage1/RUN_NOTES.md)
and [second-stage audit](runs/p02_fixed_domain/stage2/independent_audit.json)
remain separate; the [combined audit](runs/p02_fixed_domain/complete_grid/independent_audit.json)
checks their preservation, disjointness and exact aggregation. The main
screened-dangerous evidence comes from designed controls; all twelve physical
source strata have an undefined conditional primary risk because the relevant
denominator is zero. V1, V2, stage1 and R2 are preserved.

The [physical-coverage expansion protocol](docs/PHYSICAL_COVERAGE_PROTOCOL.md)
registers all 120 source equations, uniform domain and inclusion/exclusion rules,
conservative structural groups and an exposure ledger. It is a metadata/rule
registration, not a completed execution freeze or new physical confirmation.
Eighty-six equations meet its mathematical domain criteria; physical source
contracts, numerical/high-dimensional adapters and execution registration remain
pending. Under its conservative grouping, this catalog alone supplies no
uncontested unexposed eligible H1 group. A further whole physical source catalog
must be registered before claiming new multi-mechanism physical confirmation.

## Clone, install, and restore the data

```bash
git clone https://github.com/leokynlang/form_check.git
cd form_check
python -m pip install -r requirements-benchmark.txt
python tools/download_data.py --list
python tools/download_data.py
```

Code, reports, the manuscript, and small evidence files are versioned in Git.
Large recorded data are versioned as ZIP assets on the same repository's release.
The downloader verifies the archive SHA-256 and every extracted file, preserves
existing matching files, and refuses to overwrite different local data unless
explicitly passed `--overwrite`.

Select a bundle when only part of the record is needed:

```bash
python tools/download_data.py --bundle experiment-records
python tools/download_data.py --bundle engineering-data
python tools/download_data.py --bundle pysr-search
python tools/download_data.py --verify-only --bundle experiment-records
```

The [manifest](data_release/manifest.json) records the snapshot time, exact member
paths, sizes, and hashes. It includes historical runs and pilot outputs; their
presence does not imply that they belong to the final paper's analysis pool.

## Reproduce the paper's aggregation

R3's follow-up records are separate from the historical release snapshot above.
With the local P0-1/P0-2 records present, the complete-grid report verifies its
inputs and refuses to replace a changed completed report:

```bash
python -m pip install -r requirements-validation.txt
python -m experiments.reporting.p02_complete --run runs/p02_fixed_domain
python validation/audit_p02_complete.py --run runs/p02_fixed_domain
python -m experiments.physical_coverage_v1.registration --out runs/p03_physical_protocol --verify
```

The physical registration command only verifies source metadata and frozen rule
snapshots. It does not run the proposed physical-coverage experiment. See the
execution plan for the original frozen experiment commands and the protocol for
the remaining source-contract and implementation work.

After restoring `experiment-records`:

```bash
python paper/research_R2/audit_size.py
python paper/research_R2/audit_results.py
python paper/research_R2/audit_s7.py
python paper/research_R2/reproduce_r2.py
python paper/research_R2/verify_composition.py
python paper/research_R2/threshold_decomposition.py
python paper/research_R2/deep_analysis.py
python paper/research_R2/audit_task_intervals.py
```

The first three rebuild the aggregated audit evidence behind the manuscript's
tables. The remaining five produce the task-level composition and uncertainty
analysis added on 2026-10-18: `reproduce_r2.py` re-derives the joint-silence
counts, `verify_composition.py` and `threshold_decomposition.py` decompose the
released records by task across danger thresholds, `deep_analysis.py` reports
the far-threshold sensitivity of the same counts together with the per-gate
contribution of the PySR counterpart experiment, and `audit_task_intervals.py`
resamples the 17 eligible equations to report task-bootstrap and
leave-one-task-out intervals for Ld/S, Ld/Sd, and Ld/Da.

All eight scripts resolve the repository root from their own locations. They read
the original caches and write JSON evidence and figures only under
`paper/research_R2/`. Re-aggregation is separate from regenerating the original
model-fitting experiments. `deep_analysis.py` additionally reads the PySR
counterpart caches under `runs/pysr_gate/wide_svd/`, which the same bundle
contains.

The PDF is supplied ready to read and is versioned with the manuscript in this
repository. Its optional ReportLab builder requires Windows fonts (Times New Roman,
Consolas, SimHei, and Segoe UI Symbol); the statistics scripts do not require those
fonts. R3's `verify_revision.py` checks evidence hashes, preserved R2 tables,
audit certificates and PDF bounds/glyphs, and requires a visual-review record
bound to the final PDF hash. Its PDF checks use the PDF-capable build runtime.

## Run experiments and tests

```bash
python run_experiment.py --tag pilot --seeds 50
python run_srsd.py --tag srsd_pilot
python analyze.py pilot
python -m pytest -q
```

The SRSD loader fetches the upstream configuration repository on first use.
Generating the full suites can be expensive; restoring the published caches is
sufficient for the aggregation commands above.

The vendored `dimreflex/` package supplies the original diagnostic implementation.
The inherited `pyproject.toml` describes the historical DimReflex package layout;
use the requirements file above for this benchmark checkout. Framework documentation
is preserved in [docs/DIMREFLEX_LEGACY_README.md](docs/DIMREFLEX_LEGACY_README.md).

## Data attribution and redistribution

See [THIRD_PARTY_DATA.md](data_release/THIRD_PARTY_DATA.md),
[dataset_sources.json](data_release/dataset_sources.json), and the
[engineering registry](data_release/engineering_manifest.json).

The release mirrors licensed UCI snapshots, NIST StRD, USGS, Open-Meteo/ERA5, and
DOE/EPA FuelEconomy.gov data with source credit and transformation notes. WebBook
SRD 69, Rdatasets, and OpenML tables with unresolved redistribution terms are
source-download-only; their registry and acquisition scripts remain available.
Their raw/clean tables are not included in the public ZIPs.

```bash
python real_data/fetch_data.py
python real_data/prepare_data.py
```

These commands contact the original services and are subject to their terms.
Live responses may differ from the historical snapshots.
The root heat-exchanger example is synthetic, not an experimental measurement.

## License

Repository-authored code and documentation retain the [MIT license](LICENSE).
Third-party data retain their own terms and attribution requirements.
