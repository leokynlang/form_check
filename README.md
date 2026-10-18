# form_check

A research benchmark for specification gates, diagnostic blind spots, and
extrapolation-related joint release in physical model assessment.

The study compares residual-permutation curvature (`form_check`), RESET variants,
HC3 inference, and nested-quadratic model comparisons on prespecified power-law
candidates. It records false alarms, detection, failed/declined evaluations, and
the models released by a fit screen plus diagnostic panel.

## Current manuscript and evidence

- [Research manuscript R2 (PDF)](paper/research_R2/Silent_Acceptance_research_R2_2026-10-18.pdf)
- [Editable manuscript](paper/research_R2/manuscript.md)
- [Record-level evidence](paper/research_R2/evidence/)
- [Revision notes](paper/research_R2/REVISION_NOTES.txt)
- [Public artifact and licensing corrections](docs/PUBLICATION_NOTES.md)
- [Data snapshot](https://github.com/leokynlang/form_check/releases/tag/data-2026-10-18)

The current manuscript is a working research paper. Its numerical interpretation
uses the audited task pools: 18 eligible H1 equations for the S2 power comparison
and 17 with finite extrapolation outputs for joint release. Older reports and
scripts are retained as history; their counts, calibration claims, and goals may
differ from the R2 manuscript. In particular, legacy filenames containing
`no_causal_link` do not establish statistical independence or a causal conclusion.

The PySR outputs are exploratory snapshots with a separate protocol. They do not
extend the main benchmark into a validated test of arbitrary search-selected
expressions. Consult their reports and the manuscript's scope.

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

After restoring `experiment-records`:

```bash
python paper/research_R2/audit_size.py
python paper/research_R2/audit_results.py
python paper/research_R2/audit_s7.py
python paper/research_R2/reproduce_r2.py
python paper/research_R2/verify_composition.py
python paper/research_R2/threshold_decomposition.py
python paper/research_R2/deep_analysis.py
```

The first three rebuild the aggregated audit evidence behind the manuscript's
tables. The last four produce the task-level composition analysis added on
2026-10-18: `reproduce_r2.py` re-derives the joint-silence counts,
`verify_composition.py` and `threshold_decomposition.py` decompose the released
records by task across danger thresholds, and `deep_analysis.py` reports the
far-threshold sensitivity of the same counts together with the per-gate
contribution of the PySR counterpart experiment.

All seven scripts resolve the repository root from their own locations. They read
the original caches and write JSON evidence and figures only under
`paper/research_R2/`. Re-aggregation is separate from regenerating the original
model-fitting experiments. `deep_analysis.py` additionally reads the PySR
counterpart caches under `runs/pysr_gate/wide_svd/`, which the same bundle
contains.

The PDF is supplied ready to read and is versioned with the manuscript in this
repository. Its optional ReportLab builder requires Windows fonts (Times New Roman,
Consolas, SimHei, and Segoe UI Symbol); the statistics scripts do not require those
fonts.

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