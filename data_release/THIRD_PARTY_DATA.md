# Third-party data in the public release

The repository's MIT license applies to repository-authored software and documentation.
It does not replace the original terms of third-party data. Generated numerical experiment
outputs are distributed with their configurations and provenance. Cite the manuscript and
the originating sources when using the data.

The release contains original snapshots and cleaned tables for the source families below.
The per-dataset source URLs, original hashes, transformations, and inclusion decisions are
recorded in dataset_sources.json and engineering_manifest.json. Cleaned tables were
modified by real_data/prepare_data.py (column selection, finite/positive-row filtering,
and the unit conversions documented in the manifest). Neither this project nor its authors
are endorsed by the source providers.

| Source | Release treatment | Credit and terms |
|---|---|---|
| UCI | Raw/clean tables included only when the dataset page confirms CC BY 4.0 | Credit each creator and dataset DOI in the retained UCI metadata; [UCI](https://archive.ics.uci.edu/), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| NIST StRD | Original certification files and cleaned tables | National Institute of Standards and Technology; [non-SRD data terms](https://www.nist.gov/open/license); supplied AS IS, with no warranty |
| USGS NWIS / FDSN | Original API snapshots and cleaned tables | U.S. Geological Survey, Department of the Interior; [copyrights and credits](https://www.usgs.gov/information-policies-and-instructions/copyrights-and-credits) |
| Open-Meteo / ERA5 | Original API snapshots and cleaned wind tables | Weather data by [Open-Meteo.com](https://open-meteo.com/), based on Copernicus ERA5 reanalysis; [CC BY 4.0 and source terms](https://open-meteo.com/en/licence) |
| DOE/EPA FuelEconomy.gov | Original download and cleaned table | U.S. Department of Energy / Environmental Protection Agency; [source data](https://www.fueleconomy.gov/feg/download.shtml) |

NIST non-SRD data notice: This data/work was created by employees of the National
Institute of Standards and Technology, an agency of the Federal Government. Pursuant
to title 17 United States Code Section 105, works of NIST employees are not subject to
copyright protection in the United States. This data/work may be subject to foreign
copyright. The data/work is provided by NIST as a public service and is expressly
provided AS IS. NIST makes no warranty concerning accuracy, fitness, or usefulness.
See the linked complete NIST terms for all disclaimers and attribution requirements.

## Source-download-only families

- NIST Chemistry WebBook / SRD 69: the site carries an explicit SRD copyright notice.
  Its raw pages and extracted tables are not mirrored in this release.
- Rdatasets: the mirror's maintainer states that the data rows' original licenses
  could not be determined definitively. Raw and cleaned tables are not mirrored.
- OpenML Ailerons/cpu_act/elusage: metadata marks them Public but does not supply
  an unambiguous standard redistribution license. They are not mirrored here.

These exclusions affect raw/clean tables, not this project's aggregate benchmark
measurements. The complete source registry and the existing fetch/prepare scripts
remain available. Obtain excluded snapshots from the originating source under its
terms; rerunning live downloads may not reproduce a historical snapshot byte for byte.

## Other supplied data

S2/S3 use configurations and formulas from [SRSD](https://github.com/omron-sinicx/srsd-benchmark);
the upstream code is MIT and the upstream datasets are CC BY 4.0. The release's run
records are this project's generated samples, predictions, and diagnostic outputs.
World Bank snapshots retain World Bank indicator and source metadata; see
[World Bank dataset terms](https://www.worldbank.org/en/about/legal/terms-of-use-for-datasets).
The root airfoil_self_noise.csv is UCI Airfoil Self-Noise (dataset 291), CC BY 4.0;
credit Thomas Brooks, D. Stuart Pope, and Michael Marcolini (NASA).
The heat-exchanger example is synthetic and is not an experimental measurement.

Terms checked for publication on 2026-10-18. Dataset availability is distinct from
a license to redistribute it; this file records the publication choices made here.
