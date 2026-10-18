# Public artifact notes - 2026-10-18

The public repository and its data release postdate the frozen R2 manuscript.
README.md and data_release/manifest.json are the current artifact access record.
The original Git history and earlier paper/report snapshots are retained for provenance.

The publication copies of the three R2 audit scripts resolve this checkout's root
from their file locations. They were rerun on the local caches and reproduce the
audited numerical results. This portability change adds no new model-fitting data.

The frozen manuscript's reference [17] describes NIST WebBook as public domain.
That licensing description is incorrect: WebBook is copyrighted Standard Reference
Database 69. The public data bundle therefore excludes WebBook raw pages and
extracted tables. Rdatasets and ambiguous-license OpenML tables are likewise
source-download-only. Per-dataset decisions in data_release/dataset_sources.json
take precedence over historical blanket license descriptions.

The inherited root README was DimReflex framework documentation; it is preserved
in docs/DIMREFLEX_LEGACY_README.md. The current README describes this benchmark.
The package metadata and some historical builders still contain inherited settings
or machine-specific paths; the supported publication entry points are the README's
requirements, data-restoration command, and paper/research_R2 audit scripts.
