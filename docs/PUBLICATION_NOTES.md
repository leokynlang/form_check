# Public artifact notes - 2026-10-18

The public repository and its data release carry the R2 manuscript together with
the audited evidence behind it. README.md and data_release/manifest.json are the
current artifact access record. The original Git history and earlier paper/report
snapshots are retained for provenance.

The publication copies of the R2 scripts resolve this checkout's root from their
file locations. They were rerun on the local caches and reproduce the audited
numerical results. This portability change adds no new model-fitting data.

The manuscript received a task-level revision on 2026-10-18. Three passages were
added to Section 4.10 and a note to Table 15. They report that the released
dangerous records are concentrated in two tasks at every danger threshold, and
that six of the 17 pooled tasks contribute 30,351 dangerous records of which none
is released. No reported aggregate (S, Sd, Ld, Da), table row, or figure changed,
and no new experiment was run. Four scripts were added under paper/research_R2/
to reproduce the analysis. The PDF was rebuilt and renamed from
Silent_Acceptance_research_R2_2026-10-17.pdf to
Silent_Acceptance_research_R2_2026-10-18.pdf.

Reference [17] of the manuscript describes NIST WebBook as public domain. That
licensing description is incorrect: WebBook is copyrighted Standard Reference
Database 69. The public data bundle therefore excludes WebBook raw pages and
extracted tables. Rdatasets and ambiguous-license OpenML tables are likewise
source-download-only. Per-dataset decisions in data_release/dataset_sources.json
take precedence over historical blanket license descriptions.

The inherited root README was DimReflex framework documentation; it is preserved
in docs/DIMREFLEX_LEGACY_README.md. The current README describes this benchmark.
The package metadata and some historical builders still contain inherited settings
or machine-specific paths; the supported publication entry points are the README's
requirements, data-restoration command, and paper/research_R2 scripts.
