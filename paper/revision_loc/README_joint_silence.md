# Joint silence (§4.9 in v5) — provenance

Question: when the fitted form is wrong, how often does the whole
falsification panel stay quiet — a "beautiful certificate" on a wrong
function? Added to the v5 lineage 2026-09-09 as new §4.9 and an abstract
sentence.

Scripts (this directory):
- `joint_silence_v5.py` — four-gate joint silence on S1 (runs/full/) and
  S2 (runs/cd_svd/s2_cd_svd.csv joined to runs/srsd_full/raw_srsd.parquet
  on task/n/rel/seed; sentinel_ok filter; feynman-i.30.5 excluded per
  App. B). Writes runs/joint_silence/.
- `joint_silence_toolbox_v5.py` — adds the S5 model-comparison gates
  (CV-quad, AICc-quad) on the same seeds (runs/gate_cmp/), thresholds from
  runs/gate_cmp/summary.json.

Outputs:
- runs/joint_silence/joint_silence_synthetic.csv (per scenario × n × rel)
- runs/joint_silence/joint_silence_srsd_per_task.csv / _pool.csv
- runs/joint_silence/conditional_extrapolation.csv
- runs/joint_silence/toolbox_silence.csv (+ _h0_reference.csv)
- runs/joint_silence/summary.json, REPORT.md (working notes, English)

Event definition (fixed before first run, recorded in the docstrings):
joint silence = no gate rejects at threshold AND no gate declined
evaluation (fail-registration is an honest refusal, not a certificate).
Two thresholds: nominal α = .05 and H0-pool-calibrated (5th-percentile p
per gate per bench × n). Levels agree to ~1 point everywhere; the paper
quotes α.

Disclosures baked in: rel = 0 cells emit no certificates (exact fits,
fail-registered) and are excluded; withdrawn task excluded; noise-layer
H1s excluded (wrong error model, not wrong form); C/D are SVD-path;
gate_cmp covers rel = 3% only, n ∈ {50, 200}; extrapolation columns are
median/far-field relative error at ×3 the training log-range boundary
(half a decade), defined in run_srsd.py; the H0 silence baseline (70%)
comes from the same pools and is in REPORT.md.
