# Extension: From Single-Form Rejection to Family Viability (battery)

**Annotated extension v0 — to be merged into §5 (as Q5) with a forward pointer
from §6 Limitations. Status: probe-complete (200 seeds/cell), formal 500-seed
run and SRSD-side validation in progress; numbers below cite the probe runs
(`runs/battery_probe`, `runs/battery_probe_wild`) and will be refreshed from
`runs/battery_full` before any submission.**

📝 **Placement craft.** This is not a fifth method and not a new claim in the
abstract — it is the *constructive response* to the benchmark's own most
uncomfortable finding (every test evaluates one H0; a rejection says "not this
form" and stops). Reviewers of falsification-side benchmarks reliably ask "so
what should a pipeline *do* with the alarm?"; this section is the answer, and
its honest scope statement (probe-scale, single-feature) is what keeps it from
reading as scope creep. Merge as Q5 so the analysis section keeps its
one-question-one-paragraph rhythm.

## Q5. Can the same machinery answer "not a power law — then what?" (family viability)

The four tests of §3.2 share one property: each evaluates a *single* null form.
A rejection is actionable only as an escalation (§4.4's portfolio argument),
but the benchmark's own family-stratified results (§4.3) hint the machinery
could do more: each test is a curvature check *in some linearizing
coordinates*, and several physical families linearize in known coordinates —
power law (log x, log y), exponential (x, log y), logarithmic (log x, y),
saturating rational (1/x, 1/y). We therefore evaluate a **family battery**:
for each family F, run the permutation curvature test of A in F's coordinates;
report the *viability set* V = {F : p_F ≥ α} — the families not rejected at
α = 0.05 — instead of a single verdict.

Two protocol points carry the epistemics:

**(a) The null distribution must be re-derived per family, and the classical
permutation null does not transfer.** A's residual-*permutation* null is
exact only under exchangeable errors in the transformed coordinates. Log
transforms homoskedasticize multiplicative noise — that is why the log-log
design was chosen — but the same noise is *heteroskedastic* in (x, log y) and
(1/x, 1/y) coordinates. The prediction is testable and was preregistered
(`battery_probe.py` docstring, written before the runs): under multiplicative
noise, the permutation null should hold for power/exponential coordinates
(v = log y) and break for logarithmic/rational coordinates (v = y, v = 1/y).
Measured at 200 seeds/cell on S1 scenarios: correct. The rational battery
member rejects its own true family — a data-generating process of
y = x/(0.2 + 0.1x) with multiplicative noise, i.e., the *best case* for the
rational test — at rates of **.34/.30** (rel = 1%/3%): a false-escalation
rate of 30%+ caused purely by null misspecification, not by any property of
rational laws. Replacing the permutation draws with a **wild bootstrap**
(Rademacher signs on residuals; Liu 1988; Davidson & Flachaire 2008) — a
one-line change to the null construction, no change to statistic or design —
drops the rational false-escalation to **.045/.090** while leaving the
power/exponential members' size unchanged (Table B1).

| Table B1 — battery member size/containment, 200 seeds/cell, n = 200 | perm null | wild null |
|---|---|---|
| H0 power truth: size of power member | .070/.075 | .075/.060 |
| H0 power truth: false escalation (V = ∅) | .070/.075 | .075/.060 |
| rational truth: containment (rational ∈ V) | **.66/.70** | **.955/.91** |
| rational truth: false escalation (V = ∅) | **.34/.30** | **.045/.09** |
| exponential truth: containment | .935/.945 | .92/.945 |
| logarithmic truth: containment | .87/.925 | .885/.95 |
| threshold truth (no battery family): escalation — *correct outcome* | 1.00/1.00 | 1.00/1.00 |

**(b) Viability is a non-rejection semantics, so multiplicity does not
inflate error.** V is "families each not rejected at their own α"; the
first-type error of interest (truth excluded from V) is a per-family property
and needs no cross-family correction. What multiplicity *does* buy, and we
report it as a finding rather than suppressing it: near-neighbor families
overlap. Under the logarithmic truth the modal outcome (70–79%) is
V = {logarithmic, rational} — the two saturating families are close
neighbours under a ×100 domain, and the battery reports both rather than
arbitrating. Arbitration belongs to the discovery layer (compare fitted
candidates), not to a specification gate; the battery's contract is
"never silently route a truth out of V", and its failure mode under (a) was
exactly that contract breach.

📝 **Findings craft.** The permutation→wild result is the section's
load-bearing sentence, and its structure deserves care: *the validity of a
null distribution is a coordinate-dependent property*. That is a general
lesson about specification batteries (any test re-used across transformed
coordinates inherits the transform's noise geometry), not a defect of A —
and it retroactively explains why the log-log choice was load-bearing all
along. State the two errors of the battery separately (false exclusion =
missed form, false inclusion = inconsistent shortlist) and give the cost
asymmetry in one sentence: exclusion costs an escalation to the expensive
layer (piecewise/auto); inclusion costs a downstream selection mistake.
The table keeps the benchmark's honesty conventions: the battery's worst
pre-fix number is bolded, and the threshold row is labeled *correct outcome*,
not a win.

## Scope and what this section does not claim

- **Probe scale.** 200 seeds/cell, four S1 scenarios, single-feature
  instances; the formal run (19 scenarios × 500 seeds, SRSD H1 pool) is in
  progress and this section's numbers will be refreshed; directional claims
  only until then.
- **No selection claim.** The battery returns V; it does not rank members of
  V. Model *selection* between viable families (AICc on the family-fitted
  candidates) is the discovery layer's job and is not benchmarked here.
- **The exponential blind spot is not fixed by the battery** — it was never a
  rejection-side blind spot. The §4.3 miss (exponential truths with time
  constants ≫ the sampled range produce near-linear log-log curves) is an
  *information* limit inside the training domain; the battery's exponential
  member has exactly the same in-domain information, so its power there
  matches A's, and dimensional normalization (x/T as a dimensionless group)
  remains the only lever that adds information. What the battery *does*
  remove is the routing error: an exponential truth now lands in
  V = {exponential} instead of triggering an uninformative "power law
  rejected" escalation.
- **Known synthetic-only blind spots carry over.** OPPOSING curvature
  (cancelling per-feature curvature) blinds the power member exactly as it
  blinds A (§4.2, Table 4); the formal run quantifies false *inclusion*
  there (wrong families remaining in V). A battery member for interaction
  structures (cross-terms in u) is straightforward to add but is *not*
  included in the frozen configuration — flagged as future work rather than
  retrofitted.

📝 **Scope craft.** Four bullets, each killing one anticipated reviewer move
before it's typed: "is this validated?" (probe scale, formal run in flight —
never ship a number and hide its seed count); "so it does model selection?"
(no — and the division of labor is the point); "didn't §4.3 say exponentials
are missed?" (yes, by *information*, not by routing — the one sentence that
prevents 'contradiction' reviews); "what about the OPPOSING blind spot?"
(carry-over, quantified, member for interactions declared as future work
*before* results are seen — the same preregistration discipline as P3).

## References to add (beyond the draft's list)

- Liu, R. Y. (1988). Bootstrap procedures under some non-i.i.d. models. *Annals of Statistics*, 16(4), 1696–1708.
- Davidson, A. C., & Flachaire, E. (2008). The wild bootstrap, tamed at last. *Journal of Econometrics*, 146(1), 162–169.

📝 **Citation craft.** Two citations, both for the wild bootstrap — the
battery's *mechanism* is A's own (mean-centered quadratic column, residual
null, pinv-row vectorization), so the section cites methods, not neighbors.
Adding a "specification batteries" citation would create the obligation to
differentiate from it; if such prior art surfaces in review, the honest
position is already in the text: this is a constructive extension of the
benchmark's own measurement axis, preregistered and probe-validated, not a
claimed-new test.
