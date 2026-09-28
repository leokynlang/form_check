# -*- coding: utf-8 -*-
"""Build the unified falsification-side paper (Silent Acceptance + HX) as a Word document."""
import os
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_ALIGN_VERTICAL
from docx.enum.section import WD_SECTION
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "Falsification_Side_of_Physical_Model_Discovery.docx")

doc = Document()

# ---------- page setup ----------
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
for a in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
    setattr(sec, a, Cm(2.5))

BLACK = RGBColor(0, 0, 0)
GRAY_HDR = "D9D9D9"

def _set_fonts(style, cn="宋体", en="Arial", size=12, bold=False, color=BLACK):
    f = style.font
    f.name = en
    f.size = Pt(size)
    f.bold = bold
    f.color.rgb = color
    rpr = style.element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.append(rfonts)
    rfonts.set(qn("w:eastAsia"), cn)
    rfonts.set(qn("w:ascii"), en)
    rfonts.set(qn("w:hAnsi"), en)

# ---------- styles ----------
st = doc.styles
_set_fonts(st["Normal"], size=12)
st["Normal"].paragraph_format.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
st["Normal"].paragraph_format.line_spacing = 1.5
st["Normal"].paragraph_format.first_line_indent = Pt(24)
st["Normal"].paragraph_format.space_before = Pt(0)
st["Normal"].paragraph_format.space_after = Pt(0)

h1 = st["Heading 1"]; _set_fonts(h1, cn="黑体", size=16, bold=True)
h1.paragraph_format.space_before = Pt(14); h1.paragraph_format.space_after = Pt(6)
h1.paragraph_format.line_spacing = 1.5; h1.paragraph_format.first_line_indent = Pt(0)
h2 = st["Heading 2"]; _set_fonts(h2, cn="黑体", size=14, bold=True)
h2.paragraph_format.space_before = Pt(11); h2.paragraph_format.space_after = Pt(5)
h2.paragraph_format.line_spacing = 1.5; h2.paragraph_format.first_line_indent = Pt(0)
h3 = st["Heading 3"]; _set_fonts(h3, cn="黑体", size=12, bold=True)
h3.paragraph_format.space_before = Pt(9); h3.paragraph_format.space_after = Pt(4)
h3.paragraph_format.line_spacing = 1.5; h3.paragraph_format.first_line_indent = Pt(0)

cap = st["Caption"]; _set_fonts(cap, size=10.5)
cap.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
cap.paragraph_format.first_line_indent = Pt(0)
cap.paragraph_format.space_before = Pt(4); cap.paragraph_format.space_after = Pt(6)

ref = st["Normal"]
# ---------- helpers ----------
def P(text, style=None, align=None, indent=True, size=None, bold=False, italic=False):
    p = doc.add_paragraph(style=style)
    r = p.add_run(text)
    if size: r.font.size = Pt(size)
    r.bold = bold; r.italic = italic
    if align is not None: p.paragraph_format.alignment = align
    if not indent: p.paragraph_format.first_line_indent = Pt(0)
    return p

def H(level, text):
    return doc.add_heading(text, level=level)

def CAPTION(text):
    return doc.add_paragraph(text, style="Caption")

def BULLET(text):
    p = doc.add_paragraph(style="List Bullet")
    p.add_run(text)
    p.paragraph_format.first_line_indent = Pt(0)
    p.paragraph_format.line_spacing = 1.5
    return p

def TABLE(headers, rows, widths=None, align_cols=None, header_bg=GRAY_HDR):
    t = doc.add_table(rows=1 + len(rows), cols=len(headers))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = False
    if widths:
        for i, w in enumerate(widths):
            for cell in t.columns[i].cells:
                cell.width = Cm(w)
    # header
    for j, htext in enumerate(headers):
        cell = t.cell(0, j)
        cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
        p = cell.paragraphs[0]
        p.text = ""
        r = p.add_run(htext)
        r.bold = True
        r.font.size = Pt(10.5)
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.first_line_indent = Pt(0)
        p.paragraph_format.line_spacing = 1.0
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear"); shd.set(qn("w:color"), "auto"); shd.set(qn("w:fill"), header_bg)
        cell._tc.get_or_add_tcPr().append(shd)
        # repeat header row on page breaks
        trPr = cell._tc.getparent().get_or_add_trPr()
        tblHeader = OxmlElement("w:tblHeader")
        trPr.append(tblHeader)
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = t.cell(i + 1, j)
            cell.vertical_alignment = WD_ALIGN_VERTICAL.CENTER
            p = cell.paragraphs[0]
            p.text = ""
            r = p.add_run(str(val))
            r.font.size = Pt(10.5)
            txt = str(val)
            if align_cols and align_cols[j] == "r":
                p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            elif align_cols and align_cols[j] == "c":
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            else:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if len(txt) > 12 else WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.first_line_indent = Pt(0)
            p.paragraph_format.line_spacing = 1.0
    return t

def PAGE_BREAK():
    doc.add_page_break()

def TOC_BLOCK():
    P("目录", size=16, bold=True, indent=False)
    p = doc.add_paragraph()
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), 'TOC \\o "1-2" \\h \\z \\u')
    run = OxmlElement("w:r")
    t = OxmlElement("w:t")
    t.text = "（在 Word 中按 F9 或右键“更新域”刷新目录与页码）"
    run.append(t)
    fld.append(run)
    p._p.append(fld)

# =================== CONTENT ===================
C = []

def h1(t): C.append(("h1", t))
def h2(t): C.append(("h2", t))
def h3(t): C.append(("h3", t))
def p(t): C.append(("p", t))
def b(t): C.append(("b", t))
def cap(t): C.append(("cap", t))
def tbl(headers, rows, widths=None, aligns=None): C.append(("tbl", (headers, rows, widths, aligns)))
def br(): C.append(("br", None))
def refitem(t): C.append(("ref", t))

# ---------------- Title block ----------------
h1("Falsifying Proposed Physical Laws: A Specification-Test Benchmark and a Hybrid Directional Test for Log-Linear Models")
p("Unified manuscript combining two companion studies — the falsification-side benchmark and the HX method — with full code and data evidence.")
p("")
p("Anonymous manuscript — working draft, not submitted.")
p("6 September 2026")
C.append(("toc", None))
br()

# ---------------- Abstract ----------------
h2("Abstract")
p("Symbolic regression and physical model discovery pipelines certify candidate models with goodness-of-fit metrics that measure how well a model fits — never whether an alternative model family would fit strictly better. A model can score near-perfectly on R², cross-validation, and information criteria while having the wrong functional form, and fail catastrophically at extrapolation. This paper treats the missing falsification side end-to-end, in two connected parts. Part I builds the falsification-side benchmark that the discovery community lacks: three task suites (1,157 grid cells, 394,260 fitted datasets plus 472,500 held-out evaluations) put four classical specification tests on one protocol — a common significance level, a common permutation null, power compared only at a matched empirical 5% false-positive rate, and explicit registration of every degenerate evaluation. The benchmark's findings are specific: heteroskedasticity-robust (HC3) RESET variants carry 42–56% false-positive rates on true physical monomials, flat in sample size; the classical-F RESET that dominates synthetic quadratics is nearly blind to rational laws (0.13–0.19); and no test wins everywhere — the failure profiles are complementary. A fourth suite (94,000 datasets) extends the machinery from rejection to family routing and shows that null validity is a coordinate-dependent property: a permutation null exact in log-log coordinates produces 30–34% false escalations in the rational family's (1/x, 1/y) coordinates, repaired to 4.5–9% by a wild bootstrap. Part II turns the benchmark's blind-spot map into a method: HX, a hybrid directional specification test that combines structural directions (quadratic, cubic, interaction, oscillatory) with fitted-value directions (powers of the restricted fitted response), removes redundant columns by rank-aware modified Gram–Schmidt, and evaluates both a familywise maxT statistic and a joint partial-F statistic combined by a Bonferroni minimum. In a post-development randomized stress test of 1,000 groups (393 null, 607 alternative), HX attains false-positive rate 0.0509 (95% Wilson CI 0.0332–0.0773) and power 0.9522 (0.9322–0.9665), versus 0.4942 for a single-direction form check, 0.8913 for a structural-only library, and 0.7512 for RESET — while the benchmark's size results establish that every one of these comparisons is meaningful only when read against each method's empirical null calibration. Both parts share a reproducibility discipline we make explicit: fixed seed rules, content-hashed configurations, typed failure registrations, and release of every generator, per-cell result, and analysis script. The code, the data, and the test evidence are catalogued in Appendix D.")
h2("Keywords")
p("model specification; log-linear physical laws; permutation inference; RESET; heteroskedasticity-robust inference; symbolic regression; falsification; wild bootstrap; maxT; dimensionless modeling")

# ---------------- 1 Introduction ----------------
h1("1 Introduction")
p("A regression model can be simultaneously excellent on every certification metric a discovery pipeline computes and wrong in its functional form. A power law fitted to Prandtl–Karman friction data reaches R² = 0.99, passes k-fold cross-validation, wins its information-criterion comparison, and reports small standard errors — while extrapolating with roughly −20% error one decade outside the training domain. The failure was reproduced as a controlled scenario in our laboratory notebook and is scenario S1 of the benchmark in Part I, where all four specification tests detect it at power ≥ 0.99. Fit quality and form correctness are different axes, and current discovery practice measures only the first.")
p("The statistics literature has owned the second axis for over fifty years. Ramsey's RESET [1] tests whether powers of the fitted values carry additional signal; built-regressor variants use powers of the regressors; White [2] and MacKinnon–White [3] provide heteroskedasticity-robust inference. Yet these tests are absent from symbolic regression (SR) benchmarking. SRBench [10] and its interpretability-aware successor SRBench++ [11] measure whether discovery methods find the right expression across hundreds of problems; SRSD [12] redesigns Feynman-derived tasks with realistic sampling ranges and discovery-oriented metrics. In neither framework — nor in any other discovery-side benchmark — does anything measure whether an attached falsification layer would reject a wrong form. The asymmetry has practical consequences: a pipeline with no working falsification layer silently accepts misfit forms, and its downstream use — extrapolation, engineering design, scientific claim-making — inherits the error silently.")
p("Measuring the falsification side faithfully requires resolving a protocol problem that, to our knowledge, no prior work addresses for this test family: the tests reject on incommensurable bases. A bootstrap-t threshold and an asymptotic F p-value measure different things; comparing them directly measures threshold placement, not statistical power. Our first pilot did exactly this and produced inverted rankings that vanished under a common basis. The protocol of Section 3 therefore rests on three commitments: one significance level for all tests; one permutation null shared by all tests; and power compared only at a common empirical false-positive rate estimated on true-monomial tasks.")
p("One question the benchmark's own results force into the open: every test evaluated in Part I answers a single-form question — \u201cis this a power law?\u201d — and a rejection terminates the inquiry. What should the pipeline do next? Two constructive answers follow from the benchmark's blind-spot map. First, a family-viability battery (Section 4.7) runs the same curvature statistic in the linearizing coordinates of four families and returns a viability set instead of a single verdict — and in doing so uncovers the coordinate-dependence of null validity. Second, and this is Part II of the paper, the map of complementary blind spots suggests a method: a hybrid test that augments the restricted model with both interpretable structural directions and omnibus fitted-value powers, protects the search with a maxT correction, evaluates joint signal with a partial-F statistic, and exposes every failure state instead of coercing ill-posed fits into a p-value.")
p("This paper contributes: (i) the falsification-side benchmark itself — four suites totalling 1,345 grid cells and 488,260 fitted datasets, with released generators under fixed seed rules that make any single cell independently reproducible; (ii) the protocol — matched-basis comparison with empirical-FPR calibration and explicit failure registration; (iii) the resulting failure-mode map — including two previously undocumented findings, the flat-in-n 42–56% false-positive rates of HC3-Wald RESET variants on real physical monomials and the near-blindness of fitted-F RESET to rational laws, plus the battery's coordinate-dependent null result; and (iv) HX, a hybrid directional test whose randomized stress-test record (size 0.0509, power 0.9522) is reported with exact denominators, uncertainty intervals, family-level results, and paired discordance counts — and whose honest limitations (frozen response-adaptive directions; a stress test, not external validation) are stated as precisely as its strengths.")

# ---------------- 2 Related Work ----------------
h1("2 Related Work")
h2("2.1 Specification testing")
p("RESET [1] augments the regression with powers of the fitted values and tests the joint significance of the added coefficients; the exog variant uses powers of the regressors instead. Heteroskedasticity-consistent covariance estimation [2] and its HC3 leverage-weighted refinement [3, 4] are the standard small-sample corrections and are defaults in several mainstream econometrics packages. These results are derived for linear models under fixed or quasi-fixed designs; none of the classical treatments contemplates log-transformed physical laws whose regressors span only a few orders of magnitude — precisely the regime our benchmark shows breaks HC3's small-sample behavior (Section 4.3).")
h2("2.2 Resampling nulls under non-i.i.d. errors")
p("Freedman and Lane proposed a residual-based permutation construction for general linear models [5]; Winkler et al. systematized permutation strategies, exchangeability blocks, nuisance adjustment, and error control [6]. The wild bootstrap [7, 9] was developed for inference under heteroskedasticity, and its permutation cousin requires exchangeability; Westfall and Young's maxT framework motivates simultaneous adjustment across correlated statistics [8]. The battery's mechanism finding (Section 4.7) is a compact new instance of the classical dichotomy: a log transform homoskedasticizes multiplicative noise, so residual permutation is valid in coordinates ending in log y and invalid in coordinates ending in y or 1/y — the null's validity inherits the noise geometry of the coordinates the test runs in.")
h2("2.3 Symbolic regression benchmarks and physics-informed discovery")
p("SRBench [10] established the discovery-side evaluation: 14 SR methods and 7 ML baselines on 252 problems, scored by exact recovery and held-out accuracy. SRBench++ [11] extends it in the interpretability direction: it incorporates domain-expert evaluations of interpretability on a domain-specific task and characterizes algorithm performance on distinct sub-tasks of regression — feature selection, avoidance of local minima, and others — across 12 modern SR algorithms. It is the most interpretability-aware discovery-side benchmark to date, and it still contains no rejection mechanism: an algorithm that returns a wrong expression and one that flags it score identically. SRSD [12] reconceived the task distribution itself — realistic per-variable ranges reviewed against the physics — and is the source of our task configurations. AI Feynman [13] and data-driven dimensional analysis [14] use physical units to constrain the hypothesis space; PySR [15] represents the modern search side. Dimensional correctness is necessary but not sufficient for form correctness: every ground-truth law in our suites is dimensionally self-consistent, yet each non-monomial law is a form misspecification for a power-law model. The falsification side complements these constraints — units narrow what to search; specification tests check what was found.")
p("Equation discovery proposes; specification testing challenges. Existing SR benchmarks — SRBench, SRBench++, SRSD, AI Feynman, PySR — report no analogue of a false-positive rate because they contain no explicit rejection mechanism to calibrate. HX addresses the complementary stage: after a candidate equation has been proposed and transformed to a restricted log-linear representation, does an explicit battery expose remaining structure?")

# ---------------- 3 Problem and Protocol ----------------
h1("3 Problem Formulation and Protocol")
p("Let yᵢ > 0 be a response and xᵢⱼ > 0 be p predictors for observations i = 1,…,n. Define ℓyᵢ = log yᵢ and Xℓ,ᵢⱼ = log xᵢⱼ. The restricted model contains an intercept and standardized log predictors Z; with B = [1, Z] the null hypothesis is that the log-log mean function is exactly linear, E[ℓy | X] ∈ span(B). The alternative is deliberately composite: at least one admissible nonlinear direction explains residual variation beyond span(B). All tests operate in log-log space; they differ only in what they add to this design and how they test the addition. Rejection means that the specified log-linear form is incompatible with at least one tested departure at the chosen level; non-rejection means only that the battery found insufficient evidence.")
p("Five protocol commitments are shared by every evaluation in this paper. P1 — common rejection basis: all tests decide at the same significance level α = 0.05 on p-values; the permutation tests' empirical p-values are computed from B = 999 resamples. P2 — power only at matched empirical FPR: headline power comparisons calibrate each method's decision threshold to a 5% empirical false-positive rate estimated on pooled true-monomial tasks within the same suite, sample-size stratum, and noise regime; raw-α results are reported alongside, because raw tables are what current practice would see. P3 — preregistered labels: H0 versus H1 is fixed from symbolic structure of the ground truth before any test runs, and no label is revised after seeing results. P4 — explicit failure registration: every evaluation that cannot produce a valid decision is registered with a typed reason and excluded from rate numerators while remaining in the denominators' audit trail. P5 — official data configurations: tasks reuse the official per-variable sampling objects of the source registry; tasks whose inputs or outputs are not strictly positive under the official configuration are marked not-applicable and excluded, and sensitivity to the official ranges is measured explicitly rather than assumed.")

# ---------------- 4 Part I: Benchmark ----------------
h1("4 Part I — The Falsification-Side Benchmark")
h2("4.1 Methods under test")
p("Four tests are compared under the protocol of Section 3. Method A is a permutation curvature test (form check): it adds a single mean-centered quadratic column q = meanⱼ((lxⱼ − lx̄ⱼ)²) — for one feature, the squared centered log-feature itself; for several, their row-wise mean — and takes the added coefficient as statistic, deriving an empirical p-value from B = 999 residual permutations with a statistic-mismatch sentinel comparing wrapper and native values on every call. Method B is RESET fitted, classical F [1]: it adds ŷ² and ŷ³, where ŷ is the fitted value from the restricted model, and tests their joint significance with the classical F statistic, powers (2, 3) throughout. Method C is RESET fitted, HC3-Wald: the same augmented design with W = β_u′[V_HC3]⁻¹β_u on the added coefficients, where V_HC3 = (X′X)⁻¹[Σᵢ eᵢ²/(1−hᵢ)² xᵢxᵢ′](X′X)⁻¹, referred to F(q, n − k). Method D is RESET exog, HC3-Wald: built regressors are the powers of each log-regressor (lxⱼ², lxⱼ³ for every feature j) rather than of the fitted values; HC3-Wald as in C. Configurations were frozen before the main runs — powers (2, 3); B = 999; α = 0.05; the grids below — and are identical across all three suites. The same harness invokes all four tests on every dataset, so method comparisons share data, seeds, and failure handling exactly.")
h2("4.2 Task suites S1–S3")
p("S1 — synthetic control (19 scenarios × 500 seeds). Three monomial H0 scenarios (univariate; two-variable; two-variable with ρ ≈ 0.9 correlated design), each with strictly positive log-uniform features on [0.1, 10] × [0.05, 2]. Eleven structural H1 scenarios: quadratic curvature; cubic; multiplicative interaction; threshold; exponential law; logarithmic law; rational (saturating) law; omitted variable; dilution; opposing curvature; collinear-quadratic. Five noise-layer H1 scenarios have strict-monomial mean structure but non-Gaussian noise — Student-t(3); variance growing with the log prediction; 2% point outliers; 6% group outliers; AR(1) correlation — rejections there are false alarms by construction and are never pooled into structural power. Grid: n ∈ {50, 200, 1000}; relative Gaussian noise σ ∈ {0, 0.5%, 1%, 3%, 5%} injected additively in log space; 500 seeds per cell via a seed rule that embeds the grid position.")
p("S2 — SRSD-Feynman tasks (32 tasks × 500 seeds). From the 120-equation registry, 32 tasks pass the positivity screen under official configurations: 13 strict-monomial H0 (including physical-constant domains at 10⁻³⁴–10⁻²³ scale) and 19 non-monomial H1. Each training fit (n ∈ {50, 200, 1000}) is accompanied by a same-distribution test set of 2,000 samples and a far-field probe: the first regressor extrapolated ×3 beyond both training domain edges, other regressors at training medians, scored against the closed-form ground truth.")
p("S3 — family-stratified Feynman equations (30 equations × 30 seeds). Structural family classification from symbolic form: exact power law (12 H0); fractional power (1 H0, exponent 1/3); exponential (3 H1); rational (6 H1); trigonometric (5 H1); plus multivariable-additive, logarithmic, and other families represented only by not-applicable members. Each H1 equation runs over relative noise {0.5%, 1%, 3%, 5%} × {homoskedastic, heteroskedastic} × sampling-range scale {official, ×0.5, ×2}; H0 equations run the noise grid at official range. Total across suites: 1,157 cells, 394,260 fitted datasets, plus 472,500 held-out evaluations (240,000 same-distribution test sets and 232,500 extrapolation probes; the remainder exclude one equation whose ground truth is not finite on the probe domain, registered as not-a-number rather than silently scored). With the battery suite S4 (Section 4.7) the grand total is 1,345 cells and 488,260 fitted datasets (Appendix A).")
h2("4.3 Size: false-positive rates on true monomials")
p("On synthetic Gaussian data, all four tests are calibrated (Table 1): every rate falls within or near its nominal band, and no method is distinguishable from another. On real physical laws the picture splits (Table 2). Methods A and B remain calibrated — rates .049–.053, intervals overlapping 0.05, flat in n. Both HC3 variants break, and break identically at every sample size: C stabilizes at .417–.422, D at .563–.565, with 95% intervals that do not move from n = 50 to n = 1000. This is not a small-sample artifact that more data cures; it is a structural miscalibration.")
tbl(["Method", "n=50", "n=200", "n=1000"],
    [["A — permutation curvature", ".047 [.041–.053]", ".037 [.032–.043]", ".041 [.036–.047]"],
     ["B — RESET fitted F", ".045 [.039–.051]", ".041 [.036–.047]", ".038 [.033–.044]"],
     ["C — RESET fitted HC3", ".048 [.042–.054]", ".045 [.039–.051]", ".038 [.033–.044]"],
     ["D — RESET exog HC3", ".054 [.048–.061]", ".042 [.036–.048]", ".040 [.034–.046]"]],
    widths=[5.5, 3.4, 3.4, 3.4], aligns=["l", "c", "c", "c"])
cap("Table 1. S1 synthetic H0, false-positive rate at α = 0.05 (Gaussian noise pooled over σ ∈ {0.5%…5%}; 4,500 fits per cell; Wilson 95% intervals).")
tbl(["Method", "n=50", "n=200", "n=1000"],
    [["A — permutation curvature", ".049 [.046–.052]", ".053 [.050–.055]", ".052 [.049–.054]"],
     ["B — RESET fitted F", ".050 [.047–.053]", ".050 [.047–.052]", ".053 [.050–.055]"],
     ["C — RESET fitted HC3", ".422 [.416–.428]", ".417 [.411–.423]", ".418 [.412–.424]"],
     ["D — RESET exog HC3", ".564 [.558–.570]", ".563 [.557–.569]", ".565 [.559–.571]"]],
    widths=[5.5, 3.4, 3.4, 3.4], aligns=["l", "c", "c", "c"])
cap("Table 2. S2 true physical monomials (13 tasks), false-positive rate at α = 0.05 (26,000 fits per cell; Wilson 95% intervals). This is the benchmark's central negative finding: HC3-Wald variants carry 42–56% false positives, flat in sample size.")
p("The mechanism is collinearity of the built columns. Physical constants and engineered quantities span few log-decades; on such designs ŷ² and ŷ³ are nearly linear functions of the fitted values themselves (for C), and the built power columns for D are nearly constant or near-collinear with the design — S3 records the limiting case, a ⅓-power monomial whose built columns are effectively constant and where D's false-positive rate reaches 1.00. The sandwich estimator's asymptotic guarantee presumes a well-behaved design as n grows; here the design degeneracy is exact and does not wash out, so the Wald statistic's reference distribution is simply wrong. H0 exact-power-law cells of S3 reproduce the same signature at independent scale (C .52, D .36). Noise-layer scenarios (Section 4.6) confirm what does not cause this: heavy tails, growing variance, outliers, and autocorrelation leave all methods' sizes at .03–.07. The design, not the noise, is what specification tests are fragile to.")
h2("4.4 Power: raw, then calibrated")
p("Raw power at α = 0.05 on S2 H1 (n = 200, pooled over rel ≥ 0.5%) is A .52, B .28, C .63, D .94 — a naive reading crowns D. Table 2 already priced that crown: D pays a .563 false-positive rate for it. Calibrating every method to a common 5% empirical FPR on the H0 pool reorders the field (Table 3). At deployment-relevant sample sizes (n ≤ 200), A delivers roughly twice B's calibrated power and 2–10× C's and D's. At n = 1000, D overtakes (.65): with enough data, exog powers resolve the very structures that the mean-square average compresses, and its calibration failure becomes less damaging relative to its information gain. Even so, no method exceeds .65 — most of S2's H1 mass is hard for every test.")
tbl(["Method", "n=50", "n=200", "n=1000"],
    [["A — permutation curvature", ".41", ".51", ".59"],
     ["B — RESET fitted F", ".19", ".28", ".42"],
     ["C — RESET fitted HC3", ".05", ".28", ".41"],
     ["D — RESET exog HC3", ".19", ".37", ".65"]],
    widths=[5.5, 3.4, 3.4, 3.4], aligns=["l", "c", "c", "c"])
cap("Table 3. S2 power at matched empirical FPR 5%.")
p("On S1 structural H1 (raw α, rel = 1%, n = 200), the same pattern appears in sharper form. Quadratic, threshold, exponential-shape, logarithmic, rational-single, omitted-variable, and collinear-quadratic scenarios are caught by all methods at ≥ .96. The informative cells are the blind spots (Table 4). Three mechanisms, each verified by construction: (i) A's statistic is driven by the mean of squared centered columns, so when two variables carry equal and opposite curvature on an exchangeable design that mean has exactly zero population signal (OPPOSING: A .02 at n = 200, .00 at n = 1000 — the statistic doing precisely what it is defined to do); (ii) a quadratic-projection statistic loses power on odd-symmetric deviations (A's cubic power never approaches B's saturation, because the cubic's second-order projection is small relative to its third-order content); (iii) B's fitted-value powers are a local Taylor surrogate for any smooth alternative — which is why B saturates on interaction (1.00) where A (.21) and D (.10) struggle — and, as S3 will show, the same surrogate property near-blinds B on rational laws.")
tbl(["Scenario", "A", "B", "C", "D"],
    [["cubic (standardized odd term)", ".40", "1.00", "1.00", "1.00"],
     ["interaction (lx₁·lx₂ cross term)", ".21", "1.00", "1.00", ".10"],
     ["opposing curvature (+c, −c, exchangeable)", ".02", ".98", ".98", "1.00"],
     ["dilution (one curved variable of two)", "1.00", ".88", ".87", "1.00"]],
    widths=[6.4, 2.0, 2.0, 2.0, 2.0], aligns=["l", "c", "c", "c", "c"])
cap("Table 4. S1 blind-spot cartography (raw α, n = 200, rel = 1%).")
h2("4.5 Family-stratified results and sensitivity")
tbl(["Family (class)", "A", "B", "C", "D"],
    [["exact power law (H0, 13 eq. → size)", ".04", ".04", ".52", ".36"],
     ["fractional power ⅓ (H0 → size)", ".06", ".09", ".10", "1.00"],
     ["exponential (H1, 3 eq.)", ".11", ".05", ".70", ".68"],
     ["rational (H1, 6 eq.)", ".66", ".13", ".68", ".95"],
     ["trigonometric (H1, 5 eq.)", ".69", ".72", ".81", "1.00"]],
    widths=[6.4, 2.0, 2.0, 2.0, 2.0], aligns=["l", "c", "c", "c", "c"])
cap("Table 5. S3 family map, official sampling domains, homoskedastic noise pooled over rel (raw α).")
tbl(["Family", "A", "B", "C", "D"],
    [["exponential", ".11", ".05", ".63", ".04"],
     ["rational", ".66", ".14", ".50", ".66"],
     ["trigonometric", ".70", ".73", ".08", ".05"]],
    widths=[6.4, 2.0, 2.0, 2.0, 2.0], aligns=["l", "c", "c", "c", "c"])
cap("Table 6. S3 calibrated power by family (5% empirical FPR).")
p("Read jointly, Tables 5 and 6 are the benchmark's central finding. B — the classical test, near-perfect on synthetic quadratics — is nearly blind to rational laws (.13 raw, .14 calibrated): the fitted-value Taylor surrogate absorbs precisely the low-order curvature by which a saturating law departs from a power law. D's raw dominance on rational (.95) and trigonometric (1.00) is false-alarm budget, not power: calibrated, its trigonometric advantage collapses to .05. C is the exponential family's only reliable detector (.63 calibrated) while being worthless on trigonometric (.08). A misses exponentials (.11) — an exponential law's log-log curvature is modest at official scales — but is otherwise the most even performer. Every method's strength sits on another method's blind spot; the four columns of Table 6 have no row-major winner.")
p("Sensitivity: scaling every sampling domain ×2 moves exactly one family-method cell by more than 15 percentage points (rational × A: .66 → .84); ×0.5 moves none. The heteroskedastic replication of every S3 cell agrees with the homoskedastic one within .02 on every family-method cell. Sampling-range choice and variance model — the two design knobs a practitioner might suspect of driving the conclusions — do not drive them.")
h2("4.6 Noise robustness and extrapolation alarms")
tbl(["Scenario", "A", "B", "C", "D"],
    [["Student-t(3) noise", ".034", ".036", ".030", ".030"],
     ["variance growing with prediction", ".048", ".058", ".050", ".050"],
     ["AR(1) errors, ρ = 0.6", ".038", ".040", ".042", ".042"],
     ["2% point outliers (+8σ)", ".040", ".066", ".036", ".036"],
     ["6% group outliers (+5σ)", ".038", ".060", ".054", ".054"]],
    widths=[6.4, 2.0, 2.0, 2.0, 2.0], aligns=["l", "c", "c", "c", "c"])
cap("Table 7. S1 noise-layer H1 rejection rates at α = 0.05 — false alarms attributable to the noise assumption (n = 200, rel = 1%, 500 seeds).")
p("Every method, every non-Gaussian regime: .030–.066. This table's function is exclusionary: it rules out the convenient explanation that the miscalibrations of Table 2 come from noisy, ugly, real-world data. The data in Table 2 are clean by construction; what differs from Table 1 is the design.")
p("Do alarms predict extrapolation failure? Each S2 fit (H0 and H1 tasks alike — the deployment situation, since a pipeline does not know which candidates are right-formed) is scored on 2,000 same-distribution test samples and on a ×3-domain far-field probe. Base rate of far-field relative error above 10%: .198. Specification-test alarms enrich this substantially: an A alarm raises the probability of severe extrapolation failure from .198 to .39 (2.0×); a B alarm to .56 (2.8×). Of the true-positive catches, 10,560 are A-only and 5,447 B-only, against 12,575 shared — the two tests' alarm sets overlap barely more than chance, so a portfolio of both captures substantially more failures than either.")
tbl(["Method", "precision", "recall"],
    [["A — permutation curvature", ".39", ".63"],
     ["B — RESET fitted F", ".56", ".49"]],
    widths=[7.5, 3.0, 3.0], aligns=["l", "c", "c"])
cap("Table 8. Alarm precision and recall for far-field error > 10% (S2, all tasks, all n, rel ≥ 0.5%).")
p("Alarms are useful as ordinal evidence about form; they are not certificates of magnitude. At the task level, the rank correlation between median in-domain p-value and median far-field error is moderate for A (Spearman ρ = −0.33) and stronger for B (ρ = −0.68), yet both orderings fail exactly where it hurts most: the suite's two worst extrapolators (both at 98% median far-field error) sit at the weak end of every method's evidence because their in-domain curvature is near-zero while the failure lives beyond the domain edge. Conversely, the benchmark's most-rejected benign extrapolator is rejected by A essentially always yet extrapolates benignly (median far-field error 7.7%): an alarm means form evidence exists; a non-alarm means nothing about the field beyond the training domain, and an alarm does not price the error either. No rejection is evidence of extrapolation safety.")
h2("4.7 Suite S4 — the family-viability battery")
p("The four tests share one property: each evaluates a single null form. A rejection is actionable only as an escalation, but the family-stratified results hint the machinery could do more: each test is a curvature check in some linearizing coordinates, and several physical families linearize in known coordinates. The battery registers four families with their linearizing coordinate pairs (u(x), v(y)): power (log x, log y); exponential (x, log y); logarithmic (log x, y); rational (1/x, 1/y). For each family, the member test is \u201cthe conditional mean is linear in the family's coordinates\u201d: the design is [1, U, q] where U are the standardized transformed features and q = meanⱼ(Uⱼ²) is the single quadratic column. Columns are standardized before squaring so that features with different log-variances (physical constants at 10⁻³⁴ next to meter-scale quantities) contribute equally. The battery's output is the viability set V = {F : p_F ≥ α}, not a winner; three errors are defined and reported separately: false exclusion, false inclusion, and — for structures outside all four families — escalation, which is the correct outcome.")
p("Null validity is a coordinate-dependent property. A residual-permutation null is exact only under exchangeable errors in the transformed coordinates. Log transforms homoskedasticize multiplicative noise — which is why the log-log design was chosen — but the same noise is heteroskedastic in (x, log y), (log x, y), and (1/x, 1/y) coordinates. This prediction was preregistered before any run: permutation should be valid for families ending in v = log y (power, exponential) and break for families ending in v = y or v = 1/y (logarithmic, rational). The probe measured exactly that (Table 9): under multiplicative noise, the permutation-null rational member rejects its own true family at 30–34% — false escalations caused purely by null misspecification — while power and exponential members stay calibrated. Replacing the permutation draws with a wild bootstrap [7, 9] (Rademacher signs on restricted residuals; B = 999, vectorized through the pseudo-inverse row) drops the rational false escalation to 4.5–9% and leaves the log-coordinate members' size unchanged. One line of null construction, no change to statistic or design.")
tbl(["Truth", "perm null", "wild null"],
    [["power (size of power member)", ".070 / .075", ".075 / .060"],
     ["rational: false escalation (V = ∅)", ".34 / .31", ".045 / .090"],
     ["rational: containment (rational ∈ V)", ".66 / .70", ".96 / .91"],
     ["exponential: containment", ".94 / .95", ".92 / .95"],
     ["logarithmic: containment", ".87 / .93", ".89 / .95"],
     ["threshold (no battery family): escalation — correct outcome", "1.00 / 1.00", "1.00 / 1.00"]],
    widths=[6.4, 3.2, 3.2], aligns=["l", "c", "c"])
cap("Table 9. Probe: permutation vs wild null by family (n = 200, 200 seeds; entries rel = 1% / 3%).")
p("The rational row is the mechanism in one line: with the permutation null, the battery's own rational member rejects the best case for the rational test — data generated as y = x/(0.2 + 0.1x) with multiplicative noise — a third of the time. Nothing about rational laws is at fault; the null is exchangeability-invalid in (1/x, 1/y) coordinates. The wild bootstrap repairs exactly that and nothing else, confirming the coordinate-dependence reading rather than a blanket \u201cwild beats permutation.\u201d")
tbl(["Truth family", "Containment (truth ∈ V)", "False escalation (V = ∅)"],
    [["power — size, 8 scenario groups incl. all noise layers", "—", ".0496 [.046–.053] pooled; per-scenario .043–.054"],
     ["exponential", ".95", ".05"],
     ["logarithmic", ".94", ".005"],
     ["rational", ".93", ".066 pooled (≤ .082 worst rel)"]],
    widths=[6.4, 3.2, 4.4], aligns=["l", "c", "c"])
cap("Table 10. Battery on synthetic truths (500 seeds per cell, n = 200, wild null throughout; pooled over rel ∈ {0.5–5%}; 2,000 datasets per truth family per cell).")
p("Under logarithmic truth the modal viable set is {logarithmic, rational} — the rational member co-survives in 91% of logarithmic-truth datasets. This is false inclusion in the strict sense, and it is the multiplicity property the battery promised to report rather than suppress: over a ×100 domain, a logarithmic and a saturating rational law are near-neighbours the battery cannot and should not arbitrate between. Arbitration belongs to the discovery layer; the battery's contract is \u201cnever silently route a truth out of V,\u201d and after the wild fix that contract holds: containment .93–.95 against a 5% per-family budget, false escalation at or below 8.2% in the worst cell.")
tbl(["Scenario", "Correct escalation", "False-inclusion rate", "Modal false V"],
    [["quadratic", "1.00", ".00", "—"],
     ["threshold", "1.00", ".00", "—"],
     ["omitted variable", "1.00", ".00", "—"],
     ["dilution (1 of 2 features curved)", ".92", ".08", "rational"],
     ["cubic (odd term)", ".24", ".76", "power (.68)"],
     ["interaction (cross term)", ".05", ".95", "power (.95)"],
     ["opposing curvature", ".00", "1.00", "power+rational (.84)"],
     ["collinear-quadratic", ".03", ".97", "rational (.61)"]],
    widths=[6.0, 2.6, 2.6, 3.2], aligns=["l", "c", "c", "l"])
cap("Table 11. Battery on out-of-family structures (truth ∈ no battery family; escalation is the correct outcome; 2,000 datasets per scenario).")
p("Table 11 is the blind-spot cartography of Table 4 re-expressed as routing errors, and the correspondence is exact: the power member inherits A's documented blind spots (opposing curvature, odd terms, interactions), so those structures stay in V under the power label instead of escalating. Two rows are new information: collinear-quadratic — where A alone had full power — fails to escalate because the curvature hides behind the degenerate design for every member simultaneously; and dilution, A's synthetic win, escalates correctly 92% of the time. On the SRSD side, the battery's routing on the 14 preregistered H1 equations is descriptive: surface-form family labels are not linearizable form — Planck-type x³/(e^{c/x} − 1) \u201cexponentials\u201d are correctly rejected by the exponential member — so no containment score is computed; median |V| is 1–2 for 11 of 14 equations, and the benchmark's most-rejected benign extrapolator escalates on 67%, the battery correctly refusing to house it anywhere. A symbolically exact linearizability registry is future work.")
h2("4.8 What the benchmark establishes")
p("Four findings stand. First, aggregate rankings do not transfer across task families: the classical-F RESET variant that dominates synthetic quadratic benchmarks is nearly blind to rational ground truths, because its fitted-value powers locally absorb exactly the curvature by which rational laws deviate from power laws. Second, HC3-Wald variants, which appear strongest in raw tables, carry 42–56% false-positive rates on true monomial physical laws, flat in n — so their apparent advantage disappears once methods are compared at matched empirical false-positive rates. Third, no test is uniformly best; each method's failure profile is complementary to the others', and a portfolio of two tests at under 30 ms per dataset dominates any single one. Fourth, the battery shows that the same statistic can route between families — but only after its null distribution is rebuilt per family, because null validity is a coordinate-dependent property. These findings set the design agenda for Part II: a falsification screen should combine complementary direction families, protect the search multiplicity, calibrate its null honestly, and expose its failures.")

# ---------------- 5 Part II: HX ----------------
h1("5 Part II — HX: A Hybrid Directional Specification Test")
h2("5.1 Design motivation")
p("The benchmark of Part I maps a complementary failure landscape. A single-direction form check is calibrated but blind to cubic, interaction, and opposing departures; the structural library that repairs those blind spots is itself blind to dilution and opposing alternatives in its fitted-value form; and no test is calibrated in every regime. HX turns that complementarity into one auditable procedure: it augments the restricted log-linear regression with two complementary direction families — structural directions (quadratic, cubic, interaction, oscillatory terms in standardized log predictors), which target interpretable departures, and fitted-value directions (powers of the restricted fitted response), which provide an omnibus, response-sensitive probe related to Ramsey's RESET. The two branches answer different questions: a concentrated defect may align strongly with one interpretable direction, whereas a diffuse defect may be weak in every coordinate but large jointly. HX therefore retains two statistics rather than forcing one representation to serve both roles.")
h2("5.2 Core algorithm")
h3("5.2.1 Restricted fit and preprocessing")
p("Each log predictor is centered and scaled using the analysis sample. Constant or numerically degenerate columns are rejected before fitting. The restricted ordinary least-squares fit yields β̂ = B⁺ℓy, ŷ = Bβ̂, and r = ℓy − ŷ, where B⁺ denotes the Moore–Penrose inverse. The residual degrees of freedom must exceed a declared floor. A degenerate-fit guard is expressed as the translation-invariant ratio RSSr/TSSℓy, not as RSSr divided by the raw squared norm of ℓy.")
h3("5.2.2 Candidate direction library")
p("For each standardized predictor zⱼ, the structural library may contain centered quadratic zⱼ², centered cubic zⱼ³, and an oscillatory probe sin(2πzⱼ); pairwise interactions zⱼzₖ are included for j < k. The fitted-value library contains centered powers of a standardized restricted fitted vector, typically ŷs² and ŷs³, where ŷs = (ŷ − mean ŷ)/sd(ŷ). Standardizing before taking powers makes the library insensitive to an additive shift or nonzero scalar rescaling of ℓy, apart from a harmless sign change in odd powers. The names attached to retained columns are part of the diagnostic record: they support mechanistic follow-up but do not convert the selected coordinate into a separately confirmed discovery.")
h3("5.2.3 Base residualization and rank-aware selection")
p("Every candidate direction is residualized against the restricted model: D⊥ = (I − PB)D, with PB = BB⁺. A modified Gram–Schmidt process selects a numerically independent set: at each step, the candidate with the largest residual norm after projection on the already selected W is retained if the norm exceeds a tolerance. If a hard direction budget is active, selection uses a relative residual norm, ‖dⱼ,res‖/max(‖dⱼ,⊥‖,τ), to avoid a scale-dependent preference. In the reported implementation no hard budget bound the evaluated low-dimensional cases; the relative criterion is the normative rule for larger libraries.")
h3("5.2.4 Directional maxT branch")
p("For a retained orthonormal direction wⱼ, the score wⱼᵀr measures alignment between the restricted residual and that departure; it is studentized using the restricted residual variance estimate. HX records the maximum absolute score across q directions. For each resample, the entire maximum is recomputed; the resulting pmax is therefore familywise adjusted over the retained direction family. The unadjusted coordinate-level p-value, when returned, is labeled exploratory and must not be interpreted as confirmatory after selection.")
h3("5.2.5 Joint augmentation branch")
p("The unrestricted design U = [B, W] yields residual sum of squares RSSu and residual degrees of freedom dfu. The partial-F statistic Fjoint captures distributed misspecification whose energy is spread across several individually moderate coordinates — the case where maxT, by construction, is weak.")
h3("5.2.6 Residual resampling")
p("The development implementation uses a frozen-direction reduced-model residual permutation: it forms ℓy*(b) = ŷ + P(b)r, where P(b) permutes residual entries, and evaluates both statistics with B and the selected W held fixed. With B-fixed structural directions this is a standard reduced-model residual strategy whose accuracy depends on the exchangeability of residuals under the null. Because ŷ² and ŷ³ were constructed from the observed response and frozen across permutations, the fitted-value branch is response-adaptive; its finite-sample exactness is not guaranteed, and the empirical p-values should be read as calibrated approximations. A confirmatory HX implementation should use one of two declared modes: refit mode reconstructs the restricted fit, fitted-value directions, pivot selection, and both statistics inside every resample; split mode constructs response-adaptive directions on training folds and tests them on disjoint observations. Under heteroskedasticity, unrestricted residual permutation should be replaced by a wild bootstrap — the repair the battery of Section 4.7 quantifies — and clustered or repeated observations require exchangeability blocks or cluster-level sign schemes.")
h3("5.2.7 Combined decision and output contract")
p("Let pmax be the adjusted maxT p-value and pF the resampling p-value for Fjoint. HX uses a Bonferroni minimum across these two branches: pHX = min(1, 2·min(pmax, pF)). This combination does not require independence between branches. A result object contains pHX, pmax, pF, q, selected direction names, the maxT winner, observed statistics, sample size, ranks, degrees of freedom, conditioning diagnostics, resample count, seed, and a failure code. Returning only pHX is insufficient for scientific audit. The complete algorithm, failure taxonomy, and minimum reporting checklist appear in Appendices A–C.")
h2("5.3 Statistical interpretation and safeguards")
p("The maxT branch controls search across the retained coordinate family when its resampling reference is valid; the final Bonferroni step controls the union of the maxT and joint-F branch decisions when both component p-values are valid. Neither statement repairs a misspecified resampling law: in the frozen fitted-value implementation, nominal size is evaluated empirically and the p-value is best described as approximate. Four invariance requirements are enforced: translation invariance (adding a constant to ℓy changes only the intercept and must not trigger a different degeneracy classification); response-scale invariance (multiplying ℓy by a nonzero constant leaves standardized direction geometry and p-values unchanged up to numerical tolerance); predictor-scale invariance (affine recoding of log predictors is absorbed by standardization); and column-order stability (a declared deterministic tie rule makes retained direction names reproducible when residual norms are numerically equal). Guard conditions return explicit named failures — bad_shape_X, length_mismatch, nonfinite_input, domain_error, rank_deficient_base, ill_conditioned_base, constant_response, degenerate_zero_residual, no_valid_direction, insufficient_residual_df, invalid_exchangeability, resampling_failure — never a silent NaN or an unclassified exception (Appendix B).")
h2("5.4 Evaluation design")
p("Four methods were compared at α = 0.05: (i) a single-direction permutation form check (method A of Part I); (ii) MD4, the structural-only direction library using maxT; (iii) HX, which adds fitted-value directions and the joint-F branch; and (iv) fitted-value RESET using its classical parametric F reference. The earlier companion benchmark also considered heteroskedasticity-consistent RESET variants; the focused HX experiments used the classical fitted RESET comparator. All methods shared datasets, seeds, and decision thresholds; they did not share one common null generator — the permutation procedures used their declared empirical references, whereas RESET used its F reference. Power comparisons should therefore be interpreted jointly with each method's empirical null calibration, exactly as the protocol P2 of Section 3 requires.")
p("The controlled development benchmark comprised three null scenarios and six structural alternatives — quadratic, cubic, interaction, dilution, opposing, and threshold departures — at n = 50 and n = 200, with selected cubic and opposing cases also evaluated at n = 1,000; each principal cell used 50 seeded replications with 1% relative noise, and additional NIST-style analytic functions served as interpretable challenge cases. This suite was used during method development and is not an untouched final test set. The post-development randomized stress test used a master seed of 20261904 to generate 1,000 groups: 393 null cases and 607 alternatives, with sample size, predictor count, noise, signal-to-noise ratio, exponents, and departure family varying across groups. The 12 alternative families were cross-cubic, cubic, exponential deformation, interaction, logistic bump, mixed cubic–oscillatory, opposing quadratic, quadratic, quartic, one-cycle sinusoid, two-cycle sinusoid, and threshold. The generator expands coverage, but its family design remained informed by development experience; it is therefore a stress test rather than independent external validation.")
p("For null cases the primary metric was false-positive rate with a 95% Wilson interval; for alternatives, rejection rate (power) with the same interval; overall accuracy used the generated class label only as a summary, not as a substitute for separate size and power. Paired differences on the same cases were evaluated with exact two-sided McNemar tests, equivalent to a binomial test on discordant pairs.")
h2("5.5 Results")
h3("5.5.1 Controlled benchmark")
p("At n = 200, the pooled null rejection rates were 0.033 for the single-direction form check, 0.027 for MD4, 0.033 for HX, and 0.053 for RESET. In the six controlled alternative families, HX rejected all 50 replications at both n = 50 and n = 200 (Table 12). The main gains were complementary: RESET was comparatively weak for dilution and opposing alternatives, while the single-direction check was weak for cubic, interaction, and opposing alternatives — exactly the blind-spot map of Part I.")
tbl(["Family", "n", "form_check", "MD4", "HX", "RESET"],
    [["Quadratic", "50 / 200", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00"],
     ["Cubic", "50 / 200", "0.50 / 0.36", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00"],
     ["Interaction", "50 / 200", "0.24 / 0.18", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00"],
     ["Dilution", "50 / 200", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00", "0.44 / 0.88"],
     ["Opposing", "50 / 200", "0.00 / 0.02", "1.00 / 1.00", "1.00 / 1.00", "0.42 / 0.96"],
     ["Threshold", "50 / 200", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00", "1.00 / 1.00"]],
    widths=[3.4, 2.4, 2.8, 2.8, 2.8, 2.8], aligns=["l", "c", "c", "c", "c", "c"])
cap("Table 12. Controlled alternative rejection rates by sample size (50 seeded replications, 1% relative noise).")
h3("5.5.2 Randomized calibration and power")
p("In 1,000 randomized groups (Table 13), HX achieved the highest power (0.9522) and overall accuracy (0.951) while its estimated false-positive rate (0.0509) remained close to the nominal level, and its 95% Wilson interval contained 0.05. The single-direction form check, with a comparable false-positive rate, attained power 0.4942 — the cost of committing to one direction. MD4, the structural-only library, attained 0.8913; RESET 0.7512 at a false-positive estimate of 0.0611 whose interval still included 0.05. The relevant conclusion is therefore not that HX had a demonstrably lower size than every comparator, but that it delivered substantially greater power without evidence of size inflation in this generator. Read through the lens of Part I, every one of these power numbers is meaningful only because the matched-FPR discipline prices each method's false-alarm budget.")
tbl(["Method", "False positives", "FPR [95% CI]", "True positives", "Power [95% CI]", "Accuracy"],
    [["form_check", "21/393", "0.0534 [0.0352, 0.0803]", "300/607", "0.4942 [0.4546, 0.5339]", "0.672"],
     ["MD4", "16/393", "0.0407 [0.0252, 0.0651]", "541/607", "0.8913 [0.8640, 0.9136]", "0.918"],
     ["HX", "20/393", "0.0509 [0.0332, 0.0773]", "578/607", "0.9522 [0.9322, 0.9665]", "0.951"],
     ["RESET", "24/393", "0.0611 [0.0414, 0.0893]", "456/607", "0.7512 [0.7153, 0.7840]", "0.825"]],
    widths=[3.0, 2.4, 3.4, 2.4, 3.4, 1.8], aligns=["l", "c", "c", "c", "c", "c"])
cap("Table 13. Randomized stress test results with 95% Wilson intervals (1,000 groups: 393 null, 607 alternative).")
h3("5.5.3 Alternative-family coverage")
p("HX reached 1.00 rejection rate in ten of twelve randomized families (Table 14). HX and MD4 both detected 0.97 of the two-cycle sinusoid cases, while threshold departures remained difficult for every method — the threshold family is the one blind spot the benchmark also found irreducible. Cross-cubic alternatives illustrate the fitted-value contribution: MD4 detected 0.29, whereas HX and RESET both detected 1.00. Opposing quadratic alternatives illustrate the structural contribution: HX and MD4 detected all cases, but RESET detected only 0.05.")
tbl(["Family", "n", "form_check", "MD4", "HX", "RESET"],
    [["cross_cubic", "52", "0.10", "0.29", "1.00", "1.00"],
     ["cubic", "56", "0.41", "1.00", "1.00", "0.82"],
     ["exp_deform", "44", "1.00", "1.00", "1.00", "0.91"],
     ["interact", "47", "0.21", "1.00", "1.00", "1.00"],
     ["logistic_bump", "49", "0.10", "1.00", "1.00", "0.78"],
     ["mixed_cubic_osc", "54", "0.52", "1.00", "1.00", "0.87"],
     ["opposing_quad", "58", "0.88", "1.00", "1.00", "0.05"],
     ["quad", "58", "1.00", "1.00", "1.00", "0.98"],
     ["quartic", "48", "1.00", "1.00", "1.00", "1.00"],
     ["sin1", "53", "0.00", "1.00", "1.00", "0.66"],
     ["sin2", "38", "0.18", "0.97", "0.97", "0.63"],
     ["threshold", "50", "0.42", "0.44", "0.44", "0.38"]],
    widths=[3.4, 1.6, 2.8, 2.8, 2.8, 2.8], aligns=["l", "c", "c", "c", "c", "c"])
cap("Table 14. Alternative rejection rate by randomized family.")
h3("5.5.4 Paired discordance")
p("On the 607 randomized alternatives (Table 15), HX's power advantage over every comparator is not a small-parts construction: HX uniquely rejected 278, 37, and 123 cases that form_check, MD4, and RESET respectively missed, against 0, 0, and 1 case rejected only by the comparator. Exact two-sided McNemar tests give p = 4.12 × 10⁻⁸⁴, 1.46 × 10⁻¹¹, and 1.18 × 10⁻³⁵. The MD4-versus-HX comparison is a structural ablation: adding fitted-value directions and the joint branch produced 37 HX-only versus no MD4-only rejections. On null cases, HX versus form_check had 11 versus 12 discordant rejections (p = 1.000), HX versus MD4 had 6 versus 2 (p = 0.289), and HX versus RESET had 9 versus 13 (p = 0.524): the paired evidence supports a power difference under this alternative generator but not a meaningful size difference among the methods.")
tbl(["Comparator", "HX only", "Comparator only", "Exact two-sided p"],
    [["form_check", "278", "0", "4.12 × 10⁻⁸⁴"],
     ["MD4", "37", "0", "1.46 × 10⁻¹¹"],
     ["RESET", "123", "1", "1.18 × 10⁻³⁵"]],
    widths=[4.4, 2.4, 2.8, 3.4], aligns=["l", "c", "c", "c"])
cap("Table 15. Exact paired comparisons on the 607 randomized alternatives (McNemar).")
h3("5.5.5 Direction attribution and analytic challenges")
p("Among randomized alternatives with a winning direction (Table 16), the fitted-value cubic direction was the most common winner (35.6%), followed by structural cubic (27.2%), quadratic (14.2%), and oscillatory (11.4%) directions. This distribution supports the hybrid design: no single family accounted for all detections. It is descriptive, not a post-selection confidence statement about the true mechanism.")
tbl(["Direction", "Count", "Share"],
    [["yhat3", "216", "35.6%"],
     ["cubic_x1", "165", "27.2%"],
     ["quad_x1", "86", "14.2%"],
     ["osc_x1", "69", "11.4%"],
     ["yhat2", "42", "6.9%"],
     ["quad_x2", "26", "4.3%"],
     ["osc_x2", "2", "0.3%"],
     ["cubic_x2", "1", "0.2%"]],
    widths=[4.4, 3.4, 3.4], aligns=["l", "c", "c"])
cap("Table 16. Most frequent winning HX directions among randomized alternatives.")
p("In the analytic challenge set, HX and RESET both alarmed on several nonlinear functions. A notable case was F3, where structural-only MD4 missed while fitted-value terms enabled HX to reject. Conversely, F10 produced a rank-deficiency failure for HX while RESET returned an alarm, demonstrating why failure states must remain visible. The L1 case had n = 16, only three unrestricted residual degrees of freedom, and a condition number near 4.9 × 10⁹; HX rejected while RESET did not, but the numerical regime makes the result inconclusive rather than a clean method victory.")
h3("5.5.6 Runtime")
p("At n = 200 and p = 2, median development-run time was approximately 4.7 ms for HX, 4.6 ms for MD4, and 0.5 ms for parametric RESET on the recorded machine — implementation- and hardware-dependent. Residual resampling dominates HX cost; confirmatory refit mode will be slower. Matrix-free projection, batched resampling, and streaming maxima avoid unnecessary n × n projection matrices and full permutation-score storage. At the scale of a discovery pipeline that already fits thousands of candidates, the entire falsification layer is effectively free; the bottleneck resource is not compute but decision discipline.")
h2("5.6 Mechanistic interpretation")
p("Structural directions supply locality and interpretability: a quadratic winner indicates curvature in one log predictor; an interaction winner indicates non-additivity; an oscillatory winner indicates a residual periodic pattern. These probes are weak when the omitted component lies outside the library or combines across predictors in a way that cancels a marginal projection. Fitted-value powers pool information along the restricted model's predicted response and can capture broad departures without naming their source; they are correspondingly less diagnostic and response-adaptive. The maxT branch is suited to sparse departures — one direction carries most of the omitted signal; the joint-F branch is suited to dense departures — several coordinates carry moderate signal. Combining them spends at most a factor of two in p-value but protects against committing to either geometry. This is a two-branch portfolio inside one test, and it is the same complementarity the benchmark of Part I documents at the level of whole tests: structural and omnibus diagnostics are complementary, and a falsification screen that commits to either alone inherits the other's blind spots.")
p("In a model-discovery workflow, HX is a challenge layer, not a certificate: generate candidate equations using physical constraints, symbolic regression, or expert hypotheses; reserve an untouched validation partition before tuning direction libraries or thresholds; fit each candidate on training data and apply HX under a declared confirmatory resampling mode on validation data; discard or revise candidates that HX falsifies, inspecting winning directions to formulate the next hypothesis; and for candidates not rejected, report power diagnostics and the exact battery used — without stating that the equation has been proven.")

# ---------------- 6 Discussion ----------------
h1("6 Discussion")
p("Read together, the two parts of this paper state one argument. The falsification side of physical model discovery is both necessary and constructible. Part I shows that the classical tests are not interchangeable: their rankings invert across task families and invert again under matched false-positive rates, and their failure modes are structural rather than noise-driven. Part II shows that the benchmark's blind-spot map is actionable: a hybrid test that combines interpretable structural directions with omnibus fitted-value powers, protects the search with maxT, evaluates joint signal with a partial-F statistic, and repairs its null by the resampling scheme appropriate to the noise geometry, attains near-nominal size and substantially higher power than any single-strategy comparator in a randomized stress test.")
p("Three design principles unify both parts and should govern any falsification layer for model discovery. First, directions must be visible: a test that cannot name what it looked at cannot support scientific follow-up. Second, multiplicity must be priced: search across a direction family is familywise adjusted by maxT; the two-branch portfolio is Bonferroni-combined; and cross-family multiplicity in the battery is reported as a finding rather than suppressed. Third, resampling assumptions must be declared and checked: the permutation null that calibrates in log-log coordinates is silently wrong in (1/x, 1/y) coordinates, and the repair — a wild bootstrap — is the same repair a confirmatory HX implementation prescribes under heteroskedasticity. A credible model check must expose its tested directions, its multiplicity correction, its resampling assumptions, and its failures; every result object in this paper does.")
p("The honest limits of the evidence are equally unifying. The battery's SRSD-side routing is descriptive pending a symbolically exact linearizability registry; HX's randomized suite is a stress test informed by the same development context, not independent external validation; and HX's frozen fitted-value directions make its residual-permutation calibration approximate. The empirical record is strong enough to support further evaluation, not unconditional deployment — and the reproducible record makes that further evaluation possible by anyone.")

# ---------------- 7 Limitations ----------------
h1("7 Limitations")
p("Single candidate family (Part I). The four benchmark tests evaluate a power law against alternatives; nothing there benchmarks model selection among families, only rejection of one. The battery addresses routing, not selection: it returns a viability set and does not rank its members. Fixed augmentation: built-regressor powers are fixed at (2, 3); adaptive or higher-order constructions may shrink B's rational blind spot at a calibration cost we have not measured. Battery scope: the battery runs at n = 200 only; it has no member for interaction or trigonometric structures; the collinear-quadratic false-inclusion row is registered as an open failure mode of member behavior on degenerate designs. Log-space scope: all tests operate in log-log space by design; additive-error models are out of scope and marked not-applicable rather than log-shifted.")
p("Evidence status (Part II). The controlled suite informed method development, and the randomized family generator reflects the same problem-solving context; the strong power estimates may therefore be optimistic for unrelated scientific domains. The reported HX p-values use frozen fitted-value directions under reduced-model residual permutation, so calibration is approximate because the directions are response-adaptive; empirical null calibration was satisfactory in the available generator, but confirmatory use should reconstruct adaptive steps within resamples, use sample splitting, or justify an alternative resampling law. Unrestricted residual permutation assumes exchangeability; heteroskedastic, clustered, temporal, or spatial errors violate that condition, and wild bootstrap, exchangeability blocks, or domain-specific dependence models are required. The fixed library has blind spots: the threshold family remained difficult, and no finite battery can cover all misspecification. Small residual degrees of freedom and ill-conditioned designs can make both directional selection and F statistics unstable; the winning direction is post-selection and exploratory, and causal or physical interpretation requires a separately specified follow-up analysis.")
p("Algorithmic novelty is narrower than scientific validity. HX combines established ideas — RESET-style fitted-value powers, structured basis augmentation, maxT adjustment, joint testing, and residual resampling — into an explicit falsification workflow; its value must be judged by transparent error control, robustness, and out-of-domain performance rather than by naming the combination. The benchmark's contribution is likewise a measured map of failures, not a new test. Statistical precision: headline S1/S2 rates rest on 500 seeds (Wilson intervals ±.003–.006 at p ≈ .05); S3 family cells rest on 30 seeds and support directional claims only; the battery probe comparisons rest on 200 seeds per cell, with the formal run (500 seeds) carrying the headline rates.")

# ---------------- 8 Reproducibility ----------------
h1("8 Reproducibility, Code, and Data")
p("Every generator in this paper is seeded from a documented rule embedding the grid position (suite, scenario, n, noise level, seed index; the battery's rule is given in Section 4.7), so any single cell of any table regenerates bit-identically. The harness records content hashes of the code under test — including the exact form-check implementation — and of the harness itself, in every run's environment file; a statistic-mismatch sentinel verifies on every call that the wrapper's recomputation agrees with the native implementation, so drift in the tested code cannot pass silently. The battery is implemented at the experiment layer over the unmodified core code path, so the audited core-code hashes remain valid. All cells' raw results (parquet and CSV), per-equation tables, failure registrations with typed reasons, selection and exclusion tables, analysis scripts, and figure sources are released. Zero-noise grid cells deliberately exercise the failure-registration path; the battery's 94,000 datasets contain zero data-layer failures, which is itself a registered outcome, not an omission.")
p("A reproducible HX report should archive the raw positive variables, log-transformation rules, excluded observations, standardized design, restricted and augmented ranks, condition numbers, direction names, selection tolerance, residual-degree floor, α, resampling scheme, exchangeability blocks or bootstrap weights, resample count, random seed, observed Tmax and Fjoint, component and combined p-values, and failure code; the benchmark datasets and seeds must be immutable once confirmatory evaluation begins. For a submission-quality replication: at least 9,999 resamples near the decision boundary, independent seed reruns, verification of the four invariance properties of Section 5.3, regression tests for malformed shapes, length mismatch, constant response, exact fit, rank deficiency, and ill conditioning, and comparison of matrix-free and reference implementations on small problems to numerical tolerance. Every attempted case, including failures, is reported rather than deleted.")
h2("8.1 Code inventory")
p("The release is organized as a Python package plus two experiment trees. The core package (dimreflex/, 4,081 lines) implements the dimensionally-aware discovery pipeline whose form-check primitive is the statistic under test in Part I. The benchmark experiment tree (experiments/form_check_vs_reset/) implements the four-test harness, the SRSD-Feynman task registry, the family-stratified selection, and the battery; the method tree (experiments/fc_multidirection/) implements HX, the structural-only MD4 comparator, and the randomized stress-test generator. Table 17 lists the principal modules; Appendix D gives the full inventory including data files and test evidence.")
tbl(["Module", "Lines", "Role"],
    [["dimreflex/core.py", "438", "dimensional analysis, unit registry, form-check statistic"],
     ["dimreflex/discover.py", "518", "monomial discovery and candidate surface construction"],
     ["dimreflex/families.py", "650", "family templates and linearizing transformations"],
     ["dimreflex/_ols.py / scenarios.py / anomaly.py / cli.py / monitoring.py", "1,735", "least squares, sampling scenarios, anomaly layer, CLI, monitoring"],
     ["form_check_vs_reset/common.py", "302", "four-test harness: A (permutation curvature), B (RESET-F), C (RESET-HC3), D (exog-HC3); failure registry"],
     ["form_check_vs_reset/scenarios.py", "319", "S1 synthetic scenario generators with frozen seed rules"],
     ["form_check_vs_reset/run_experiment.py / run_srsd.py / run_families.py", "619", "suite runners (S1–S3)"],
     ["form_check_vs_reset/battery.py / battery_probe.py / run_battery.py", "514", "family-viability battery, permutation-vs-wild probe, formal run"],
     ["form_check_vs_reset/srsd_data.py / analyze.py / srsd_analyze.py / families_analyze.py / battery_analyze.py / make_report.py", "1,339", "task registry and analysis scripts"],
     ["fc_multidirection/hybrid.py", "225", "HX: direction library, rank-aware selection, maxT + joint-F, Bonferroni combination"],
     ["fc_multidirection/multidir.py / run_bench.py / compare_vs_reset.py / random_ten.py", "965", "MD4 structural-only comparator, controlled and randomized runners, NIST challenge cases"],
     ["tests/ (11 files)", "1,936", "144 pytest test functions across core and experiment code"]],
    widths=[7.0, 1.8, 6.2], aligns=["l", "r", "l"])
cap("Table 17. Principal code modules and their roles (line counts measured on the release tree).")
h2("8.2 Data inventory")
p("Table 18 summarizes the released run data. All headline rates in this paper were computed from these artifacts; the per-cell counts permit recomputation of every Wilson interval and every calibrated-power number in Sections 4–5. The randomized stress-test data (1,000 groups) are released as a single CSV with per-group p-values for all four methods, the generated class label, the alternative family, and the winning HX direction, enabling exact reproduction of Tables 13–16.")
tbl(["Artifact", "Size", "Contents"],
    [["runs/full_500/raw.csv + raw.parquet + extrap.parquet", "53.3 MB + 16.5 MB + 5.2 MB", "S1: 142,500 fits, per-cell p-values and failure codes; 5,000 extrapolation probes"],
     ["runs/srsd_full/raw_srsd.parquet", "33.0 MB", "S2: 240,000 fits, 480 cells, min 500 seeds/cell; task registry"],
     ["runs/feynman_families/raw_families.csv + family_table_all120.csv", "1.5 MB", "S3: 11,760 fits, 120-equation family table, selection table"],
     ["runs/battery_full + battery_probe + battery_probe_wild", "337 KB + 23 KB + 8 KB", "S4: 94,000 battery datasets (zero data-layer failures); probe runs for Table 9"],
     ["runs/full_r0–r4 + srsd_full_r0–r4 + logs", "≈ 80 MB", "sharded raw outputs and run logs for audit"],
     ["fc_multidirection/runs/random_1000.csv", "134 KB", "1,000 randomized groups: p-values, class labels, families, winning directions"],
     ["fc_multidirection/runs/bench_summary.csv + argmax_directions.csv + compare_vs_reset.md", "217 KB", "controlled-benchmark summary, direction attribution, NIST comparisons"]],
    widths=[6.6, 3.4, 5.0], aligns=["l", "r", "l"])
cap("Table 18. Released run data (sizes measured on the release tree).")
h2("8.3 Test evidence")
p("The release carries 144 pytest test functions across 11 test files (Table 19). In addition to the project test suite, three independent verification passes were run on the released artifacts before this manuscript was assembled. First, every headline number in Tables 1–11 was re-derived from the released summary JSONs and compared digit-by-digit with the manuscript tables; the comparison matched at the reported precision (e.g., S2 HC3 false-positive rates .422/.417/.418 and the battery's wild-null repair .045/.090). Second, Tables 13–16 were recomputed from the released random_1000.csv using the rejection rule p ≤ 0.05 (the rule that reproduces the published counts 300/607 and 578/607; boundary p-values at exactly 0.050 account for the one-case differences under the alternative rule p < 0.05); all Wilson intervals, family-level rates, McNemar p-values, and direction shares reproduced. Third, the four-test harness itself was executed on fresh synthetic data (30 true-power-law datasets, n = 200): all four tests returned valid decisions with false-positive rates ≈ .03–.05, matching the S1 calibration conclusion of Table 1 and confirming that the HC3 miscalibration of Table 2 is a property of real physical designs, not of the implementation.")
tbl(["Test file", "Test functions", "Coverage"],
    [["test_families.py", "30", "family templates, linearizing transformations"],
     ["test_discover.py", "23", "discovery pipeline invariants"],
     ["test_anomaly.py", "20", "anomaly layer"],
     ["test_cli.py", "16", "command-line paths"],
     ["test_core.py", "12", "dimensional analysis core"],
     ["test_monitoring.py", "11", "monitoring"],
     ["test_autodiscovery.py", "10", "autodiscovery"],
     ["test_scenarios.py", "9", "scenario generators"],
     ["test_form_check.py", "7", "form-check statistic and null behavior"],
     ["test_uncertainty.py", "6", "uncertainty module"],
     ["Total", "144", "core + experiment layer"]],
    widths=[5.2, 3.0, 6.8], aligns=["l", "c", "l"])
cap("Table 19. Test suite inventory.")

# ---------------- 9 Conclusion ----------------
h1("9 Conclusion")
p("This paper treats the falsification side of physical model discovery end-to-end. Part I built the benchmark the field lacks: four suites, 488,260 fitted datasets, one protocol that puts four classical specification tests on a common basis, and full registration of every failure. Its findings are specific and actionable. Robust-covariance RESET variants — the field's default \u201csafe\u201d choice — are structurally miscalibrated on real physical laws (42–56% false alarms, flat in n), an inversion of their textbook reputation documented at a scale no single-dataset study could reach. Raw power rankings invert under matched false-positive rates, and invert again across task families. No test wins everywhere; the failure profiles are complementary; and the epistemic boundary is now priced: alarms enrich severe-extrapolation probability 2–3×, silence certifies nothing. The battery shows the way past the dead end that a rejection is: route, don't just refuse — but rebuild the null per family first, because a null that is exact in one coordinate system is silently wrong in another.")
p("Part II turned that map into a method. HX combines interpretable structural directions with fitted-value powers, removes redundant directions, evaluates concentrated and distributed departures, records numerical failure states, and — in the available randomized stress test — paired near-nominal null rejection (FPR 0.0509, Wilson CI 0.0332–0.0773) with power 0.9522 (0.9322–0.9665), substantially above a single-direction form check (0.4942), a structural-only library (0.8913), and classical RESET (0.7512), with paired discordance tests supporting the power difference at p < 10⁻¹¹. The result supports further evaluation, not unconditional deployment. The central methodological lesson is that a credible model check must expose its tested directions, its multiplicity correction, its resampling assumptions, and its failures — and that the benchmark and the method, read together, convert \u201cdid the pipeline find the right equation?\u201d into the falsifiable question \u201cwould a working falsification layer have caught it if it didn't?\u201d The code, data, and tests that answer that question for every claim in this paper are catalogued in Appendix D and released under fixed seed rules.")

# ---------------- References ----------------
h1("References")
refitem("[1] Ramsey, J. B. (1969). Tests for specification errors in classical linear least-squares regression analysis. Journal of the Royal Statistical Society, Series B, 31(2), 539–551.")
refitem("[2] White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. Econometrica, 48(4), 817–838.")
refitem("[3] MacKinnon, J. G., & White, H. (1985). Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties. Journal of Econometrics, 29(3), 305–325.")
refitem("[4] Long, J. S., & Ervin, L. H. (2000). Using heteroskedasticity consistent standard errors in the linear regression model. The American Statistician, 54(3), 217–224.")
refitem("[5] Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1, 292–298.")
refitem("[6] Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. NeuroImage, 92, 381–397.")
refitem("[7] Davidson, R., & Flachaire, E. (2008). The wild bootstrap, tamed at last. Journal of Econometrics, 146(1), 162–169.")
refitem("[8] Westfall, P. H., & Young, S. S. (1993). Resampling-Based Multiple Testing: Examples and Methods for P-Value Adjustment. Wiley.")
refitem("[9] Liu, R. Y. (1988). Bootstrap procedures under some non-i.i.d. models. Annals of Statistics, 16(4), 1696–1708.")
refitem("[10] La Cava, W., Orzechowski, P., Burlacu, B., et al. (2021). Contemporary symbolic regression methods and their relative performance. Proceedings of NeurIPS Datasets and Benchmarks Track.")
refitem("[11] de França, F. O., Virgolin, M., Kommenda, M., Majumder, M. S., Cranmer, M., Espada, G., Ingelse, L., Fonseca, A., Landajuela, M., Petersen, B. K., Glatt, R., Mundhenk, T. N., Lee, C. S., Hochhalter, J. D., Randall, D. L., Kamienny, P., Zhang, H., Dick, G., Simon, A., Burlacu, B., Kasak, J., Machado, M., Wilstrup, C., & La Cava, W. G. (2024). SRBench++: Principled benchmarking of symbolic regression with domain-expert interpretation. IEEE Transactions on Evolutionary Computation. https://doi.org/10.1109/TEVC.2024.3423681")
refitem("[12] Mundhenk, T. N., Landajuela, M., Glatt, R., et al. (2022). Rethinking symbolic regression datasets and benchmarks for scientific discovery. arXiv:2206.10540.")
refitem("[13] Udrescu, S.-M., & Tegmark, M. (2020). AI Feynman: A physics-inspired method for symbolic regression. Science Advances, 6(16), eaay2631.")
refitem("[14] Bakarji, J., Callaham, J. L., Brunton, S. L., & Kutz, J. N. (2022). Dimensionally consistent learning with Buckingham Pi. Nature Computational Science, 2, 834–844.")
refitem("[15] Cranmer, M. (2023). Interpretable machine learning for science with PySR and SymbolicRegression.jl. arXiv:2305.01582.")
refitem("[16] Fisher, R. A. (1935). The Design of Experiments. Oliver & Boyd.")

# ---------------- Appendix A ----------------
br()
h1("Appendix A — Suite Composition and Algorithm Specification")
h2("A.1 Suite composition")
tbl(["Suite", "Cells", "Fits per cell", "Total fits", "Extra evaluations"],
    [["S1 synthetic", "19 × 3n × 5σ = 285", "500", "142,500", "—"],
     ["S2 SRSD-Feynman", "32 × 3n × 5σ = 480", "500", "240,000", "240,000 same-distribution test sets; 232,500 extrapolation probes (7,500 not finite on probe domain, registered)"],
     ["S3 family-stratified", "392", "30", "11,760", "—"],
     ["S4 battery", "188 (76 synth + 112 SRSD)", "500", "94,000", "—"],
     ["Total", "1,345", "—", "488,260", "472,500"]],
    widths=[3.4, 4.4, 2.4, 2.4, 6.4], aligns=["l", "l", "c", "r", "l"])
cap("Table A.1. Grand totals across the four suites.")
h2("A.2 HX algorithm (declared resampling mode)")
b("1. Validate that Xℓ is finite and two-dimensional, ℓy is finite and aligned, n exceeds the residual-degree floor, and requested parameters are valid.")
b("2. Standardize columns of Xℓ, form B = [1, Z], check rank and conditioning, and fit the restricted model.")
b("3. Compute ŷ and r; reject constant-response or degenerate-exact-fit cases using TSS-based criteria.")
b("4. Construct named structural directions and standardized fitted-value powers.")
b("5. Residualize candidates against B without materializing an n × n projector.")
b("6. Use modified Gram–Schmidt with deterministic, scale-aware pivoting to retain q independent coordinates W.")
b("7. Compute observed Tmax and Fjoint; record the winning direction as exploratory.")
b("8. Generate R resamples under the declared mode; for confirmatory response-adaptive inference, rebuild every adaptive step inside each resample or use a split-sample construction.")
b("9. Compute pmax and pF with the plus-one correction, then pHX = min(1, 2·min(pmax, pF)).")
b("10. Return the complete result object or a named failure state; never return a silent NaN or an unclassified exception.")

# ---------------- Appendix B ----------------
br()
h1("Appendix B — Failure Registration Taxonomy")
tbl(["Reason", "Meaning", "Exercised by"],
    [["degenerate_zero_residual", "Zero-noise cell on exact ground truth; no residual variation → test undefined", "S1/S2 σ = 0 cells; battery degenerate draws"],
     ["rank_deficient", "Design matrix numerically rank-deficient", "collinear designs"],
     ["insufficient_df", "n − k < 1 after augmentation", "n = 50 multi-feature cells"],
     ["hc3_leverage_degenerate", "1 − hᵢ ≤ 10⁻¹²; HC3 undefined", "degenerate built columns"],
     ["hc3_singular_covariance", "V_HC3 numerically singular", "degenerate built columns"],
     ["statistic_mismatch", "Wrapper statistic ≠ native statistic (drift sentinel)", "continuous audit"],
     ["checked_false_*", "Native test declines to run (non-positive data, n < 12, non-finite)", "not-applicable tasks"],
     ["nonpositive / transform_failed / transform_nonfinite / degenerate_null", "Battery member precondition failures", "battery members on edge-case data"],
     ["bad_shape_X / length_mismatch / nonfinite_input / domain_error", "HX input and domain guards", "HX guard conditions"],
     ["rank_deficient_base / ill_conditioned_base / constant_response / no_valid_direction / insufficient_residual_df", "HX base-model and augmentation guards", "HX guard conditions"],
     ["invalid_exchangeability / resampling_failure", "Declared null generator unjustified or unexecutable", "HX resampling guards"]],
    widths=[5.2, 5.8, 4.0], aligns=["l", "l", "l"])
cap("Table B.1. Failure reasons registered across the benchmark harness and HX, with the paths that exercise them.")
p("All registered failures remain in the audit trail; rates' denominators are documented per table. Zero-noise cells intentionally trigger degenerate_zero_residual for all four methods — the registration path is tested by the benchmark itself. The battery's formal run registered zero data-layer failures across 94,000 datasets.")

# ---------------- Appendix C ----------------
br()
h1("Appendix C — Headline Wilson 95% Intervals")
p("S2 H0, n = 200 (26,000 fits per cell): A .053 [.050–.055]; B .050 [.047–.052]; C .417 [.411–.423]; D .563 [.557–.569]. S1 H0, n = 200 (4,500 fits per cell): A .037 [.032–.043]; B .041 [.036–.047]; C .045 [.039–.051]; D .042 [.036–.048]. Battery power-family size (S4, 16,000 member tests): .0496 [.046–.053]. Randomized stress test (393 null / 607 alternative): HX FPR 0.0509 [0.0332, 0.0773], power 0.9522 [0.9322, 0.9665]. Intervals for all other cells are computable from the released per-cell counts.")

# ---------------- Appendix D ----------------
br()
h1("Appendix D — Code, Data, and Test Inventory")
p("This appendix is the complete release manifest for the claims made in this paper. The release tree contains the core package, two experiment trees, all run data, and the test suite described in Section 8. Every number in Tables 1–16 is computed from the files below; the analysis scripts in the two experiment trees regenerate the summary JSONs and CSV tables from the raw parquet/CSV artifacts, so any claim is one command away from its evidence.")
h2("D.1 Core package (dimreflex/)")
tbl(["File", "Lines", "Role"],
    [["__init__.py", "89", "package surface, version, exports"],
     ["core.py", "438", "dimensional analysis, unit registry, form-check statistic (the primitive under test in Part I)"],
     ["discover.py", "518", "monomial discovery, candidate surface construction, screening"],
     ["families.py", "650", "family templates, linearizing transformations, registry"],
     ["_ols.py", "85", "least-squares kernel with rank handling"],
     ["autodiscovery.py", "322", "automatic pipeline orchestration"],
     ["uncertainty.py", "167", "uncertainty propagation"],
     ["anomaly.py", "463", "anomaly detection layer"],
     ["scenarios.py", "383", "sampling scenarios and seed rules"],
     ["templates.py", "141", "report templates"],
     ["cli.py", "326", "command-line interface"],
     ["monitoring.py", "477", "run monitoring and content hashing"],
     ["launcher.py", "22", "launcher entry"]],
    widths=[6.0, 1.8, 7.2], aligns=["l", "r", "l"])
cap("Table D.1. Core package modules (total 4,081 lines).")
h2("D.2 Benchmark experiment tree (experiments/form_check_vs_reset/)")
tbl(["File", "Lines", "Role"],
    [["common.py", "302", "four-test harness (A/B/C/D), permutation and wild nulls, failure registry, environment snapshot"],
     ["scenarios.py", "319", "S1 synthetic scenario generators, frozen seed rule"],
     ["run_experiment.py", "187", "S1 runner"],
     ["run_srsd.py", "184", "S2 runner over the SRSD-Feynman registry"],
     ["run_families.py", "248", "S3 family-stratified runner"],
     ["run_battery.py", "170", "S4 battery formal runner (500 seeds)"],
     ["battery.py", "132", "battery members and viability semantics"],
     ["battery_probe.py", "212", "permutation-vs-wild probe (Table 9)"],
     ["battery_probe2.py", "147", "extended probe"],
     ["srsd_data.py", "217", "SRSD task registry and official sampling objects"],
     ["analyze.py / srsd_analyze.py / families_analyze.py / battery_analyze.py", "910", "summary generation, Wilson intervals, calibrated power"],
     ["make_report.py", "212", "markdown report assembly"]],
    widths=[6.4, 1.8, 6.8], aligns=["l", "r", "l"])
cap("Table D.2. Benchmark experiment tree.")
h2("D.3 Method experiment tree (experiments/fc_multidirection/)")
tbl(["File", "Lines", "Role"],
    [["hybrid.py", "225", "HX: direction library, rank-aware modified Gram–Schmidt, maxT + joint-F, Bonferroni combination, failure states"],
     ["multidir.py", "219", "MD4 structural-only comparator (maxT over structural directions)"],
     ["run_bench.py", "260", "controlled development benchmark runner"],
     ["compare_vs_reset.py", "199", "NIST challenge cases and RESET comparison"],
     ["random_ten.py", "287", "randomized stress-test generator (master seed 20261904, 1,000 groups)"]],
    widths=[6.0, 1.8, 7.2], aligns=["l", "r", "l"])
cap("Table D.3. Method experiment tree (HX development).")
h2("D.4 Run data")
tbl(["Path", "Size", "Contents"],
    [["runs/full_500/raw.csv", "53.3 MB", "S1 raw per-cell results (142,500 fits)"],
     ["runs/full_500/raw.parquet", "16.5 MB", "S1 parquet mirror"],
     ["runs/full_500/extrap.parquet", "5.2 MB", "S1/S2 extrapolation probes"],
     ["runs/full_500/cell_stats.csv", "101 KB", "per-cell counts and failure registrations"],
     ["runs/full_500/summary.json, config.json, env.json, report.md, win_matrix.csv", "22 KB", "S1 summaries, frozen configuration, environment hash, report"],
     ["runs/srsd_full/raw_srsd.parquet", "33.0 MB", "S2 raw results (240,000 fits)"],
     ["runs/srsd_full/task_registry.csv", "44 KB", "SRSD task registry with per-task metadata"],
     ["runs/srsd_full/summary_srsd.json, h1_power_per_task.csv, config.json, env.json", "24 KB", "S2 summaries and configuration"],
     ["runs/feynman_families/raw_families.csv", "1.5 MB", "S3 raw results (11,760 fits)"],
     ["runs/feynman_families/family_table_all120.csv, per_equation.csv, selection_table.csv, summary_families.json, report_families.md", "67 KB", "S3 family tables, selection, summaries"],
     ["runs/battery_full/raw_battery.parquet + summary_battery.json + report_battery.md", "352 KB", "S4 formal run (94,000 datasets, zero data-layer failures)"],
     ["runs/battery_probe/ + battery_probe_wild/", "33 KB", "probe runs for Table 9 (permutation vs wild)"],
     ["runs/full_r0–r4/, srsd_full_r0–r4/, full/, pilot/, pilot_v2/, _smoke_fam/", "≈ 80 MB", "sharded and pilot runs for audit"],
     ["fc_multidirection/runs/random_1000.csv", "134 KB", "1,000 randomized groups: p_fc, p_md4, p_hx, p_reset, class labels, families, argmax_hx"],
     ["fc_multidirection/runs/bench_raw.csv + bench_summary.csv + argmax_directions.csv", "217 KB", "controlled benchmark raw and summary"],
     ["fc_multidirection/runs/compare_vs_reset.md + compare_nist_md4.csv + random_10.* + report.md", "14 KB", "NIST challenge results and development reports"]],
    widths=[7.6, 2.4, 5.0], aligns=["l", "r", "l"])
cap("Table D.4. Run data artifacts released with this paper.")
h2("D.5 Test evidence")
p("The release carries 144 pytest test functions across 11 test files (Table 19). Three independent verification passes were run on the released artifacts before assembly of this manuscript, all reproducible from the release: (i) digit-by-digit re-derivation of Tables 1–11 from the released summary JSONs; (ii) recomputation of Tables 13–16 from random_1000.csv under the rejection rule p ≤ 0.05 (the rule reproducing the published counts; see Section 8.3); (iii) fresh execution of the four-test harness on 30 synthetic true-power-law datasets (n = 200), yielding false-positive rates ≈ .03–.05 consistent with Table 1 and confirming that the HC3 miscalibration of Table 2 is design-driven rather than implementation-driven. These passes are themselves part of the released record, and the analysis scripts used are listed in Tables D.2–D.3.")

# =================== MAIN ===================
first_toc = False
for kind, payload in C:
    if kind == "h1": H(1, payload)
    elif kind == "h2": H(2, payload)
    elif kind == "h3": H(3, payload)
    elif kind == "p": P(payload)
    elif kind == "b": BULLET(payload)
    elif kind == "cap": CAPTION(payload)
    elif kind == "tbl":
        headers, rows, widths, aligns = payload
        TABLE(headers, rows, widths, aligns)
    elif kind == "br": PAGE_BREAK()
    elif kind == "toc":
        TOC_BLOCK()
        PAGE_BREAK()
    elif kind == "ref":
        p = doc.add_paragraph()
        p.paragraph_format.first_line_indent = Pt(0)
        p.paragraph_format.left_indent = Cm(1.0)
        p.paragraph_format.first_line_indent = Cm(-1.0)  # hanging indent
        p.paragraph_format.line_spacing = 1.5
        p.add_run(payload)

doc.save(OUT)
print("saved:", OUT, os.path.getsize(OUT), "bytes")
