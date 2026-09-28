import sys
sys.path.insert(0, "D:/form_check")
import srsd_data
import pandas as pd
from sympy import sympify, Add, Symbol, Float, Mul
from sympy.functions.elementary.trigonometric import sin, cos, tan

reg = srsd_data.load_feynman()

def to_class_name(task):
    if task.startswith("bonus"):
        return "FeynmanBonus" + task.split(".")[1]
    vol, ch, eq = task.split(".")
    roman = {"i": "I", "ii": "II", "iii": "III"}[vol]
    return f"Feynman{roman}Ch{ch}Eq{eq}"

def is_tiny(node):
    # any Float coefficient under 1e-5 anywhere in this subexpression
    fl = node.atoms(Float)
    return any(abs(float(f)) < 1e-5 for f in fl)

def classify(expr):
    # 1. trig terms with O(1) coefficient -> visible oscillating curvature
    for t in expr.atoms(sin, cos, tan):
        # walk up not possible in sympy; instead: check each Mul containing a trig has non-tiny
        # numeric coefficient OR appears as a bare arg. Simplest honest proxy: if the trig's
        # enclosing multiplicative coefficient includes a tiny float -> suppressed.
        pass
    trigs = expr.atoms(sin, cos, tan)
    if trigs:
        # find the minimal subtrees containing each trig and a coefficient
        suppressed = True
        for t in trigs:
            found = False
            for a in expr.atoms(Add):
                for arg in a.args:
                    if t in arg.free_symbols or t in arg.atoms(sin, cos, tan):
                        if not is_tiny(arg):
                            suppressed = False
            # trig added directly at top level with no Add -> O(1)
        # simpler: if EVERY Add-arg containing trig is tiny-coefficiented, trig is suppressed
        has_visible = False
        for t in trigs:
            in_any_arg = False
            for a in expr.atoms(Add):
                for arg in a.args:
                    if arg.has(t):
                        if not is_tiny(arg):
                            has_visible = True
            # if trig not inside any Add at all, it is structural -> visible
            if not any(arg.has(t) for a in expr.atoms(Add) for arg in a.args):
                has_visible = True
        if has_visible:
            return False, "trig_O1"
        return True, "suppressed_trig"
    # 2. cross-term
    for a in expr.atoms(Add):
        for arg in a.args:
            if len(arg.atoms(Symbol)) >= 2:
                return True, "cross_term"
    # 3. mild correction
    if any(abs(float(f)) < 1e-5 for f in expr.atoms(Float)):
        return True, "small_correction"
    return False, "plain"

tasks = ["iii.4.33","ii.11.28","iii.4.32","ii.11.27","i.41.16","bonus.15",
         "i.10.7","i.48.2","ii.13.23","ii.13.34","bonus.5","i.27.6","i.37.4",
         "i.39.11","i.43.43","bonus.20","i.26.2","ii.38.14"]
measured = {"iii.4.33":.849,"ii.11.28":.830,"iii.4.32":.827,"ii.11.27":.813,
            "i.41.16":.797,"bonus.15":.734,"i.10.7":.385,"i.48.2":.373,
            "ii.13.23":.373,"ii.13.34":.355,"bonus.5":.164,"i.27.6":.162,
            "i.37.4":.090,"i.39.11":.030,"i.43.43":.029,"bonus.20":.000,
            "i.26.2":.000,"ii.38.14":.000}
rows = []
for t in tasks:
    e = sympify(str(reg[to_class_name(t)]().sympy_eq))
    blind, why = classify(e)
    rows.append((t, blind, why, measured[t]))
d = pd.DataFrame(rows, columns=["task","pred_blind","rule","measured_silence"])
d["hit"] = (d.pred_blind == (d.measured_silence > .5))
print(d.to_string(index=False))
print("\nhits:", d.hit.sum(), "/", len(d), "=", round(d.hit.mean(),3))
from scipy.stats import pointbiserialr
rr = pointbiserialr(d.pred_blind.astype(int), d.measured_silence)
print("point-biserial r:", round(rr.statistic,3), " p:", round(rr.pvalue,5))
tp = int(((d.pred_blind) & (d.measured_silence>.5)).sum()); fp = int(((d.pred_blind) & (d.measured_silence<=.5)).sum())
fn = int(((~d.pred_blind) & (d.measured_silence>.5)).sum()); tn = int(((~d.pred_blind) & (d.measured_silence<=.5)).sum())
print(f"TP={tp} FP={fp} FN={fn} TN={tn}")
d.to_csv("D:/form_check/runs/joint_silence/blind_test_per_task.csv", index=False)
print("written")
