# ============================================================================
#  fig1_evidence_base.py  —  Section-1 figure (Python / stdlib SVG)
#     "A young, low-quality, redundant and uneven evidence base"
#
#  Four panels, all with axes (a scientific figure, not an infographic):
#   a  review quality vs publication year (OLS fit + 95% CI band; no improvement)
#   b  evidence coverage: plant organ × NCP category, coloured by primary studies
#   c  redundancy: distribution of primary studies shared between review pairs
#   d  evidence depth: distribution of primary studies per trait–NCP relationship
#
#  INPUT   data/data_model.csv, data/meta-data2.csv, data/study_scores.csv,
#          data/sigma.csv
#  OUTPUT  outputs/fig1_evidence_base.svg (+ .png if cairosvg installed)
#  RUN     python3 figures/fig1_evidence_base.py
#
#  ADDED:  4 seperate figures for poster presentation (see below)
# ============================================================================
import csv, os, sys, math, statistics as st
from collections import defaultdict

#HERE = os.path.dirname(os.path.abspath(__file__))
HERE = os.getcwd() #for RStudio environment
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == "figures" else HERE
D = lambda f: os.path.join(ROOT, "data", f)
DM, META, SS, SIG = D("data_model.csv"), D("meta-data2.csv"), D("study_scores.csv"), D("sigma.csv")
OUT = os.path.join(ROOT, "outputs", "fig1_evidence_base.svg")
for p in (DM, META, SS, SIG):
    if not os.path.exists(p): sys.exit("Missing %s" % p)

dm = list(csv.DictReader(open(DM)))
def num(x):
    try: return float(x)
    except Exception: return 0.0
revs = sorted(set(r["source_id"] for r in dm))

# meta (first row per review), quality
meta = {}
for r in csv.DictReader(open(META, encoding="latin1"), delimiter=";"):
    sid = r["source_id"].strip()
    if sid not in meta: meta[sid] = r
ss = {r["source_id"]: float(r["prop_score"]) for r in csv.DictReader(open(SS))}

# a) quality vs year
qy = []
for s in revs:
    try: qy.append((int(float(meta[s]["year"])), ss[s]))
    except Exception: pass

# b) redundancy: shared studies per review pair (recovered from sigma)
rows = list(csv.reader(open(SIG))); mat = [[float(x) for x in r] for r in rows[1:]]; nrev = len(rows[0])
nj = [1 / mat[i][i] if mat[i][i] > 0 else 0 for i in range(nrev)]
shared = [round(mat[i][j] * nj[i] * nj[j]) for i in range(nrev) for j in range(i + 1, nrev)]
sb = {"0": 0, "1": 0, "2": 0, "3–5": 0, "6+": 0}
for x in shared:
    sb["0" if x == 0 else "1" if x == 1 else "2" if x == 2 else "3–5" if x <= 5 else "6+"] += 1
SPAIRS = len(shared)

# c) coverage organ × ncp-category
organs = ["leaf", "whole plant", "root", "stem", "flower", "reproductive", "litter", "life history"]
cats = ["regulation", "material", "non-material"]
cov = defaultdict(float)
for r in dm: cov[(r["trait_category"], r["ncp_category"])] += num(r["ntot"])

# d) depth: studies per relationship
ps = defaultdict(float)
for r in dm: ps[r["ts"]] += num(r["ntot"])
db = {"1": 0, "2": 0, "3–5": 0, "6–10": 0, ">10": 0}
for n in ps.values():
    db["1" if n <= 1 else "2" if n == 2 else "3–5" if n <= 5 else "6–10" if n <= 10 else ">10"] += 1
NPAIR = len(ps)

def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
F = "font-family='system-ui,-apple-system,\"Segoe UI\",sans-serif'"
BLUE, MUTE, INK, SUB, ORANGE = "#2a78d6", "#c8cbd1", "#0b0b0b", "#52514e", "#eb6834"
W, H = 1080, 820
PX, PY, PW, PH = [60, 610], [135, 490], 410, 225
s = [f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">',
     f'<rect width="{W}" height="{H}" fill="#ffffff"/>'#, #here change color to either #ffffff for a full white background or #fcfcfb for a beige/kind greyish background
     #f'<text x="44" y="46" {F} font-size="21" font-weight="800" fill="{INK}">A young, low-quality, redundant and uneven evidence base</text>',
     #f'<text x="44" y="70" {F} font-size="12.5" fill="{SUB}">22 systematic reviews · 474 trait–NCP relationships · 171 traits × 13 NCPs · 2,203 primary-study outcomes.</text>'
     ]
def axes(x, y, w, h):
    s.append(f'<line x1="{x}" y1="{y+h}" x2="{x+w}" y2="{y+h}" stroke="#9a9a94" stroke-width="1.2"/>')
    s.append(f'<line x1="{x}" y1="{y}" x2="{x}" y2="{y+h}" stroke="#9a9a94" stroke-width="1.2"/>')

# ---- a) quality vs year (OLS + 95% CI band) ----
x0, y0 = PX[0], PY[0]
s.append(f'<text x="{x0}" y="{y0-12}" {F} font-size="13.5" font-weight="700" fill="{INK}">a Reviews are not becoming more rigorous</text>')
yrs = [a for a, b in qy]
axes(x0, y0, PW, PH)
ax = lambda v: x0 + PW * (v - 2017) / 7
ay = lambda v: y0 + PH * (1 - v / 0.8)
for g in (0, .2, .4, .6, .8):
    s.append(f'<line x1="{x0}" y1="{ay(g):.0f}" x2="{x0+PW}" y2="{ay(g):.0f}" stroke="#f0efe9"/>')
    s.append(f'<text x="{x0-8}" y="{ay(g)+3:.0f}" {F} font-size="9" fill="{SUB}" text-anchor="end">{g:.1f}</text>')
for yr in range(2017, 2025):
    s.append(f'<text x="{ax(yr):.0f}" y="{y0+PH+16}" {F} font-size="8.5" fill="{SUB}" text-anchor="middle">{yr}</text>')
mx = st.mean(yrs); my = st.mean([b for a, b in qy]); nn = len(qy)
sxx = sum((a - mx) ** 2 for a in yrs); b1 = sum((a - mx) * (b - my) for a, b in qy) / sxx; b0 = my - b1 * mx
sse = sum((b - (b0 + b1 * a)) ** 2 for a, b in qy); s2 = sse / (nn - 2); tc = 2.086
top, bot = [], []
for gx in [2017 + .25 * k for k in range(29)]:
    fit = b0 + b1 * gx; sef = math.sqrt(s2 * (1 / nn + (gx - mx) ** 2 / sxx))
    top.append((gx, fit + tc * sef)); bot.append((gx, fit - tc * sef))
poly = " ".join(f"{ax(gx):.1f},{ay(min(.8,max(0,v))):.1f}" for gx, v in top)
poly += " " + " ".join(f"{ax(gx):.1f},{ay(min(.8,max(0,v))):.1f}" for gx, v in reversed(bot))
s.append(f'<polygon points="{poly}" fill="{ORANGE}" fill-opacity="0.10"/>')
s.append(f'<line x1="{ax(2017):.0f}" y1="{ay(b0+b1*2017):.0f}" x2="{ax(2024):.0f}" y2="{ay(b0+b1*2024):.0f}" stroke="{ORANGE}" stroke-width="1.8"/>')
for a, b in qy:
    s.append(f'<circle cx="{ax(a):.0f}" cy="{ay(b):.0f}" r="3.6" fill="none" stroke="{BLUE}" stroke-width="1.5"/>')
ci_lo, ci_hi = b1 - tc * math.sqrt(s2 / sxx), b1 + tc * math.sqrt(s2 / sxx)
tval = b1 / math.sqrt(s2 / sxx)
# two-sided p via a small t-CDF (Hill's approximation), df = nn-2
def tp(t, df):
    t = abs(t); a = df / 2; x = df / (df + t * t)
    # incomplete beta via continued fraction
    def betacf(x, a, b):
        MAXIT = 200; EPS = 3e-12; qab = a + b; qap = a + 1; qam = a - 1
        c = 1; d = 1 - qab * x / qap; d = 1 / (d if abs(d) > 1e-30 else 1e-30); h = d
        for m in range(1, MAXIT):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d
            d = 1 / (d if abs(d) > 1e-30 else 1e-30); c = 1 + aa / (c if abs(c) > 1e-30 else 1e-30); h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d
            d = 1 / (d if abs(d) > 1e-30 else 1e-30); c = 1 + aa / (c if abs(c) > 1e-30 else 1e-30)
            de = d * c; h *= de
            if abs(de - 1) < EPS: break
        return h
    bt = math.exp(math.lgamma(a + b if False else 0) ) if False else 0
    lb = math.lgamma(a + 0.5) - math.lgamma(a) - math.lgamma(0.5) + a * math.log(x) + 0.5 * math.log(1 - x)
    ib = math.exp(lb) * betacf(x, a, 0.5) / a
    return ib  # two-sided p
pval = tp(tval, nn - 2)
#s.append(f'<text x="{x0+PW-6:.0f}" y="{y0+15}" {F} font-size="10" fill="{INK}" text-anchor="end">no improvement with year</text>')
s.append(f'<text x="{x0+PW-6:.0f}" y="{y0+30}" {F} font-size="9.5" fill="{SUB}" text-anchor="end">β = {b1:.3f} yr⁻¹ (95% CI {ci_lo:.3f}, {ci_hi:.3f}); p = {pval:.2f}; n = {nn}</text>')
s.append(f'<text x="{x0-42}" y="{y0+PH/2:.0f}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle" transform="rotate(-90 {x0-42} {y0+PH/2:.0f})">review quality score</text>')

# ---- b) redundancy ----
x0, y0 = PX[1], PY[0]
s.append(f'<text x="{x0}" y="{y0-12}" {F} font-size="13.5" font-weight="700" fill="{INK}">b  Few primary studies are shared between reviews</text>')
axes(x0, y0, PW, PH); order = ["0", "1", "2", "3–5", "6+"]; mxs = max(sb.values()); bw = PW / len(order)
for i, k in enumerate(order):
    v = sb[k]; bh = 0.9 * PH * v / mxs; col = BLUE if k == "0" else MUTE
    s.append(f'<rect x="{x0+i*bw+bw*0.16:.0f}" y="{y0+PH-bh:.0f}" width="{bw*0.68:.0f}" height="{bh:.0f}" fill="{col}"/>')
    s.append(f'<text x="{x0+i*bw+bw/2:.0f}" y="{y0+PH-bh-6:.0f}" {F} font-size="10" font-weight="{"700" if k=="0" else "400"}" fill="{INK if k=="0" else SUB}" text-anchor="middle">{round(100*v/SPAIRS)}%</text>')
    s.append(f'<text x="{x0+i*bw+bw/2:.0f}" y="{y0+PH+16}" {F} font-size="9.5" fill="{SUB}" text-anchor="middle">{esc(k)}</text>')
s.append(f'<text x="{x0+PW/2:.0f}" y="{y0+PH+32}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle">primary studies shared between a pair of reviews (n = {SPAIRS} pairs)</text>')
s.append(f'<text x="{x0-42}" y="{y0+PH/2:.0f}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle" transform="rotate(-90 {x0-42} {y0+PH/2:.0f})">% of review pairs</text>')
#s.append(f'<text x="{x0+PW-6:.0f}" y="{y0+16}" {F} font-size="10" font-weight="700" fill="{BLUE}" text-anchor="end">{round(100*sb["0"]/SPAIRS)}% share none · 10.7% CCA (overall redundancy)</text>')

# ---- c) coverage heatmap ----
x0, y0 = PX[0], PY[1]+12
s.append(f'<text x="{x0}" y="{y0-24}" {F} font-size="13.5" font-weight="700" fill="{INK}">c  Evidence concentrated on leaf × regulating</text>')
olab = {"leaf":"Leaf","whole plant":"Whole-plant","root":"Root","stem":"Stem","flower":"Flower","reproductive":"Reprod.","litter":"Litter","life history":"Life hist."}
clab = {"regulation":"Regulating","material":"Material","non-material":"Non-mat."}
maxv = max(cov.values()); cw = (PW - 120) / len(cats); ch = (PH-10) / len(organs); gx0 = x0 + 120
for j, cat in enumerate(cats):
    s.append(f'<text x="{gx0+j*cw+cw/2:.0f}" y="{y0-6}" {F} font-size="9.5" fill="{SUB}" text-anchor="middle">{clab[cat]}</text>')
for i, o in enumerate(organs):
    s.append(f'<text x="{gx0-8}" y="{y0+i*ch+ch/2+3:.0f}" {F} font-size="9.5" fill="{INK}" text-anchor="end">{olab[o]}</text>')
    for j, cat in enumerate(cats):
        v = cov.get((o, cat), 0); t = math.log1p(v) / math.log1p(maxv) if maxv else 0
        col = "#%02x%02x%02x" % (int(255-(255-42)*t), int(255-(255-120)*t), int(255-(255-214)*t)) if v > 0 else "#f6f5f0"
        cx, cy = gx0 + j * cw, y0 + i * ch
        s.append(f'<rect x="{cx:.0f}" y="{cy:.0f}" width="{cw-2:.0f}" height="{ch-2:.0f}" fill="{col}" stroke="#fff" stroke-width="1"/>')
        if v > 0:
            text_col = "#ffffff" if t > 0.75 else "#000000"
            s.append(f'<text x="{cx+cw/2:.0f}" y="{cy+ch/2+3:.0f}" {F} font-size="8.5" fill="{text_col}" text-anchor="middle">{int(v)}</text>')
s.append(f'<text x="{gx0+(PW-120)/2:.0f}" y="{y0+PH+20}" {F} font-size="9.5" fill="{SUB}" text-anchor="middle">colour = primary-study outcomes; many combinations have none</text>')

# ---- d) depth ----
x0, y0 = PX[1], PY[1]
s.append(f'<text x="{x0}" y="{y0-12}" {F} font-size="13.5" font-weight="700" fill="{INK}">d  Many relationships rely on a single study</text>')
axes(x0, y0, PW, PH); o2 = ["1", "2", "3–5", "6–10", ">10"]; mx2 = max(db.values()); bw = PW / len(o2)
for i, k in enumerate(o2):
    v = db[k]; bh = 0.9 * PH * v / mx2; col = BLUE if k == "1" else MUTE
    s.append(f'<rect x="{x0+i*bw+bw*0.16:.0f}" y="{y0+PH-bh:.0f}" width="{bw*0.68:.0f}" height="{bh:.0f}" fill="{col}"/>')
    s.append(f'<text x="{x0+i*bw+bw/2:.0f}" y="{y0+PH-bh-6:.0f}" {F} font-size="10" font-weight="{"700" if k=="1" else "400"}" fill="{INK if k=="1" else SUB}" text-anchor="middle">{round(100*v/NPAIR)}%</text>')
    s.append(f'<text x="{x0+i*bw+bw/2:.0f}" y="{y0+PH+16}" {F} font-size="9.5" fill="{SUB}" text-anchor="middle">{esc(k)}</text>')
s.append(f'<text x="{x0+PW/2:.0f}" y="{y0+PH+32}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle">primary-study outcomes per trait–NCP relationship (n = {NPAIR})</text>')
s.append(f'<text x="{x0-42}" y="{y0+PH/2:.0f}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle" transform="rotate(-90 {x0-42} {y0+PH/2:.0f})">% of relationships</text>')

#s.append(f'<line x1="44" y1="{H-52}" x2="{W-44}" y2="{H-52}" stroke="#ecebe5"/>')
#s.append(f'<text x="44" y="{H-30}" {F} font-size="12" font-weight="600" fill="{INK}">The field is expanding fast without consolidating: reviews are low-quality (and no better over time), largely independent,</text>')
#s.append(f'<text x="44" y="{H-12}" {F} font-size="12" font-weight="600" fill="{INK}">skewed toward leaf traits and regulating services, and — relationship by relationship — very thinly evidenced.</text>')
s.append('</svg>')
open(OUT, "w").write("\n".join(s))
try:
    import cairosvg
    cairosvg.svg2png(url=OUT, write_to=OUT.replace(".svg", ".png"), output_width=1450)
    png = " (+ PNG)"
except Exception:
    png = " (SVG only; pip install cairosvg for PNG)"
print("wrote %s%s" % (OUT, png))

# or
#Inkscape is required

import subprocess
import os

png_file = OUT.replace(".svg", ".png")

try:
    subprocess.run(
        [
            r"C:\Program Files\Inkscape\bin\inkscape.exe",
            OUT,
            f"--export-filename={png_file}"
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    png = " (+ PNG)"
except Exception as e:
    print("Erreur Inkscape :", e)
    png = " (SVG only)"

print("wrote %s%s" % (OUT, png))

# ============================================================================
# fig1_evidence_base_individual.py
# Generate the four panels of Fig. 1 as INDIVIDUAL SVG files,
# then export each SVG to PNG using Inkscape.
#
# INPUT:
#   data/data_model.csv
#   data/meta-data2.csv
#   data/study_scores.csv
#   data/sigma.csv
#
# OUTPUT:
#   outputs/fig1a_quality_year.svg
#   outputs/fig1a_quality_year.png
#   outputs/fig1b_redundancy.svg
#   outputs/fig1b_redundancy.png
#   outputs/fig1c_coverage.svg
#   outputs/fig1c_coverage.png
#   outputs/fig1d_depth.svg
#   outputs/fig1d_depth.png
#
# RUN:
#   python3 figures/fig1_evidence_base_individual.py
# ============================================================================

import csv
import os
import sys
import math
import statistics as st
import subprocess
from collections import defaultdict


# ============================================================================
# PATHS
# ============================================================================
HERE = os.getcwd()
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == "figures" else HERE

D = lambda f: os.path.join(ROOT, "data", f)

DM = D("data_model.csv")
META = D("meta-data2.csv")
SS = D("study_scores.csv")
SIG = D("sigma.csv")

OUTDIR = os.path.join(ROOT, "outputs")
os.makedirs(OUTDIR, exist_ok=True)

for p in (DM, META, SS, SIG):
    if not os.path.exists(p):
        sys.exit("Missing: %s" % p)


# ============================================================================
# INKSCAPE
# ============================================================================

INKSCAPE = r"C:\Program Files\Inkscape\bin\inkscape.exe"

if not os.path.exists(INKSCAPE):
    print("WARNING: Inkscape not found at:")
    print(INKSCAPE)
    print("SVG files will still be generated.")
    print()


# ============================================================================
# DATA
# ============================================================================

dm = list(csv.DictReader(open(DM, encoding="utf-8")))


def num(x):
    try:
        return float(x)
    except Exception:
        return 0.0


revs = sorted(set(r["source_id"] for r in dm))


# ----------------------------------------------------------------------------
# Meta-data: first row per review
# ----------------------------------------------------------------------------

meta = {}

for r in csv.DictReader(
    open(META, encoding="latin1"),
    delimiter=";"
):
    sid = r["source_id"].strip()

    if sid not in meta:
        meta[sid] = r


# ----------------------------------------------------------------------------
# Review quality scores
# ----------------------------------------------------------------------------

ss = {
    r["source_id"]: float(r["prop_score"])
    for r in csv.DictReader(open(SS, encoding="utf-8"))
}


# ============================================================================
# A) QUALITY VS YEAR
# ============================================================================

qy = []

for s in revs:
    try:
        qy.append(
            (
                int(float(meta[s]["year"])),
                ss[s]
            )
        )
    except Exception:
        pass


# ============================================================================
# B) REDUNDANCY
# ============================================================================

rows = list(
    csv.reader(
        open(SIG, encoding="utf-8")
    )
)

mat = [
    [float(x) for x in r]
    for r in rows[1:]
]

nrev = len(rows[0])

nj = [
    1 / mat[i][i] if mat[i][i] > 0 else 0
    for i in range(nrev)
]

shared = [
    round(mat[i][j] * nj[i] * nj[j])
    for i in range(nrev)
    for j in range(i + 1, nrev)
]

sb = {
    "0": 0,
    "1": 0,
    "2": 0,
    "3–5": 0,
    "6+": 0
}

for x in shared:

    if x == 0:
        k = "0"
    elif x == 1:
        k = "1"
    elif x == 2:
        k = "2"
    elif x <= 5:
        k = "3–5"
    else:
        k = "6+"

    sb[k] += 1

SPAIRS = len(shared)


# ============================================================================
# C) COVERAGE
# ============================================================================

organs = [
    "leaf",
    "whole plant",
    "root",
    "stem",
    "flower",
    "reproductive",
    "litter",
    "life history"
]

cats = [
    "regulation",
    "material",
    "non-material"
]

cov = defaultdict(float)

for r in dm:
    cov[
        (
            r["trait_category"],
            r["ncp_category"]
        )
    ] += num(r["ntot"])


# ============================================================================
# D) EVIDENCE DEPTH
# ============================================================================

ps = defaultdict(float)

for r in dm:
    ps[r["ts"]] += num(r["ntot"])


db = {
    "1": 0,
    "2": 0,
    "3–5": 0,
    "6–10": 0,
    ">10": 0
}

for n in ps.values():

    if n <= 1:
        k = "1"
    elif n == 2:
        k = "2"
    elif n <= 5:
        k = "3–5"
    elif n <= 10:
        k = "6–10"
    else:
        k = ">10"

    db[k] += 1

NPAIR = len(ps)


# ============================================================================
# GRAPHICS SETTINGS
# ============================================================================

F = (
    "font-family='system-ui,-apple-system,"
    "\"Segoe UI\",sans-serif'"
)

BLUE = "#2a78d6"
MUTE = "#c8cbd1"
INK = "#0b0b0b"
SUB = "#000000" #changed from grey "#52514e" to be black
ORANGE = "#eb6834"

BACKGROUND = "#ffffff"

W = 1080
H = 620


# ============================================================================
# HELPERS
# ============================================================================

def esc(x):
    return (
        str(x)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def new_svg(title=None):

    s = [
        f'<svg width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" '
        f'xmlns="http://www.w3.org/2000/svg">',

        f'<rect width="{W}" height="{H}" '
        f'fill="{BACKGROUND}"/>'
    ]

    if title:
        s.append(
            f'<text x="50" y="45" {F} '
            f'font-size="21" font-weight="800" '
            f'fill="{INK}">{esc(title)}</text>'
        )

    return s


def close_svg(s):
    s.append("</svg>")


def axes(s, x, y, w, h):

    s.append(
        f'<line x1="{x}" y1="{y+h}" '
        f'x2="{x+w}" y2="{y+h}" '
        f'stroke="#9a9a94" stroke-width="1.2"/>'
    )

    s.append(
        f'<line x1="{x}" y1="{y}" '
        f'x2="{x}" y2="{y+h}" '
        f'stroke="#9a9a94" stroke-width="1.2"/>'
    )


def write_svg(s, filename):

    path = os.path.join(OUTDIR, filename)

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(s))

    print("wrote:", path)

    return path


def export_png(svg_path):
    png_path = os.path.splitext(svg_path)[0] + ".png"

    cmd = [
        INKSCAPE,
        svg_path,
        "--export-type=png",
        "--export-filename=" + png_path,
        "--export-dpi=300"
    ]

    print("Exporting:", os.path.basename(svg_path))

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    returncode = process.wait()

    if returncode != 0:
        raise RuntimeError(
            f"Inkscape failed for {svg_path} "
            f"(return code {returncode})"
        )

    print("wrote:", png_path)


# ============================================================================
# A — QUALITY VS YEAR
# ============================================================================

def make_quality():

    s = new_svg()

    x0 = 110
    y0 = 90
    PW = 850
    PH = 420

    s.append(
        f'<text x="{x0}" y="{y0-25}" {F} '
        f'font-size="34" font-weight="700" '
        f'fill="{INK}">'
        f'Booming, but not getting more rigorous'
        f'</text>'
    )

    axes(s, x0, y0, PW, PH)

    # ------------------------------------------------------------------------
    # axes
    # ------------------------------------------------------------------------

    ax = lambda v: x0 + PW * (v - 2017) / 7

    ay = lambda v: y0 + PH * (
        1 - v / 0.8
    )

    for g in (0, .2, .4, .6, .8):

        s.append(
            f'<line x1="{x0}" y1="{ay(g):.0f}" '
            f'x2="{x0+PW}" y2="{ay(g):.0f}" '
            f'stroke="#f0efe9"/>'
        )

        s.append(
            f'<text x="{x0-10}" y="{ay(g)+4:.0f}" '
            f'{F} font-size="24" fill="{SUB}" '
            f'text-anchor="end">{g:.1f}</text>'
        )

    for yr in range(2017, 2025):

        s.append(
            f'<text x="{ax(yr):.0f}" '
            f'y="{y0+PH+25}" {F} '
            f'font-size="24" fill="{SUB}" '
            f'text-anchor="middle">{yr}</text>'
        )

    # ------------------------------------------------------------------------
    # OLS
    # ------------------------------------------------------------------------

    yrs = [a for a, b in qy]
    vals = [b for a, b in qy]

    mx = st.mean(yrs)
    my = st.mean(vals)

    nn = len(qy)

    sxx = sum(
        (a - mx) ** 2
        for a in yrs
    )

    b1 = sum(
        (a - mx) * (b - my)
        for a, b in qy
    ) / sxx

    b0 = my - b1 * mx

    sse = sum(
        (
            b -
            (b0 + b1 * a)
        ) ** 2
        for a, b in qy
    )

    s2 = sse / (nn - 2)

    tc = 2.086

    top = []
    bot = []

    for gx in [
        2017 + .25 * k
        for k in range(29)
    ]:

        fit = b0 + b1 * gx

        sef = math.sqrt(
            s2 *
            (
                1 / nn +
                (gx - mx) ** 2 / sxx
            )
        )

        top.append(
            (
                gx,
                fit + tc * sef
            )
        )

        bot.append(
            (
                gx,
                fit - tc * sef
            )
        )

    poly = " ".join(
        f"{ax(gx):.1f},"
        f"{ay(min(.8, max(0, v))):.1f}"
        for gx, v in top
    )

    poly += " " + " ".join(
        f"{ax(gx):.1f},"
        f"{ay(min(.8, max(0, v))):.1f}"
        for gx, v in reversed(bot)
    )

    s.append(
        f'<polygon points="{poly}" '
        f'fill="{ORANGE}" '
        f'fill-opacity="0.10"/>'
    )

    s.append(
        f'<line '
        f'x1="{ax(2017):.0f}" '
        f'y1="{ay(b0+b1*2017):.0f}" '
        f'x2="{ax(2024):.0f}" '
        f'y2="{ay(b0+b1*2024):.0f}" '
        f'stroke="{ORANGE}" '
        f'stroke-width="4"/>'
    )

    # ------------------------------------------------------------------------
    # points
    # ------------------------------------------------------------------------

    for a, b in qy:

        s.append(
            f'<circle '
            f'cx="{ax(a):.0f}" '
            f'cy="{ay(b):.0f}" '
            f'r="8" fill="none" '
            f'stroke="{BLUE}" '
            f'stroke-width="1.7"/>'
        )

    # ------------------------------------------------------------------------
    # statistics
    # ------------------------------------------------------------------------

    ci_lo = (
        b1 -
        tc * math.sqrt(s2 / sxx)
    )

    ci_hi = (
        b1 +
        tc * math.sqrt(s2 / sxx)
    )

    s.append(
        f'<text x="{x0+PW}" '
        f'y="{y0+25}" {F} '
        f'font-size="26" '
        f'fill="{SUB}" '
        f'text-anchor="end">'
        f'β = {b1:.3f} yr⁻¹ '
        f'(95% CI {ci_lo:.3f}, {ci_hi:.3f}); '
        f'n = {nn}'
        f'</text>'
    )

    s.append(
        f'<text x="{x0-55}" '
        f'y="{y0+PH/2:.0f}" {F} '
        f'font-size="26" fill="{SUB}" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 '
        f'{x0-55} {y0+PH/2:.0f})">'
        f'review methodological quality score'
        f'</text>'
    )

    s.append(
        f'<text x="{x0+PW/2}" '
        f'y="{y0+PH+55}" {F} '
        f'font-size="26" fill="{SUB}" '
        f'text-anchor="middle">'
        f'publication year'
        f'</text>'
    )

    close_svg(s)

    return write_svg(
        s,
        "fig1a_quality_year.svg"
    )


# ============================================================================
# B — REDUNDANCY
# ============================================================================

def make_redundancy():

    s = new_svg()

    x0 = 110
    y0 = 90
    PW = 850
    PH = 420

    s.append(
        f'<text x="{x0}" y="{y0-25}" {F} '
        f'font-size="18" font-weight="700" '
        f'fill="{INK}">'
        f'b  Reviews are largely independent'
        f'</text>'
    )

    axes(s, x0, y0, PW, PH)

    order = [
        "0",
        "1",
        "2",
        "3–5",
        "6+"
    ]

    mxs = max(sb.values())

    bw = PW / len(order)

    for i, k in enumerate(order):

        v = sb[k]

        bh = (
            0.9 *
            PH *
            v /
            mxs
        )

        col = (
            BLUE
            if k == "0"
            else MUTE
        )

        x = (
            x0 +
            i * bw +
            bw * 0.16
        )

        y = (
            y0 +
            PH -
            bh
        )

        s.append(
            f'<rect '
            f'x="{x:.0f}" '
            f'y="{y:.0f}" '
            f'width="{bw*0.68:.0f}" '
            f'height="{bh:.0f}" '
            f'fill="{col}"/>'
        )

        s.append(
            f'<text '
            f'x="{x0+i*bw+bw/2:.0f}" '
            f'y="{y-10:.0f}" '
            f'{F} font-size="12" '
            f'font-weight="'
            f'{"700" if k=="0" else "400"}" '
            f'fill="{INK if k=="0" else SUB}" '
            f'text-anchor="middle">'
            f'{round(100*v/SPAIRS)}%'
            f'</text>'
        )

        s.append(
            f'<text '
            f'x="{x0+i*bw+bw/2:.0f}" '
            f'y="{y0+PH+28}" '
            f'{F} font-size="12" '
            f'fill="{SUB}" '
            f'text-anchor="middle">'
            f'{esc(k)}'
            f'</text>'
        )

    s.append(
        f'<text x="{x0+PW/2:.0f}" '
        f'y="{y0+PH+60}" {F} '
        f'font-size="12" fill="{SUB}" '
        f'text-anchor="middle">'
        f'primary studies shared between a pair of reviews '
        f'(n = {SPAIRS} pairs)'
        f'</text>'
    )

    s.append(
        f'<text x="{x0-55}" '
        f'y="{y0+PH/2:.0f}" {F} '
        f'font-size="12" fill="{SUB}" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 '
        f'{x0-55} {y0+PH/2:.0f})">'
        f'% of review pairs'
        f'</text>'
    )

    close_svg(s)

    return write_svg(
        s,
        "fig1b_redundancy.svg"
    )


# ============================================================================
# C — COVERAGE
# ============================================================================

def make_coverage():

    s = new_svg()

    x0 = 110
    y0 = 115
    PW = 850
    PH = 390

    s.append(
        f'<text x="{x0}" y="{y0-45}" {F} '
        f'font-size="18" font-weight="700" '
        f'fill="{INK}">'
        f'c  Evidence concentrated on leaf × regulating'
        f'</text>'
    )

    olab = {
        "leaf": "Leaf",
        "whole plant": "Whole-plant",
        "root": "Root",
        "stem": "Stem",
        "flower": "Flower",
        "reproductive": "Reprod.",
        "litter": "Litter",
        "life history": "Life hist."
    }

    clab = {
        "regulation": "Regulating",
        "material": "Material",
        "non-material": "Non-mat."
    }

    maxv = max(cov.values())

    label_width = 160

    cw = (
        PW -
        label_width
    ) / len(cats)

    ch = PH / len(organs)

    gx0 = x0 + label_width

    for j, cat in enumerate(cats):

        s.append(
            f'<text '
            f'x="{gx0+j*cw+cw/2:.0f}" '
            f'y="{y0-15}" {F} '
            f'font-size="12" fill="{SUB}" '
            f'text-anchor="middle">'
            f'{clab[cat]}'
            f'</text>'
        )

    for i, o in enumerate(organs):

        s.append(
            f'<text '
            f'x="{gx0-12}" '
            f'y="{y0+i*ch+ch/2+4:.0f}" '
            f'{F} font-size="12" '
            f'fill="{INK}" '
            f'text-anchor="end">'
            f'{olab[o]}'
            f'</text>'
        )

        for j, cat in enumerate(cats):

            v = cov.get(
                (o, cat),
                0
            )

            t = (
                math.log1p(v) /
                math.log1p(maxv)
                if maxv
                else 0
            )

            if v > 0:

                col = (
                    "#%02x%02x%02x" %
                    (
                        int(
                            255 -
                            (255-42) * t
                        ),
                        int(
                            255 -
                            (255-120) * t
                        ),
                        int(
                            255 -
                            (255-214) * t
                        )
                    )
                )

            else:

                col = "#f6f5f0"

            cx = (
                gx0 +
                j * cw
            )

            cy = (
                y0 +
                i * ch
            )

            s.append(
                f'<rect '
                f'x="{cx:.0f}" '
                f'y="{cy:.0f}" '
                f'width="{cw-3:.0f}" '
                f'height="{ch-3:.0f}" '
                f'fill="{col}" '
                f'stroke="#ffffff" '
                f'stroke-width="1"/>'
            )

            if v > 0:

                txtcol = (
                    "#ffffff"
                    if t > 0.5
                    else SUB
                )

                s.append(
                    f'<text '
                    f'x="{cx+cw/2:.0f}" '
                    f'y="{cy+ch/2+4:.0f}" '
                    f'{F} font-size="11" '
                    f'fill="{txtcol}" '
                    f'text-anchor="middle">'
                    f'{int(v)}'
                    f'</text>'
                )

    s.append(
        f'<text '
        f'x="{gx0+(PW-label_width)/2:.0f}" '
        f'y="{y0+PH+35}" {F} '
        f'font-size="12" fill="{SUB}" '
        f'text-anchor="middle">'
        f'number of primary studies'
        f'</text>'
    )

    close_svg(s)

    return write_svg(
        s,
        "fig1c_coverage.svg"
    )


# ============================================================================
# D — DEPTH
# ============================================================================

def make_depth():

    s = new_svg()

    x0 = 110
    y0 = 90
    PW = 850
    PH = 420

    s.append(
        f'<text x="{x0}" y="{y0-25}" {F} '
        f'font-size="34" font-weight="700" '
        f'fill="{INK}">'
        f'Most relationships rest on a single study'
        f'</text>'
    )

    axes(s, x0, y0, PW, PH)

    order = [
        "1",
        "2",
        "3–5",
        "6–10",
        ">10"
    ]

    mx2 = max(db.values())

    bw = PW / len(order)

    for i, k in enumerate(order):

        v = db[k]

        bh = (
            0.9 *
            PH *
            v /
            mx2
        )

        col = (
            BLUE
            if k == "1"
            else MUTE
        )

        x = (
            x0 +
            i * bw +
            bw * 0.16
        )

        y = (
            y0 +
            PH -
            bh
        )

        s.append(
            f'<rect '
            f'x="{x:.0f}" '
            f'y="{y:.0f}" '
            f'width="{bw*0.68:.0f}" '
            f'height="{bh:.0f}" '
            f'fill="{col}"/>'
        )

        s.append(
            f'<text '
            f'x="{x0+i*bw+bw/2:.0f}" '
            f'y="{y-10:.0f}" '
            f'{F} font-size="26" '
            f'font-weight="'
            f'{"700" if k=="1" else "400"}" '
            f'fill="{INK if k=="1" else SUB}" '
            f'text-anchor="middle">'
            f'{round(100*v/NPAIR)}%'
            f'</text>'
        )

        s.append(
            f'<text '
            f'x="{x0+i*bw+bw/2:.0f}" '
            f'y="{y0+PH+28}" '
            f'{F} font-size="26" '
            f'fill="{SUB}" '
            f'text-anchor="middle">'
            f'{esc(k)}'
            f'</text>'
        )

    s.append(
        f'<text x="{x0+PW/2:.0f}" '
        f'y="{y0+PH+60}" {F} '
        f'font-size="26" fill="{SUB}" '
        f'text-anchor="middle">'
        f'primary studies per trait–NCP relationship '
        f'(n = {NPAIR})'
        f'</text>'
    )

    s.append(
        f'<text x="{x0-55}" '
        f'y="{y0+PH/2:.0f}" {F} '
        f'font-size="26" fill="{SUB}" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 '
        f'{x0-55} {y0+PH/2:.0f})">'
        f'% of relationships'
        f'</text>'
    )

    close_svg(s)

    return write_svg(
        s,
        "fig1d_depth.svg"
    )


# ============================================================================
# RUN
# ============================================================================
if __name__ == "__main__":
    print()
    print("=" * 70)
    print("Generating individual figures")
    print("=" * 70)
    print()

    files = [
        make_quality(),
        make_redundancy(),
        make_coverage(),
        make_depth()
    ]

    print()
    print("=" * 70)
    print("Exporting PNG files with Inkscape")
    print("=" * 70)
    print()

    for svg in files:
        export_png(svg)

    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)
    
    
