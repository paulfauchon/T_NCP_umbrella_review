# ============================================================================
#  fig2_review_disagreement.py  —  Section-2 figure (Python / stdlib SVG)
#     "The same relationship, opposite verdicts across reviews"
#
#  The human-readable face of the between-review result: for well-covered
#  trait-NCP relationships, each dot is ONE systematic review's overall verdict
#  (negative / no clear direction / positive). Reviews contradict one another
#  pair by pair — the reported direction depends more on which review reports it
#  than on the trait or the service. (The full variance decomposition + the
#  robustness battery are reported in the text / Supplementary.)
#
#  INPUT   data/data_model.csv
#  OUTPUT  outputs/fig2_review_disagreement.svg (+ .png if cairosvg installed)
#  RUN     python3 figures/fig2_review_disagreement.py
#  PARAMS  MIN_REVIEWS (default 5), N_PAIRS (default 10)
# ============================================================================
import csv, os, sys
from collections import defaultdict

MIN_REVIEWS = 5
N_PAIRS     = 10

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == "figures" else HERE
DM   = os.path.join(ROOT, "data", "data_model.csv")
OUT  = os.path.join(ROOT, "outputs", "fig2_review_disagreement.svg")
if not os.path.exists(DM):
    sys.exit("Missing %s" % DM)

dm = list(csv.DictReader(open(DM)))
for r in dm:
    for k in ("Npos", "Nneg", "Nnull"):
        try: r[k] = float(r[k])
        except Exception: r[k] = 0.0

def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

# ---- one direction per UNIQUE review per pair -------------------------------
pair = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0, 0.0]))   # ts -> review -> [pos,neg,null]
for r in dm:
    a = pair[r["ts"]][r["source_id"]]
    a[0] += r["Npos"]; a[1] += r["Nneg"]; a[2] += r["Nnull"]

def summarise(ts):
    npos = nneg = nnull = 0
    for _, (p, n, z) in pair[ts].items():
        if   p > n: npos  += 1
        elif n > p: nneg  += 1
        else:       nnull += 1
    return npos, nneg, nnull

cand = []
for ts in pair:
    npos, nneg, nnull = summarise(ts)
    nrev = npos + nneg + nnull
    if nrev >= MIN_REVIEWS:
        cand.append((ts, nrev, npos, nneg, nnull, min(npos, nneg)))   # last = disagreement
cand.sort(key=lambda x: (-x[5], -x[1]))     # most-contested, then best-covered
rows = cand[:N_PAIRS]

SHORT_TRAIT = {"specific leaf area":"Specific leaf area","leaf nitrogen content":"Leaf nitrogen",
 "leaf dry matter content":"Leaf dry-matter content","height":"Plant height","wood density":"Wood density",
 "root diameter":"Root diameter","aboveground biomass":"Aboveground biomass","leaf area":"Leaf area",
 "leaf phosphorus content":"Leaf phosphorus","crown size":"Crown size","leaf carbon content":"Leaf carbon"}
SHORT_NCP = {"regulation of climate":"climate regulation","formation and protection of soils":"soil formation",
 "regulation of freshwater quantity":"freshwater regulation","food and feed":"food & feed",
 "regulation of hazards and extreme events":"hazard regulation","regulation of detrimental organisms":"pest regulation"}
def nicelabel(ts):
    t, n = ts.split("__")
    return "%s  →  %s" % (SHORT_TRAIT.get(t, t.capitalize()), SHORT_NCP.get(n, n))

# ---- SVG --------------------------------------------------------------------
F   = "font-family='system-ui,-apple-system,\"Segoe UI\",sans-serif'"
RED, BLUE, GRAY = "#d1495b", "#2a78d6", "#b9bcc2"
INK, SUB = "#0b0b0b", "#52514e"
X0, COLW, R, rowH = 380, 150, 6.5, 54
xn, xz, xp = X0 + COLW*0.5, X0 + COLW*1.5, X0 + COLW*2.5
W = X0 + COLW*3 + 50
H = 150 + len(rows)*rowH + 96
top = 112

s = [f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">',
     f'<rect width="{W}" height="{H}" fill="#fcfcfb"/>',
     f'<text x="44" y="46" {F} font-size="20" font-weight="700" fill="{INK}">The same relationship, opposite verdicts across reviews</text>',
     f'<text x="44" y="70" {F} font-size="12.5" fill="{SUB}">Each dot is one systematic review’s overall conclusion for that trait–NCP relationship (well-covered pairs, ≥{MIN_REVIEWS} reviews).</text>']
for cx, lab, col in [(xn,"reported NEGATIVE",RED),(xz,"no clear direction",GRAY),(xp,"reported POSITIVE",BLUE)]:
    s.append(f'<text x="{cx}" y="{top-8}" {F} font-size="11.5" font-weight="700" fill="{col}" text-anchor="middle" letter-spacing="0.3">{lab}</text>')
for bx in (X0+COLW, X0+COLW*2):
    s.append(f'<line x1="{bx}" y1="{top}" x2="{bx}" y2="{top+len(rows)*rowH}" stroke="#ecebe6" stroke-width="1"/>')

def dots(cx, cy, count, col):
    if count == 0: return
    gap = R*2 + 3; per = 6
    nlines = (count + per - 1) // per
    for i in range(count):
        li, idx = divmod(i, per)
        inrow = min(count - li*per, per)
        tw = (inrow - 1) * gap
        x = cx - tw/2 + idx*gap
        y = cy + (li - (nlines-1)/2) * (R*2 + 3)
        s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{R}" fill="{col}" stroke="#fff" stroke-width="1.1"/>')

for i, (ts, nrev, npos, nneg, nnull, _) in enumerate(rows):
    cy = top + i*rowH + rowH/2
    if i > 0:
        s.append(f'<line x1="44" y1="{top+i*rowH}" x2="{X0+COLW*3}" y2="{top+i*rowH}" stroke="#f0efe9" stroke-width="1"/>')
    s.append(f'<text x="{X0-18}" y="{cy-3}" {F} font-size="12.5" font-weight="600" fill="{INK}" text-anchor="end">{esc(nicelabel(ts))}</text>')
    s.append(f'<text x="{X0-18}" y="{cy+14}" {F} font-size="9.5" fill="{SUB}" text-anchor="end">{nrev} reviews</text>')
    dots(xn, cy, nneg, RED); dots(xz, cy, nnull, GRAY); dots(xp, cy, npos, BLUE)

fy = top + len(rows)*rowH + 38
s.append(f'<text x="44" y="{fy}" {F} font-size="12.5" font-weight="600" fill="{INK}">Reviews contradict one another pair by pair — the reported direction depends more on</text>')
s.append(f'<text x="44" y="{fy+19}" {F} font-size="12.5" font-weight="600" fill="{INK}">which review reports it than on the trait or the service involved.</text>')
s.append('</svg>')
open(OUT, "w").write("\n".join(s))
try:
    import cairosvg
    cairosvg.svg2png(url=OUT, write_to=OUT.replace(".svg", ".png"), output_width=1300)
    png = " (+ PNG)"
except Exception:
    png = " (SVG only; pip install cairosvg for PNG)"
print("wrote %s  [%d pairs]%s" % (OUT, len(rows), png))
