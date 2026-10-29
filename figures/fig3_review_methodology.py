# ============================================================================
#  fig3_review_methodology.py  —  Section-3 figure (Python / stdlib SVG)
#     "How a review is done predicts what it concludes"
#
#  Each dot is one review, placed by its positivity (share of directional
#  evidence reported as positive), split by three synthesis-rigour attributes.
#  Group means (diamonds) + the shift between them, annotated with the p-value
#  from a study-level binomial GLMM (random review & trait-NCP-pair effects) —
#  computed by figures/fig3_review_methodology.R, which writes outputs/fig3_glmm.csv.
#  If that file is absent, the figure still renders (p shown as "—").
#
#  INPUT   data/data_model.csv, data/meta-data2.csv, [outputs/fig3_glmm.csv]
#  OUTPUT  outputs/fig3_review_methodology.svg (+ .png if cairosvg installed)
#  RUN     python3 figures/fig3_review_methodology.py
# ============================================================================
import csv, os, sys, statistics as st
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))

HERE = os.getcwd()
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == "figures" else HERE
DM   = os.path.join(ROOT, "data", "data_model.csv")
META = os.path.join(ROOT, "data", "meta-data2.csv")
GLMM = os.path.join(ROOT, "outputs", "fig3_glmm.csv")
OUT  = os.path.join(ROOT, "outputs", "fig3_review_methodology.svg")
for p in (DM, META):
    if not os.path.exists(p): sys.exit("Missing %s" % p)

# ---- per-review positivity + rigour attributes ------------------------------
dm = list(csv.DictReader(open(DM)))
agg = defaultdict(lambda: [0.0, 0.0])
for r in dm:
    for k in ("Npos", "Nneg"):
        try: v = float(r[k])
        except Exception: v = 0.0
        agg[r["source_id"]][0 if k == "Npos" else 1] += v
pos = {k: p / (p + n) for k, (p, n) in agg.items() if (p + n) > 0}
meta = {}
for r in csv.DictReader(open(META, encoding="latin1"), delimiter=";"):
    sid = r["source_id"].strip()
    if sid not in meta: meta[sid] = r          # first row per review (match R)
def at(sid, c): return (meta.get(sid, {}).get(c, "") or "").strip().lower()
rows = [{"pos": pos[sid],
         "trait_def": at(sid, "Specific_Traits_Definitions"),
         "method": "vote" if "vote" in at(sid, "Statistical_Method_to_Assess_the_Effect") else "narrative",
         "eco": at(sid, "Ecological_standpoint")} for sid in pos]

# ---- p-values from the GLMM (if computed) -----------------------------------
pval = {"method": "—", "trait_def": "—", "eco": "—"}
if os.path.exists(GLMM):
    for r in csv.DictReader(open(GLMM)):
        try: p = float(r["p"])
        except Exception: continue
        pval[r["attribute"]] = ("p < 0.001" if p < 0.001 else "p = %.2f" % p) if p < 0.10 else "n.s."

def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
F = "font-family='system-ui,-apple-system,\"Segoe UI\",sans-serif'"
BLUE, ORANGE, INK, SUB = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
W = 1040
XL, XR, X0v, X1v = 470, 940, 0.50, 1.00
def X(v): return XL + (XR - XL) * (v - X0v) / (X1v - X0v)
R, rowH, blockH, blockGap, top = 5.2, 40, 80, 54, 210

def bee(vals):
    placed = []
    for v in sorted(vals):
        x = X(v); step = R * 2 + 1.2; k = 0; yoff = 0
        while True:
            for cand in ([0] if k == 0 else [k * step, -k * step]):
                if all(abs(x - px) > R * 1.7 or abs(cand - py) > R * 1.7 for px, py in placed):
                    yoff = cand; break
            else:
                k += 1; continue
            break
        placed.append((x, yoff))
    return placed

blocks = [
 ("a Synthesis method",     ("Quantitative", lambda r: r["method"] == "vote"),     ("Narrative only", lambda r: r["method"] == "narrative"),  pval["method"]),
 ("b Trait definitions",    ("Explicit",     lambda r: r["trait_def"] == "yes"),   ("None",           lambda r: r["trait_def"] == "no"),        pval["trait_def"]),
 ("c Ecological framework", ("Explicit",     lambda r: r["eco"] == "yes"),         ("None",           lambda r: r["eco"] == "no"),              pval["eco"]),
]
H = top + len(blocks) * (blockH + blockGap) + 96
def mean_of(sel): v = [r["pos"] for r in rows if sel(r)]; return st.mean(v), len(v), v

s = [f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg">',
     f'<rect width="{W}" height="{H}" fill="#fcfcfb"/>',
     f'<text x="44" y="52" {F} font-size="23" font-weight="800" fill="{INK}">How a review is done predicts what it concludes</text>',
     f'<text x="44" y="80" {F} font-size="13.5" fill="{SUB}">Each dot is one systematic review, placed by the share of its directional evidence reported as positive.</text>',
     f'<text x="44" y="99" {F} font-size="13.5" fill="{SUB}">Reviews using less rigorous methods report systematically more positive trait–NCP relationships. Exploratory (n = {len(rows)}).</text>']
#s.append(f'<text x="{X(0.55):.0f}" y="{top-52}" {F} font-size="11" font-weight="700" fill="{BLUE}" text-anchor="middle">← more rigorous, less positive</text>')
#s.append(f'<text x="{X(0.93):.0f}" y="{top-52}" {F} font-size="11" font-weight="700" fill="{ORANGE}" text-anchor="middle">less rigorous, more positive →</text>')
for g in range(50, 101, 10):
    gx = X(g / 100); s.append(f'<line x1="{gx:.0f}" y1="{top-24}" x2="{gx:.0f}" y2="{H-70}" stroke="#f1f0eb"/>')
    s.append(f'<text x="{gx:.0f}" y="{top-30}" {F} font-size="9.5" fill="{SUB}" text-anchor="middle">{g}%</text>')

y = top
for title, (l1, s1), (l2, s2), pv in blocks:
    m1, n1, v1 = mean_of(s1); m2, n2, v2 = mean_of(s2)
    yr = y + rowH / 2; yl = y + rowH + rowH / 2
    s.append(f'<text x="44" y="{y+rowH+2:.0f}" {F} font-size="14" font-weight="700" fill="{INK}">{esc(title)}</text>')
    x1, x2 = X(m1), X(m2); midy = (yr + yl) / 2
    s.append(f'<rect x="{min(x1,x2):.1f}" y="{yr:.1f}" width="{abs(x2-x1):.1f}" height="{yl-yr:.1f}" fill="{ORANGE}" opacity="0.06"/>')
    s.append(f'<line x1="{x1:.1f}" y1="{midy:.1f}" x2="{x2-6:.1f}" y2="{midy:.1f}" stroke="{INK}" stroke-width="1.5"/>')
    s.append(f'<path d="M{x2:.1f},{midy:.1f} L{x2-9:.1f},{midy-4:.1f} L{x2-9:.1f},{midy+4:.1f} Z" fill="{INK}"/>')
    s.append(f'<text x="{(x1+x2)/2:.1f}" y="{midy-8:.1f}" {F} font-size="11" font-weight="700" fill="{INK}" text-anchor="middle">+{round((m2-m1)*100)} pts</text>')
    s.append(f'<text x="{(x1+x2)/2:.1f}" y="{midy+16:.1f}" {F} font-size="9.5" font-weight="700" fill="{"#0b0b0b" if pv not in ("n.s.","—") else "#8a8f98"}" text-anchor="middle">{esc(pv)}</text>')
    for lab, vals, col, ry, mv, nn in [(l1, v1, BLUE, yr, m1, n1), (l2, v2, ORANGE, yl, m2, n2)]:
        s.append(f'<text x="{XL-16}" y="{ry-1:.0f}" {F} font-size="11.5" fill="{INK}" text-anchor="end">{esc(lab)}</text>')
        s.append(f'<text x="{XL-16}" y="{ry+13:.0f}" {F} font-size="8.5" fill="{SUB}" text-anchor="end">n={nn}</text>')
        for dx, dy in bee(vals):
            s.append(f'<circle cx="{dx:.1f}" cy="{ry+dy:.1f}" r="{R}" fill="{col}" fill-opacity="0.5" stroke="{col}" stroke-width="0.8"/>')
        mx = X(mv)
        s.append(f'<path d="M{mx:.1f},{ry-8:.1f} L{mx+7:.1f},{ry:.1f} L{mx:.1f},{ry+8:.1f} L{mx-7:.1f},{ry:.1f} Z" fill="{col}" stroke="#fff" stroke-width="1.4"/>')
        ly = ry - 13 if col == BLUE else ry + 20
        s.append(f'<text x="{mx:.1f}" y="{ly:.0f}" {F} font-size="10.5" font-weight="700" fill="{col}" text-anchor="middle">{round(mv*100)}%</text>')
    y += blockH + blockGap

s.append(f'<text x="{(XL+XR)/2:.0f}" y="{H-52}" {F} font-size="10.5" fill="{SUB}" text-anchor="middle">Share of a review’s directional evidence reported as positive (mean ◆ per group)</text>')
s.append(f'<text x="44" y="{H-32}" {F} font-size="12" font-weight="600" fill="{INK}">The between-review differences track how evidence is synthesised — not the ecosystem studied (biome ≈ 10%).</text>')
s.append(f'<text x="44" y="{H-14}" {F} font-size="10.5" fill="{SUB}">p from a study-level mixed-effects model (2,203 outcomes; random review &amp; trait–NCP-pair effects). Recovery + data-rich-subset checks: Supplementary.</text>')
s.append('</svg>')
open(OUT, "w").write("\n".join(s))
import subprocess


png_file = OUT.replace(".svg", ".png")

# Essayer d'abord avec CairoSVG
try:
    import cairosvg
    cairosvg.svg2png(
        url=OUT,
        write_to=png_file,
        output_width=1400
    )
    print("PNG créé :", png_file)

except ImportError:
    # Sinon, essayer Inkscape
    inkscape = r"C:\Program Files\Inkscape\bin\inkscape.exe"

    if os.path.exists(inkscape):
        subprocess.run(
            [inkscape, OUT, f"--export-filename={png_file}"],
            check=True
        )
        print("PNG créé :", png_file)
    else:
        print("CairoSVG absent et Inkscape introuvable.")
        

