# ============================================================================
#  fig_main_tradeoff.py  —  Main figure: "Same trait, opposite contributions"
#
#  Trait x NCP posterior effect-direction matrix from the PAIR-LEVEL (simple)
#  model (no cross-NCP borrowing) — the correct parameterisation for reading
#  per-pair directions. Colour = continuous P(effect > 0) on a diverging RdBu
#  ramp; cell number = primary studies. Two conventional reading aids overlaid
#  on the continuous scale:
#     bold black outline  = strong      (P > 0.90 or < 0.10)
#     thin grey outline   = suggestive  (P > 0.75 or < 0.25)
#  Cells with < 3 primary studies are faded (under-evidenced).
#
#  INPUT   outputs/reversal_check.csv   (from R/check_reversals.R)
#          data/data_model.csv          (trait/NCP hierarchy)
#  OUTPUT  outputs/fig_main_tradeoff.svg / .png        (13 curated traits)
#          outputs/fig_supp_tradeoff.svg / .png        (19 traits, wider set)
#
#  RUN     python3 figures/fig_main_tradeoff.py
#  DEPS    stdlib only for SVG. PNG is written too IF cairosvg is installed
#          (pip install cairosvg); otherwise open the SVG in any browser.
# ============================================================================
import csv, os, sys

# ---- Resolve pipeline paths (script lives in <root>/figures/) --------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE) if os.path.basename(HERE) == "figures" else HERE
OUT_DIR  = os.path.join(ROOT, "outputs")
DATA_DIR = os.path.join(ROOT, "data")
REV = os.path.join(OUT_DIR, "reversal_check.csv")
DM  = os.path.join(DATA_DIR, "data_model.csv")
if not os.path.exists(REV):
    sys.exit("Missing %s — run R/check_reversals.R first." % REV)

# ---- Load posterior direction (pair-level) + study counts ------------------
P, NS = {}, {}
for r in csv.DictReader(open(REV)):
    key = (r["trait_grouped"], r["ncp"])
    P[key]  = float(r["P_simple"])
    NS[key] = float(r["n_studies"])

# ---- Trait rows: short label -> trait_grouped, in super-groups --------------
GROUPS_MAIN = [
 ("Leaf economics & size", [("SLA","specific leaf area"),("LDMC","leaf dry matter content"),
    ("Leaf N","leaf nitrogen content"),("Leaf P","leaf phosphorus content"),
    ("Leaf area","leaf area"),("LAI","leaf area index")]),
 ("Whole-plant", [("Height","height"),("Crown size","crown size"),
    ("Wood density","wood density"),("Foliage density","foliage density")]),
 ("Root", [("Root diameter","root diameter"),("Root length","root length"),
    ("Rooting depth","rooting depth")]),
]
GROUPS_SUPP = [
 ("Leaf economics & size", [("SLA","specific leaf area"),("LDMC","leaf dry matter content"),
    ("Leaf N","leaf nitrogen content"),("Leaf P","leaf phosphorus content"),("Leaf C","leaf carbon content"),
    ("Leaf area","leaf area"),("LAI","leaf area index"),("Leaf thickness","leaf thickness"),
    ("Leaf toughness","leaf toughness")]),
 ("Whole-plant", [("Height","height"),("Crown size","crown size"),
    ("Wood density","wood density"),("Foliage density","foliage density"),
    ("Aboveg. biomass","aboveground biomass"),("Aboveg. litter","aboveground litter")]),
 ("Root", [("Root diameter","root diameter"),("Root length","root length"),("Rooting depth","rooting depth")]),
 ("Reprod.", [("Seed mass","seed mass")]),
]

# ---- Columns: short -> full NCP name, grouped by NCP category ----------------
COLS = [
 ("REGULATION", [("Soil formation","formation and protection of soils"),("Climate","regulation of climate"),
    ("Freshwater qty","regulation of freshwater quantity"),("Freshwater qual","regulation of freshwater quality"),
    ("Hazards","regulation of hazards and extreme events"),("Pest reg.","regulation of detrimental organisms"),
    ("Air quality","regulation of air quality"),("Pollination","pollination and seed dispersal"),
    ("Habitat","habitat creation and maintenance")]),
 ("MATERIAL", [("Food & feed","food and feed"),("Materials","materials and assistance"),
    ("Medicinal","medicinal and genetic ressources")]),
 ("NON-MAT.", [("Experiences","experiences")]),
]

# ---- Continuous diverging colour ramp (ColorBrewer RdBu: neg=red, pos=blue) --
ANCHORS = [(0.0,(178,24,43)),(0.1,(214,96,77)),(0.25,(244,165,130)),(0.4,(253,219,199)),
           (0.5,(238,236,230)),(0.6,(209,229,240)),(0.75,(146,197,222)),(0.9,(67,147,195)),(1.0,(33,102,172))]
def _lerp(a,b,t): return tuple(round(a[i]+(b[i]-a[i])*t) for i in range(3))
def colour(p):
    for i in range(len(ANCHORS)-1):
        x0,c0 = ANCHORS[i]; x1,c1 = ANCHORS[i+1]
        if p <= x1:
            t = (p-x0)/(x1-x0) if x1 > x0 else 0
            r,g,b = _lerp(c0,c1,t); return "#%02x%02x%02x" % (r,g,b)
    return "#2166ac"

def esc(x): return str(x).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")

CW, CH, X0, Y0, PAD, FADE = 44, 30, 210, 190, 40, "#f7f6f3"
FONT = "font-family='system-ui,-apple-system,\"Segoe UI\",sans-serif'"
GRAY = ["#f4f3ef","#faf9f6","#f4f3ef","#faf9f6"]

def build(groups, out_svg):
    flatcols = [c for _,cs in COLS for c in cs]
    nrows = sum(len(cs) for _,cs in groups); ncols = len(flatcols)
    Wd = X0 + ncols*CW + PAD; Ht = Y0 + nrows*CH + 150
    s = [f'<svg width="{Wd}" height="{Ht}" viewBox="0 0 {Wd} {Ht}" xmlns="http://www.w3.org/2000/svg">',
         f'<rect width="{Wd}" height="{Ht}" fill="#fcfcfb"/>',
         f'<text x="40" y="40" {FONT} font-size="18" font-weight="700" fill="#0b0b0b">Same trait, opposite contributions</text>',
         f'<text x="40" y="62" {FONT} font-size="11" fill="#52514e">Posterior effect direction of each well-studied trait–NCP relationship (pair-level model). Cell number = primary studies.</text>']
    # column category headers
    cx = X0
    for cat,cs in COLS:
        x1 = cx; x2 = cx + len(cs)*CW
        s.append(f'<line x1="{x1+2}" y1="90" x2="{x2-2}" y2="90" stroke="#52514e" stroke-width="3"/>')
        s.append(f'<text x="{(x1+x2)/2}" y="84" {FONT} font-size="11" font-weight="700" fill="#52514e" text-anchor="middle">{esc(cat)}</text>')
        cx = x2
    for j,(short,_) in enumerate(flatcols):
        tx = X0 + j*CW + CW/2
        s.append(f'<text x="{tx}" y="182" {FONT} font-size="10.5" fill="#52514e" text-anchor="start" transform="rotate(-45 {tx} 182)">{esc(short)}</text>')
    # rows
    ri = 0
    for gi,(gname,cs) in enumerate(groups):
        gy = Y0 + ri*CH; gh = len(cs)*CH
        s.append(f'<rect x="{X0}" y="{gy}" width="{ncols*CW}" height="{gh}" fill="{GRAY[gi%len(GRAY)]}"/>')
        s.append(f'<line x1="48" y1="{gy+4}" x2="48" y2="{gy+gh-4}" stroke="#52514e" stroke-width="2.5"/>')
        if gh >= CH*1.5:
            s.append(f'<text x="40" y="{gy+gh/2}" {FONT} font-size="10.5" font-weight="700" fill="#52514e" text-anchor="middle" transform="rotate(-90 40 {gy+gh/2})">{esc(gname)}</text>')
        for short,tg in cs:
            rY = Y0 + ri*CH
            s.append(f'<text x="202" y="{rY+CH/2+3.5}" {FONT} font-size="11.5" font-weight="600" fill="#0b0b0b" text-anchor="end">{esc(short)}</text>')
            for j,(_,ncp) in enumerate(flatcols):
                cxp = X0 + j*CW; p = P.get((tg,ncp)); n = NS.get((tg,ncp),0)
                if p is None or n < 3:
                    s.append(f'<rect x="{cxp+1}" y="{rY+1}" width="{CW-2}" height="{CH-2}" fill="{FADE}" rx="2"/>')
                else:
                    fc = colour(p)
                    if p > 0.90 or p < 0.10:   stroke = ' stroke="#1a1a1a" stroke-width="1.8"'
                    elif p > 0.75 or p < 0.25: stroke = ' stroke="#6f6f6f" stroke-width="0.9"'
                    else:                      stroke = ''
                    s.append(f'<rect x="{cxp+1}" y="{rY+1}" width="{CW-2}" height="{CH-2}" fill="{fc}" rx="2"{stroke}/>')
                    tcol = "#ffffff" if (p < 0.16 or p > 0.84) else "#52514e"
                    s.append(f'<text x="{cxp+CW/2}" y="{rY+CH/2+3.5}" {FONT} font-size="8.5" fill="{tcol}" text-anchor="middle">{int(n)}</text>')
            ri += 1
    # separators
    yend = Y0 + nrows*CH; acc = 0
    for gname,cs in groups[:-1]:
        acc += len(cs); yy = Y0 + acc*CH
        s.append(f'<line x1="{X0}" y1="{yy}" x2="{X0+ncols*CW}" y2="{yy}" stroke="#c3c2b7" stroke-width="1.5"/>')
    cxx = X0
    for cat,cs in COLS[:-1]:
        cxx += len(cs)*CW
        s.append(f'<line x1="{cxx}" y1="92" x2="{cxx}" y2="{yend}" stroke="#c3c2b7" stroke-width="1.5"/>')
    # legend ramp + keys
    ly = yend + 30; n = 21
    s.append(f'<text x="{X0}" y="{ly}" {FONT} font-size="11" fill="#52514e">Effect direction (posterior probability)</text>')
    for i in range(n):
        s.append(f'<rect x="{X0+i*17}" y="{ly+8}" width="17" height="16" fill="{colour(i/(n-1))}"/>')
    s.append(f'<text x="{X0}" y="{ly+40}" {FONT} font-size="10" fill="#b2182b">negative</text>')
    s.append(f'<text x="{X0+n*17/2}" y="{ly+40}" {FONT} font-size="10" fill="#898781" text-anchor="middle">uncertain</text>')
    s.append(f'<text x="{X0+n*17}" y="{ly+40}" {FONT} font-size="10" fill="#2166ac" text-anchor="end">positive</text>')
    kx = X0 + n*17 + 40
    s.append(f'<rect x="{kx}" y="{ly+6}" width="15" height="14" fill="#c45863" stroke="#1a1a1a" stroke-width="1.8"/>')
    s.append(f'<text x="{kx+21}" y="{ly+17}" {FONT} font-size="10" fill="#52514e">strong (P&#62;0.90)</text>')
    s.append(f'<rect x="{kx}" y="{ly+24}" width="15" height="14" fill="#e8b9bc" stroke="#6f6f6f" stroke-width="0.9"/>')
    s.append(f'<text x="{kx+21}" y="{ly+35}" {FONT} font-size="10" fill="#52514e">suggestive (P&#62;0.75)</text>')
    s.append(f'<rect x="{kx}" y="{ly+42}" width="15" height="14" fill="{FADE}" stroke="#e1e0d9"/>')
    s.append(f'<text x="{kx+21}" y="{ly+53}" {FONT} font-size="10" fill="#898781">&#60; 3 studies</text>')
    s.append('</svg>')
    open(out_svg,"w").write("\n".join(s))
    # optional PNG
    try:
        import cairosvg
        cairosvg.svg2png(url=out_svg, write_to=out_svg.replace(".svg",".png"), output_width=1400)
        png = " (+ PNG)"
    except Exception:
        png = " (SVG only; pip install cairosvg for PNG)"
    print("wrote %s  [%d traits x %d NCPs]%s" % (out_svg, nrows, ncols, png))

if __name__ == "__main__":
    build(GROUPS_MAIN, os.path.join(OUT_DIR, "fig_main_tradeoff.svg"))
    build(GROUPS_SUPP, os.path.join(OUT_DIR, "fig_supp_tradeoff.svg"))
