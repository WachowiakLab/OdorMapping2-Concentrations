"""Build the shareable figure site from the notebook source.

One self-contained HTML file: a contents page plus one page per figure, with
the human caption, a collapsible "The grisly details" carrying the LLM-register
caption, and the section's code — analysis code open, composition code folded.

    python dev/build_site.py                 # -> notebooks/site/highconc_figures.html
    python dev/build_site.py --out PATH

Everything comes from `notebooks/figures.py` (both caption registers, tagged)
and `notebooks/outputs/*.png` (the rendered figures, embedded as data URIs).
Cell stdout is lifted from the executed `figures_LLM.ipynb`, matched to source
cells by position, so the printed numbers appear under the code that produced
them without this script running anything.

Nothing here recomputes a figure. Run the notebook first; if a stem is missing
from `notebooks/outputs/` the build says which and stops, rather than shipping
a page with a hole in it.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO / "notebooks" / "figures.py"
EXECUTED = REPO / "notebooks" / "figures_LLM.ipynb"
OUTPUTS = REPO / "notebooks" / "outputs"
#: Netlify (and every other static host) serves `index.html` at the site root,
#: so the build writes that name. `notebooks/site/` is the publish directory:
#: one file, no assets, nothing else to upload.
DEFAULT_OUT = REPO / "notebooks" / "site" / "index.html"

# Figures in page order: stem -> (short label, figure name, one-line summary).
FIGURES = {
    "fig3_correlational_structure": (
        "Figure 3", "Odorant–odorant response similarity",
        "How similar the glomerular response patterns of odorant pairs are, at "
        "both concentrations, against chemical similarity."),
    "fig3_supp_threshold_robustness": (
        "Figure 3 supp.", "The recruitment matrix at five cutoffs",
        "Whether the recruitment result depends on where the assignment "
        "threshold is drawn."),
    "fig4_tuning_by_cluster": (
        "Figure 4", "The seven glomerular groups",
        "The groups seen four ways: as a correlation matrix, as raw responses, "
        "one mouse at a time, and summarised by chemical class."),
    "fig5_spatial_by_cluster": (
        "Figure 5", "The groups on the dorsal bulb",
        "Where each of the seven groups sits, all four mice overlaid."),
    "fig6_pls": (
        "Figure 6", "A linear model of the low → high change",
        "One model fitted across all odorants, predicting held-out "
        "high-concentration responses."),
    "fig7_motif_transfer": (
        "Figure 7", "Directional transfer of tuning",
        "Whether adding one group's low-concentration profile improves "
        "prediction of another's high-concentration response."),
    "supp_nmf_dim_choice": (
        "Figure S1", "How the number of groups was chosen",
        "Variance explained against shuffled data, and stability across 50 "
        "random starts."),
    "supp_per_subject": (
        "Figure S2", "The same structure within each mouse",
        "The group structure computed within each animal separately, rather "
        "than on the pooled data."),
    "supp_tuning_flow": (
        "Figure S3", "Direction of change through tuning space",
        "Whether glomeruli in the same group move in similar directions as "
        "concentration rises."),
    "supp_parity": (
        "Figure S4", "Predicted vs measured, three ways",
        "The same predictions at three levels of averaging."),
}

# --- source parsing --------------------------------------------------------

def parse_cells(text: str) -> list[dict]:
    cells: list[dict] = []
    cur: dict | None = None
    for line in text.split("\n"):
        if line.startswith("# %%"):
            if cur is not None:
                cells.append(cur)
            m = re.search(r'tags=\[([^\]]*)\]', line)
            tags = ([t.strip().strip('"\'') for t in m.group(1).split(",") if t.strip()]
                    if m else [])
            cur = {"markdown": "[markdown]" in line, "tags": tags, "lines": []}
        elif cur is not None:
            cur["lines"].append(line)
    if cur is not None:
        cells.append(cur)
    for c in cells:
        body = "\n".join(c["lines"])
        if c["markdown"]:
            body = "\n".join(l[2:] if l.startswith("# ") else l.lstrip("#")
                             for l in c["lines"])
        c["text"] = body.strip("\n")
    return cells


def cell_stdout(nb_path: Path, n_code: int) -> list[str]:
    """Printed output of each code cell, in order, from the executed notebook."""
    if not nb_path.exists():
        return [""] * n_code
    nb = json.loads(nb_path.read_text())
    out = []
    for c in nb["cells"]:
        if c["cell_type"] != "code":
            continue
        txt = "".join("".join(o.get("text", []))
                      for o in c.get("outputs", [])
                      if o.get("output_type") == "stream")
        # "saved: fig6_pls.{png,pdf,svg}" is about the build, not the analysis.
        txt = "\n".join(l for l in txt.split("\n")
                        if not l.startswith("saved: "))
        out.append(txt.rstrip())
    return out + [""] * max(0, n_code - len(out))


# --- a very small markdown subset -----------------------------------------

def md_inline(s: str) -> str:
    s = html.escape(s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![*\w])\*([^*\n]+)\*(?!\w)", r"<em>\1</em>", s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
               r'<a href="\2" rel="noreferrer">\1</a>', s)
    return s


def strip_title(text: str) -> str:
    """Remove a leading markdown heading; the page supplies its own title."""
    lines = text.split("\n")
    for i, ln in enumerate(lines):
        if not ln.strip():
            continue
        return "\n".join(lines[i + 1:]).lstrip("\n") if ln.startswith("#") \
            else text
    return text


def md(text: str) -> str:
    """Headings, paragraphs, lists, tables, rules. Enough for these captions."""
    lines = text.split("\n")
    out, i = [], 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1
            continue
        if re.match(r"^-{3,}$", ln.strip()):
            out.append("<hr>")
            i += 1
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", ln)
        if m:
            lvl = min(len(m.group(1)) + 1, 6)
            out.append(f"<h{lvl}>{md_inline(m.group(2))}</h{lvl}>")
            i += 1
            continue
        if ln.lstrip().startswith("|") and i + 1 < len(lines) and \
                re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
            head = [c.strip() for c in ln.strip().strip("|").split("|")]
            i += 2
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append([c.strip() for c in
                             lines[i].strip().strip("|").split("|")])
                i += 1
            th = "".join(f"<th>{md_inline(c)}</th>" for c in head)
            tb = "".join("<tr>" + "".join(f"<td>{md_inline(c)}</td>"
                                          for c in r) + "</tr>" for r in rows)
            out.append(f'<div class="scroll"><table><thead><tr>{th}</tr>'
                       f"</thead><tbody>{tb}</tbody></table></div>")
            continue
        if re.match(r"^\s*[-*]\s+", ln) or re.match(r"^\s*\d+\.\s+", ln):
            ordered = bool(re.match(r"^\s*\d+\.\s+", ln))
            items, cur = [], ""
            while i < len(lines):
                l2 = lines[i]
                if re.match(r"^\s*[-*]\s+", l2) or re.match(r"^\s*\d+\.\s+", l2):
                    if cur:
                        items.append(cur)
                    cur = re.sub(r"^\s*(?:[-*]|\d+\.)\s+", "", l2)
                elif l2.strip() and (l2.startswith("  ") or l2.startswith("\t")):
                    cur += " " + l2.strip()
                else:
                    break
                i += 1
            if cur:
                items.append(cur)
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{md_inline(x)}</li>"
                                            for x in items) + f"</{tag}>")
            continue
        para = []
        while i < len(lines) and lines[i].strip() and \
                not lines[i].lstrip().startswith(("|", "#", "- ", "* ")) and \
                not re.match(r"^\s*\d+\.\s+", lines[i]) and \
                not re.match(r"^-{3,}$", lines[i].strip()):
            para.append(lines[i].strip())
            i += 1
        if para:
            out.append(f"<p>{md_inline(' '.join(para))}</p>")
        else:
            i += 1
    return "\n".join(out)


#: Longest edge, in pixels, of the figure embedded in the page. The renders are
#: 2100-3100 px wide; at full size the ten of them make a ~10 MB page, which is
#: a slow load for something meant to be opened from a link. 1900 keeps the
#: heatmap cells and tick labels legible at the width the page shows them. The
#: PDF and SVG in `notebooks/outputs/` remain the publication-resolution copies.
WEB_MAX_PX = 1900


def vector_payload(stem: str) -> dict:
    """The downloadable PDF and SVG for one figure, base64 for the page.

    The PDF ships as-is; PDF is already a compressed container and gzip buys
    nothing. The SVG ships **gzipped**, because matplotlib SVG is verbose XML
    that compresses about sevenfold here (7.0 MB of SVG becomes 1.0 MB), and
    the whole page has a 16 MB ceiling. The page inflates it with
    `DecompressionStream` before handing the bytes to the viewer.
    """
    import gzip
    out = {}
    pdf = OUTPUTS / f"{stem}.pdf"
    svg = OUTPUTS / f"{stem}.svg"
    if pdf.exists():
        out["pdf"] = base64.b64encode(pdf.read_bytes()).decode()
    if svg.exists():
        out["svgz"] = base64.b64encode(
            gzip.compress(svg.read_bytes(), 9)).decode()
    return out


def data_uri(p: Path) -> str:
    """Base64 PNG, downscaled and palette-quantised for the web.

    Resampling alone makes these files BIGGER, not smaller: the renders are
    flat-shaded heatmaps with few distinct colours, and LANCZOS interpolation
    invents thousands of intermediate ones (3.9 MB of source PNG becomes
    5.6 MB). Quantising back to a 256-colour palette after the resize is what
    actually shrinks them, to about 2.2 MB, and on these figures it is visually
    lossless — checked against the diverging W matrix in Figure 6, the most
    gradient-sensitive panel in the set, with no banding.
    """
    raw = p.read_bytes()
    try:
        from PIL import Image
        import io
    except ImportError:
        return "data:image/png;base64," + base64.b64encode(raw).decode()
    im = Image.open(io.BytesIO(raw)).convert("RGB")
    if max(im.size) > WEB_MAX_PX:
        scale = WEB_MAX_PX / max(im.size)
        im = im.resize((round(im.size[0] * scale), round(im.size[1] * scale)),
                       Image.LANCZOS)
    buf = io.BytesIO()
    im.quantize(colors=256, method=Image.MEDIANCUT).save(
        buf, "PNG", optimize=True)
    if len(buf.getvalue()) < len(raw):
        raw = buf.getvalue()
    return "data:image/png;base64," + base64.b64encode(raw).decode()


# --- page assembly ---------------------------------------------------------

def build_pages(cells: list[dict], stdout: list[str]) -> list[dict]:
    """Group cells into pages: setup, then one per figure, in notebook order."""
    pages: list[dict] = []
    cur = {"stem": None, "human": [], "llm": [], "blocks": []}
    code_i = 0
    for c in cells:
        if c["markdown"]:
            # A new numbered section closes an open page that has code but no
            # figure — that is §1, the setup, which saves nothing and would
            # otherwise be swallowed into Figure 3's page.
            if re.match(r"^##\s+\d+\.", c["text"]) and cur["blocks"] \
                    and cur["stem"] is None:
                # Code with no figure of its own. Before any figure it is §1,
                # the setup, and becomes its own page. After one it is a
                # trailing check or printout — the Fig 6 spectrum assertion,
                # Fig 5's spread numbers — which belongs to the figure above
                # it, not to a page of its own and certainly not dropped.
                if pages and pages[-1]["stem"]:
                    pages[-1]["blocks"] += cur["blocks"]
                    pages[-1]["human"] += cur["human"]
                    pages[-1]["llm"] += cur["llm"]
                else:
                    pages.append(cur)
                cur = {"stem": None, "human": [], "llm": [], "blocks": []}
            if "human" in c["tags"]:
                cur["human"].append(c["text"])
            elif "llm" in c["tags"]:
                cur["llm"].append(c["text"])
            else:
                cur["human"].append(c["text"])
                cur["llm"].append(c["text"])
            continue
        code = c["text"]
        out = stdout[code_i] if code_i < len(stdout) else ""
        code_i += 1
        if out.strip():
            cur["blocks"].append(out.rstrip())
        m = re.search(r'save\(fig,\s*"([a-z0-9_]+)"', code)
        if m:
            cur["stem"] = m.group(1)
            pages.append(cur)
            cur = {"stem": None, "human": [], "llm": [], "blocks": []}
    if cur["blocks"] or cur["human"]:
        # Whatever is left at the end — the closing prose cell, and any final
        # check — belongs to the last figure.
        if pages:
            pages[-1]["blocks"] += cur["blocks"]
            pages[-1]["human"] += cur["human"]
            pages[-1]["llm"] += cur["llm"]
        else:
            pages.append(cur)
    return pages


def render(pages: list[dict], images: dict[str, str]) -> str:
    nav, panes = [], []
    setup_seen = False
    for idx, pg in enumerate(pages):
        stem = pg["stem"]
        if stem is None and not setup_seen:
            setup_seen = True
            pid, label, name = "setup", "Setup", "Data and model"
            summary = ("Loads the imaging data and fits the model the figure "
                       "sections share.")
        elif stem in FIGURES:
            label, name, summary = FIGURES[stem]
            pid = stem
        else:
            continue
        pg["_id"], pg["_label"], pg["_name"] = pid, label, name

        nav.append(
            f'<li><a class="navlink" href="#{pid}" data-page="{pid}">'
            f'<span class="navnum">{html.escape(label)}</span>'
            f'<span class="navname">{html.escape(name)}</span></a></li>')

        img = ""
        if stem and stem in images:
            img = (f'<figure class="fig">'
                   f'<a class="figlink" href="{images[stem]}" target="_blank" '
                   f'rel="noreferrer">'
                   f'<img src="{images[stem]}" alt="{html.escape(name)}" '
                   f'loading="lazy"></a>'
                   f'<figcaption><span>{html.escape(label)} '
                   f"&middot; <code>{stem}</code> &middot; click to enlarge"
                   f'</span><span class="dl" data-stem="{stem}">'
                   f'<button class="dlbtn" data-stem="{stem}" data-fmt="pdf">'
                   f"Download PDF</button>"
                   f'<button class="dlbtn" data-stem="{stem}" data-fmt="svg">'
                   f"Download SVG</button>"
                   f'<span class="dlmsg" role="status" aria-live="polite">'
                   f"</span></span></figcaption></figure>")

        # Each caption opens with its own "## 4. Figure 4 - ..." heading, which
        # the page already shows as its title. Drop it so the title is not
        # printed twice under itself.
        human = md(strip_title("\n\n".join(pg["human"]))) if pg["human"] else ""
        grisly = md(strip_title("\n\n".join(pg["llm"]))) if pg["llm"] else ""
        # The notebook's printed numbers are results, not code, so they stay
        # -- inside the details, where the rest of the provenance lives.
        numbers = ""
        if pg["blocks"]:
            joined = "\n".join(pg["blocks"])
            numbers = ('<h4 class="numhead">Numbers this section printed</h4>'
                       f'<pre class="out">{html.escape(joined)}</pre>')
        details = ""
        if grisly or numbers:
            details = (
                '<details class="grisly"><summary>The grisly details</summary>'
                f'<div class="grisly-body">{grisly}{numbers}</div></details>')

        panes.append(
            f'<section class="page" id="page-{pid}" data-page="{pid}">'
            f'<p class="eyebrow">{html.escape(label)}</p>'
            f"<h1>{html.escape(name)}</h1>"
            f'<p class="lede">{html.escape(summary)}</p>'
            f"{img}"
            f'<div class="prose">{human}</div>'
            f"{details}"
            "</section>")

    toc_rows = "".join(
        f'<a class="tocrow" href="#{p["_id"]}" data-page="{p["_id"]}">'
        f'<span class="tocnum">{html.escape(p["_label"])}</span>'
        f'<span class="toctext"><strong>{html.escape(p["_name"])}</strong>'
        f'<span>{html.escape(FIGURES[p["stem"]][2] if p["stem"] in FIGURES else "Loads the imaging data and fits the model the figure sections share.")}</span>'
        "</span></a>"
        for p in pages if p.get("_id"))

    files = {p["stem"]: vector_payload(p["stem"])
             for p in pages if p.get("stem")}
    return TEMPLATE.format(nav="".join(nav), panes="".join(panes),
                           toc=toc_rows,
                           files=json.dumps(files, separators=(",", ":")))


TEMPLATE = """<title>Concentration Broadening Atlas</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,600;1,6..72,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --ground:#faf8f6; --panel:#ffffff; --rail:#f3efec;
  --ink:#1b1614; --ink-2:#4a403c; --ink-3:#7d716c;
  --rule:#e2dbd6; --rule-2:#d3cac4;
  --blue:#0072b2; --vermilion:#d55e00; --teal:#009e73; --amber:#e69f00;
  --code-bg:#f6f2ef; --code-ink:#2a2320;
  --k:#a63d00; --s:#0a6b52; --c:#8d817b; --n:#0060a0;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground:#14100f; --panel:#1b1615; --rail:#191413;
    --ink:#f2ebe7; --ink-2:#c3b7b1; --ink-3:#8e807a;
    --rule:#2e2523; --rule-2:#3d322f;
    --blue:#4ea6de; --vermilion:#ff8a43; --teal:#3ecfa5; --amber:#f0bd52;
    --code-bg:#211a18; --code-ink:#e8ded9;
    --k:#ff9d6b; --s:#6fd7b4; --c:#8a7a74; --n:#7cc0ef;
  }}
}}
:root[data-theme="dark"] {{
  --ground:#14100f; --panel:#1b1615; --rail:#191413;
  --ink:#f2ebe7; --ink-2:#c3b7b1; --ink-3:#8e807a;
  --rule:#2e2523; --rule-2:#3d322f;
  --blue:#4ea6de; --vermilion:#ff8a43; --teal:#3ecfa5; --amber:#f0bd52;
  --code-bg:#211a18; --code-ink:#e8ded9;
  --k:#ff9d6b; --s:#6fd7b4; --c:#8a7a74; --n:#7cc0ef;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--ground); color:var(--ink);
  font-family:"IBM Plex Sans",-apple-system,system-ui,sans-serif;
  font-size:16px; line-height:1.62; -webkit-font-smoothing:antialiased;
}}
.shell {{ display:grid; grid-template-columns:19rem minmax(0,1fr); min-height:100vh; }}
/* ---- rail ---- */
.rail {{
  background:var(--rail); border-right:1px solid var(--rule);
  padding:2rem 1.25rem 3rem; position:sticky; top:0; height:100vh;
  overflow-y:auto; display:flex; flex-direction:column; gap:1.4rem;
}}
.brand a {{ text-decoration:none; color:inherit; display:block; }}
.brand h2 {{
  font-family:Newsreader,Georgia,serif; font-weight:600; font-size:1.28rem;
  line-height:1.22; margin:0 0 .3rem; text-wrap:balance;
}}
.brand p {{ margin:0; font-size:.82rem; color:var(--ink-3); }}
.railhead {{
  font-size:.7rem; letter-spacing:.13em; text-transform:uppercase;
  color:var(--ink-3); margin:0; font-weight:600;
}}
.rail ul {{ list-style:none; margin:0; padding:0; display:flex; flex-direction:column; gap:.1rem; }}
.navlink {{
  display:grid; grid-template-columns:5.4rem 1fr; gap:.5rem; align-items:baseline;
  padding:.42rem .5rem; border-radius:3px; text-decoration:none; color:var(--ink-2);
  border-left:2px solid transparent;
}}
.navlink:hover {{ background:var(--panel); color:var(--ink); }}
.navlink.active {{
  background:var(--panel); color:var(--ink); border-left-color:var(--vermilion);
}}
.navnum {{
  font-family:"IBM Plex Mono",monospace; font-size:.7rem; color:var(--ink-3);
  letter-spacing:.02em; font-variant-numeric:tabular-nums;
}}
.navlink.active .navnum {{ color:var(--vermilion); }}
.navname {{ font-size:.87rem; line-height:1.32; }}
.railfoot {{ margin-top:auto; font-size:.74rem; color:var(--ink-3); line-height:1.5; }}
.railfoot code {{ font-size:.72rem; }}
/* ---- main ---- */
main {{ padding:3.2rem clamp(1.25rem,4vw,4rem) 6rem; max-width:66rem; }}
.page {{ display:none; }}
.page.active {{ display:block; animation:rise .34s ease both; }}
@keyframes rise {{ from {{ opacity:0; transform:translateY(6px); }} to {{ opacity:1; transform:none; }} }}
@media (prefers-reduced-motion:reduce) {{ .page.active {{ animation:none; }} }}
.eyebrow {{
  font-family:"IBM Plex Mono",monospace; font-size:.73rem; letter-spacing:.14em;
  text-transform:uppercase; color:var(--vermilion); margin:0 0 .55rem; font-weight:500;
}}
h1 {{
  font-family:Newsreader,Georgia,serif; font-weight:600; font-size:clamp(1.9rem,3.6vw,2.7rem);
  line-height:1.14; margin:0 0 .7rem; letter-spacing:-.012em; text-wrap:balance;
}}
.lede {{
  font-size:1.06rem; color:var(--ink-2); margin:0 0 2.2rem; max-width:60ch;
  border-left:2px solid var(--rule-2); padding-left:1rem;
}}
/* ---- figure ---- */
.fig {{ margin:0 0 2.4rem; }}
.figlink {{ display:block; background:#fff; border:1px solid var(--rule);
  border-radius:2px; padding:1.1rem; }}
.figlink:focus-visible {{ outline:2px solid var(--blue); outline-offset:3px; }}
.fig img {{ display:block; width:100%; height:auto; }}
.fig figcaption {{
  margin-top:.6rem; font-size:.78rem; color:var(--ink-3);
  font-family:"IBM Plex Mono",monospace;
  display:flex; flex-wrap:wrap; gap:.6rem 1rem;
  align-items:center; justify-content:space-between;
}}
.dl {{ display:flex; gap:.45rem; }}
.dlbtn {{
  font-family:"IBM Plex Mono",monospace; font-size:.72rem; letter-spacing:.02em;
  color:var(--blue); background:none; border:1px solid var(--rule-2);
  border-radius:2px; padding:.28rem .6rem; cursor:pointer;
  transition:background .14s,border-color .14s;
}}
.dlbtn:hover {{ background:var(--panel); border-color:var(--blue); }}
.dlbtn:focus-visible {{ outline:2px solid var(--blue); outline-offset:2px; }}
.dlbtn[disabled] {{ opacity:.5; cursor:default; }}
.dlmsg {{ color:var(--vermilion); font-size:.72rem; max-width:30ch; }}
/* ---- prose ---- */
.prose {{ max-width:66ch; }}
.prose h2 {{ font-family:Newsreader,Georgia,serif; font-size:1.42rem; font-weight:600;
  margin:2.2rem 0 .6rem; line-height:1.25; text-wrap:balance; }}
.prose h3,.prose h4 {{ font-family:"IBM Plex Sans",sans-serif; font-size:.79rem;
  letter-spacing:.11em; text-transform:uppercase; color:var(--ink-3);
  margin:2rem 0 .5rem; font-weight:600; }}
.prose p {{ margin:0 0 1.05rem; }}
.prose ul,.prose ol {{ margin:0 0 1.15rem; padding-left:1.15rem; }}
.prose li {{ margin-bottom:.5rem; }}
.prose strong {{ font-weight:600; color:var(--ink); }}
.prose a {{ color:var(--blue); text-underline-offset:2px; }}
.prose hr {{ border:0; border-top:1px solid var(--rule); margin:2rem 0; }}
code {{ font-family:"IBM Plex Mono",monospace; font-size:.855em;
  background:var(--code-bg); padding:.1em .34em; border-radius:2px; color:var(--code-ink); }}
table {{ border-collapse:collapse; width:100%; font-size:.88rem; margin:0 0 1.2rem; }}
th,td {{ text-align:left; padding:.44rem .7rem; border-bottom:1px solid var(--rule); }}
th {{ font-size:.7rem; letter-spacing:.09em; text-transform:uppercase; color:var(--ink-3); }}
td {{ font-variant-numeric:tabular-nums; }}
.scroll {{ overflow-x:auto; }}
/* ---- disclosures ---- */
details {{ margin:0 0 .7rem; }}
summary {{ cursor:pointer; list-style:none; }}
summary::-webkit-details-marker {{ display:none; }}
summary:focus-visible {{ outline:2px solid var(--blue); outline-offset:2px; }}
.grisly {{ margin:2.2rem 0 2.6rem; max-width:66ch; }}
.grisly > summary {{
  display:flex; align-items:center; gap:.55rem; font-size:.82rem; font-weight:600;
  letter-spacing:.05em; text-transform:uppercase; color:var(--teal);
  padding:.6rem .85rem; border:1px solid var(--rule);
  border-left:2px solid var(--teal); border-radius:2px; background:var(--panel);
}}
/* CSS unicode escape, not the literal glyph: the whole file is written as
   ASCII with numeric character references, and CSS does not decode HTML
   entities -- a literal here would render as the text "&#9656;". */
.grisly > summary::before {{ content:"\\25B8"; font-size:.8em; transition:transform .16s; }}
.grisly[open] > summary::before {{ transform:rotate(90deg); }}
.grisly-body {{
  padding:1.3rem .95rem .2rem; border-left:2px solid var(--rule-2);
  margin:.7rem 0 0 .1rem; font-size:.95rem; color:var(--ink-2);
}}
.grisly-body p,.grisly-body li {{ margin-bottom:.85rem; }}
.grisly-body strong {{ color:var(--ink); }}
pre {{ margin:.3rem 0 0; }}
.numhead {{
  font-family:"IBM Plex Sans",sans-serif; font-size:.7rem; letter-spacing:.12em;
  text-transform:uppercase; color:var(--ink-3); font-weight:600;
  margin:1.8rem 0 .55rem;
}}
.out {{
  border-left:2px solid var(--amber); padding:.6rem .95rem; margin:0;
  font-family:"IBM Plex Mono",monospace; font-size:.755rem; color:var(--ink-2);
  line-height:1.6; white-space:pre-wrap; background:var(--panel);
  overflow-x:auto; font-variant-numeric:tabular-nums;
}}
/* ---- contents page ---- */
.tocrow {{
  display:grid; grid-template-columns:7rem 1fr; gap:1rem; padding:1rem .6rem;
  border-bottom:1px solid var(--rule); text-decoration:none; color:inherit;
  align-items:baseline;
}}
.tocrow:hover {{ background:var(--panel); }}
.tocrow:focus-visible {{ outline:2px solid var(--blue); outline-offset:-2px; }}
.tocnum {{
  font-family:"IBM Plex Mono",monospace; font-size:.75rem; color:var(--vermilion);
  letter-spacing:.03em;
}}
.toctext {{ display:flex; flex-direction:column; gap:.18rem; }}
.toctext strong {{ font-size:1.02rem; font-weight:600; }}
.toctext span {{ font-size:.88rem; color:var(--ink-3); line-height:1.5; }}
.pager {{
  display:flex; justify-content:space-between; gap:1rem; margin-top:3.5rem;
  padding-top:1.2rem; border-top:1px solid var(--rule); font-size:.87rem;
}}
.pager a {{ color:var(--blue); text-decoration:none; }}
.pager a:hover {{ text-decoration:underline; }}
@media (max-width:900px) {{
  .shell {{ grid-template-columns:1fr; }}
  .rail {{ position:static; height:auto; }}
  main {{ padding-top:1.6rem; }}
}}
</style>

<div class="shell">
  <nav class="rail">
    <div class="brand">
      <a href="#contents" data-page="contents">
        <h2>Concentration broadening in the olfactory bulb</h2>
        <p>Figure notebook · Castro Lab</p>
      </a>
    </div>
    <p class="railhead">Contents</p>
    <ul>
      <li><a class="navlink" href="#contents" data-page="contents">
        <span class="navnum">—</span><span class="navname">Contents</span></a></li>
      {nav}
    </ul>
    <p class="railfoot">
      Every figure is drawn by one notebook, <code>notebooks/figures.py</code>.
      Analysis lives in the <code>highconc</code> package; each panel is a
      function that takes an axes. Groups are numbered 1–7 throughout.
    </p>
  </nav>

  <main>
    <section class="page" id="page-contents" data-page="contents">
      <p class="eyebrow">Figure notebook</p>
      <h1>Concentration broadening in the olfactory bulb</h1>
      <p class="lede">Every figure in the manuscript, with the reasoning and the
        code that produced it. 376 glomeruli, 59 odorants, two concentrations,
        four mice.</p>
      <div class="prose">
        <p>Raising odorant concentration does not simply make glomerular
        responses larger. It changes which odorants evoke similar response
        patterns, and it broadens tuning selectively — particular groups of
        glomeruli take on the tuning of particular others, and not the
        reverse.</p>
        <p>Each page below carries one figure and a plain-language account of
        what it shows and how it was made. Under that, a collapsible
        <strong>grisly details</strong> block holds the full provenance: every
        parameter, every value, every judgement call, and the numbers the
        analysis printed.</p>
      </div>
      {toc}
    </section>
    {panes}
  </main>
</div>

<script id="vectors" type="application/json">{files}</script>
<script>
(function () {{
  // ---- downloads -------------------------------------------------------
  // This page has to work in two places, and they need different mechanisms.
  //
  // Inside the claude.ai artifact viewer the frame is sandboxed and blocks any
  // download the page starts itself, so an <a download> is inert there; the
  // `downloads` capability is the supported route and shows the viewer a
  // confirmation. Served from anywhere else -- a lab server, GitHub Pages, or
  // straight off disk -- there is no capability, and a blob <a download> is
  // both available and better: no prompt, no allowlist, nothing to enable.
  //
  // So the buttons are always visible, and the mechanism is chosen when one is
  // clicked: capability if this view has it, plain link otherwise.
  var FILES = JSON.parse(document.getElementById("vectors").textContent);
  var saverP = (window.claude && claude.use)
    ? claude.use("downloads").catch(function () {{ return null; }})
    : Promise.resolve(null);

  function bytes(b64) {{
    var bin = atob(b64), n = bin.length, u = new Uint8Array(n);
    for (var i = 0; i < n; i++) u[i] = bin.charCodeAt(i);
    return u;
  }}

  // SVG is stored gzipped to fit the page budget; inflate before saving.
  function inflate(u8) {{
    if (typeof DecompressionStream === "undefined") {{
      return Promise.reject({{ code: "no_inflate" }});
    }}
    var s = new Blob([u8]).stream().pipeThrough(new DecompressionStream("gzip"));
    return new Response(s).arrayBuffer();
  }}

  function payload(stem, fmt) {{
    var f = FILES[stem] || {{}};
    if (fmt === "pdf") return Promise.resolve(bytes(f.pdf));
    return inflate(bytes(f.svgz));
  }}

  var TYPES = {{ pdf: "application/pdf", svg: "image/svg+xml" }};

  function viaLink(stem, fmt, data) {{
    var url = URL.createObjectURL(new Blob([data], {{ type: TYPES[fmt] }}));
    var a = document.createElement("a");
    a.href = url;
    a.download = stem + "." + fmt;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(function () {{ URL.revokeObjectURL(url); }}, 30000);
  }}

  var MESSAGES = {{
    declined: "Cancelled.",
    extension_not_enabled: "This viewer will not accept " +
      "PDF or SVG. Open the page in a browser tab instead.",
    rejected_extension: "This viewer will not accept that file type.",
    too_large: "Too large for the viewer to save.",
    rate_limited: "A save prompt is already open.",
    unavailable: "Saving is not available in this view.",
    no_inflate: "This browser cannot unpack the SVG. Try the PDF."
  }};

  function wire() {{
    document.querySelectorAll(".dlbtn").forEach(function (b) {{
      b.addEventListener("click", function () {{
        var msg = b.parentNode.querySelector(".dlmsg");
        var stem = b.dataset.stem, fmt = b.dataset.fmt;
        b.disabled = true;
        msg.textContent = "";
        Promise.all([payload(stem, fmt), saverP]).then(function (r) {{
          var data = r[0], saver = r[1];
          if (!saver) {{ viaLink(stem, fmt, data); return; }}
          return saver.save({{ filename: stem + "." + fmt, data: data }});
        }}).catch(function (e) {{
          msg.textContent = MESSAGES[e && e.code] || "Could not save that file.";
        }}).then(function () {{ b.disabled = false; }});
      }});
    }});
  }}
  wire();

  // ---- paging ----------------------------------------------------------
  var pages = document.querySelectorAll(".page");
  var links = document.querySelectorAll("[data-page]");
  var order = Array.prototype.map.call(pages, function (p) {{
    return p.dataset.page;
  }});

  function show(id) {{
    var found = false;
    pages.forEach(function (p) {{
      var on = p.dataset.page === id;
      p.classList.toggle("active", on);
      if (on) found = true;
    }});
    if (!found) return show("contents");
    document.querySelectorAll(".navlink").forEach(function (a) {{
      a.classList.toggle("active", a.dataset.page === id);
    }});
    var i = order.indexOf(id);
    var pane = document.getElementById("page-" + id);
    var old = pane.querySelector(".pager");
    if (old) old.remove();
    var prev = i > 0 ? order[i - 1] : null;
    var next = i < order.length - 1 ? order[i + 1] : null;
    if (prev || next) {{
      var nav = document.createElement("div");
      nav.className = "pager";
      nav.innerHTML =
        (prev ? '<a href="#' + prev + '">\\u2190 previous</a>' : "<span></span>") +
        (next ? '<a href="#' + next + '">next \\u2192</a>' : "<span></span>");
      pane.appendChild(nav);
    }}
    window.scrollTo(0, 0);
  }}

  window.addEventListener("hashchange", function () {{
    show(location.hash.slice(1) || "contents");
  }});
  show(location.hash.slice(1) || "contents");
}})();
</script>
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    cells = parse_cells(SOURCE.read_text())
    n_code = sum(1 for c in cells if not c["markdown"])
    pages = build_pages(cells, cell_stdout(EXECUTED, n_code))

    missing = [s for s in FIGURES if not (OUTPUTS / f"{s}.png").exists()]
    if missing:
        print("missing renders — run the notebook first:", file=sys.stderr)
        for s in missing:
            print(f"    notebooks/outputs/{s}.png", file=sys.stderr)
        return 1

    images = {s: data_uri(OUTPUTS / f"{s}.png") for s in FIGURES}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    # Written as ASCII with numeric character references for everything else.
    # The page is one file served by hosts we do not control, and a missing
    # `charset=utf-8` turns every em-dash into mojibake. Entities render
    # correctly under any encoding and cost about 1% in size here.
    out.write_bytes(render(pages, images)
                    .encode("ascii", "xmlcharrefreplace"))
    mb = out.stat().st_size / 1048576
    print(f"{out.relative_to(REPO)}  {len(pages)} pages, "
          f"{len(images)} figures, {mb:.1f} MB\n")
    for pg in pages:
        if not pg.get("_id"):
            continue
        lines = sum(len(b.split("\n")) for b in pg["blocks"])
        print(f"  {pg['_label']:14} {pg['_name'][:38]:40} "
              f"details={'yes' if pg['llm'] else 'NO ':3}  "
              f"printed={lines:3d} lines")
    if mb > 15:
        print("  WARNING: approaching the 16 MB artifact ceiling", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
