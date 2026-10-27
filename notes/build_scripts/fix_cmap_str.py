"""panels/transfer.py: accept a colormap NAME as well as a Colormap in pvalue_matrix; then exercise both panels."""
import sys
from pathlib import Path
from types import SimpleNamespace

F = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")
p = F / "highconc" / "panels" / "transfer.py"
t = p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def sub(old, new):
    global t
    assert t.count(old) == 1, (t.count(old), old[:60])
    t = t.replace(old, new)


sub("import numpy as np\nfrom matplotlib.colors import Normalize\n",
    "import matplotlib as mpl\nimport numpy as np\nfrom matplotlib.colors import Normalize\n")
sub("    cmap = hs.SEQ if cmap is None else cmap\n    sig = t.significant(alpha)",
    "    cmap = hs.SEQ if cmap is None else cmap\n"
    "    if isinstance(cmap, str):            # a name such as \"OrRd\", not just a Colormap\n"
    "        cmap = mpl.colormaps[cmap]\n"
    "    sig = t.significant(alpha)")
sub("""    `cmap` is the hexbin colormap""" if False else "    The box, not the colour, is the claim:",
    "    `cmap` is a Colormap or its name (default `highconc.style.SEQ`).\n\n"
    "    The box, not the colour, is the claim:")
p.write_bytes(t.encode("utf-8"))
print("patched", p.name)

# ---- exercise it: a stand-in for highconc.transfer.Transfer, real matplotlib, no data needed ----
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(F))
from unittest.mock import MagicMock
for n in ("glom_explorer", "glom_explorer.topography", "glom_explorer.tuning", "pls_explorer", "plotly"):
    sys.modules.setdefault(n, MagicMock())
from highconc import panels as P

K = 7
rng = np.random.default_rng(0)
Pm = rng.uniform(0.001, 0.6, (K, K))
sig = np.zeros((K, K), bool); sig[0, 1] = sig[3, 2] = sig[5, 6] = True
off = ~np.eye(K, dtype=bool)
tt = SimpleNamespace(
    grouping=SimpleNamespace(n_groups=K), labels=[str(i + 1) for i in range(K)],
    P=Pm, p_floor=1e-4, n_perm=10000, off=off,
    significant=lambda alpha: sig,
    pct_err=rng.normal(0, 8, (K, K)), pct_var=rng.normal(0, 8, (K, K)), M=rng.normal(0, .02, (K, K)),
)
for cm in (None, "OrRd"):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 3.2))
    P.transfer.pvalue_matrix(a, tt, alpha=0.05, clip=0.1, xlabel="x", ylabel="y",
                             cmap=cm, box_color="#2b9743", text_significant_only=True)
    P.transfer.improvement_matrix(b, tt, metric="pct_err", alpha=0.05, xlabel="x", ylabel="",
                                  box_color="#2b9743", text_significant_only=True)
    P.transfer.significance_legend(b, tt, alpha=0.05, loc=(1.02, -0.16), box_color="#2b9743")
    n_text_a = sum(1 for tx in a.texts if tx.get_text())
    n_text_b = sum(1 for tx in b.texts if tx.get_text())
    print(f"cmap={cm!r}: drew OK | numbers printed: A={n_text_a}, B={n_text_b} (3 significant cells expected)")
    plt.close(fig)
