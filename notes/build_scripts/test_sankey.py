import sys, warnings
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True
sys.path.insert(0, r"C:\Users\u6068031\Documents\highconc-merge\final")
sys.path.insert(1, r"C:\Users\u6068031\Documents\highconc-merge\analysis\GlomExplorer")
from unittest.mock import MagicMock
for n in ("pls_explorer",):
    sys.modules.setdefault(n, MagicMock())
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from highconc.data import load
from highconc import folding, panels as P

d = load("nonresponders_dropped")
fl = folding.flows(d)
print(fl)
np.set_printoptions(linewidth=140)
print("counts (rows: original 1-7, then not-in-original; cols: new 1-7 by best match)")
print(fl.counts)
print("original sizes", fl.original_sizes, "| new sizes", fl.new_sizes)
print("stayed in best-matched new cluster:", np.round(fl.stayed, 2))
print("raw NMF label -> number:", fl.label_to_number)

fig, ax = plt.subplots(figsize=(7, 2.2), dpi=200)
P.sankey.flow_diagram(ax, fl)
out = r"C:\Users\u6068031\AppData\Local\Temp\claude\C--Users-u6068031-Documents-highconc-merge\05160cad-cae0-44a6-ada4-e804c5fd776b\scratchpad\sankey_test.png"
fig.savefig(out, bbox_inches="tight")
print("wrote", out)

