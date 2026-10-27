"""Observed recruitment matrix under a given package's algorithm. usage: rec_variant.py <baseline|final> <out.npz>
Observed side uses only m.low / m.high (your measured data); the stand-in model is only needed to build `m`."""
import os, sys, warnings
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True
which, out = sys.argv[1], sys.argv[2]
ROOT = r"C:\Users\u6068031\Documents\highconc-merge"
pkg = ROOT + (r"\github\original_codebase" if which == "baseline" else r"\final")
sys.path.insert(0, pkg)
sys.path.insert(1, ROOT + r"\analysis\GlomExplorer")
from unittest.mock import MagicMock
sys.modules.setdefault("pls_explorer", MagicMock())
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import KFold
import highconc
assert highconc.__file__.startswith(pkg), highconc.__file__
from highconc.data import load
from highconc import pls as _pls, recruitment as rc

d = load("nonresponders_dropped")
low = np.clip(d.low.to_numpy(float), 0, None)
high = np.clip(d.high.to_numpy(float), 0, None)
Xlog, Ylog = np.log1p(low), np.log1p(high)
Yhat = np.zeros_like(Ylog)
for tr, te in KFold(4, shuffle=True, random_state=0).split(Xlog):
    Yhat[te] = np.clip(PLSRegression(8, scale=False).fit(Xlog[tr], Ylog[tr]).predict(Xlog[te]), 0, None)
m = _pls.PLSModel(low=low, high=high, Xlog=Xlog, Ylog=Ylog, Yhat=Yhat, folds=[], chem=np.asarray(d.chem_class),
                  chem_pooled=np.asarray(d.chem_class), names=d.odorants["Name"].to_numpy(),
                  cluster_index=d.cluster_index, cluster_number=d.cluster_number, xy=d.xy,
                  subject=d.subject, odorants=d.odorants, variant=d.variant)
codes = np.asarray([s.split()[0] for s in d.odorants["Odorant"].values])

res = {}
if which == "baseline":
    r = rc.recruitment_matrix(m, level="odorant", row_rule="halfmax", frac=0.50, order="matt",
                              prediction_threshold="rate_matched")
    res["baseline halfmax 0.50"] = r
else:
    for fr in (0.50, 0.20):
        r = rc.recruitment_matrix(m, level="odorant", row_rule="halfmax", frac=fr, order="matt",
                                  within="table", codes=codes, prediction_threshold="rate_matched")
        res[f"original rule halfmax {fr:.2f}"] = r
arrs = {}
for k, r in res.items():
    O = np.asarray(r.observed, float)
    arrs[k] = O
    off = ~np.eye(O.shape[0], dtype=bool)
    v = O[off & np.isfinite(O)]
    print(f"{k:34s} shape {O.shape} | max {np.nanmax(O):.3f} | mean(off-diag) {v.mean():.3f} | "
          f"NaN cells {int(np.isnan(O).sum())} | empty rows {int(r.n_empty)} | mean glom/row {float(np.mean(r.n_rows)):.1f}")
np.savez(out, **{k.replace(" ", "_").replace("=", ""): v for k, v in arrs.items()})
