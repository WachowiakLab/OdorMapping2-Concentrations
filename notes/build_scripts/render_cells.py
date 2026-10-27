"""Execute chosen cells of final/notebooks/figures.py against the real data; renders go to final/notebooks/outputs.
usage: render_cells.py <comma-separated cell indices>
env:   PATCH="old=>new||old2=>new2"   in-memory text patches applied to each cell before running
       FAKE_M=1                       define `m` as a TEST STAND-IN model (plain sklearn PLS, NOT the paper's model)
                                      so layouts that need `m` can be rendered. Numbers are not meaningful."""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True
os.environ["MPLBACKEND"] = "Agg"
FINAL = r"C:\Users\u6068031\Documents\highconc-merge\final"
sys.path.insert(0, FINAL)
sys.path.insert(1, r"C:\Users\u6068031\Documents\highconc-merge\analysis\GlomExplorer")
from unittest.mock import MagicMock
sys.modules.setdefault("pls_explorer", MagicMock())
from pathlib import Path

if os.environ.get("FAKE_M"):
    # TEST-ONLY stand-ins for the two tiny pls_explorer helpers the renders reach.
    import types
    import numpy as _np

    def _reorder_by_clusters(W, cn):
        cn = _np.asarray(cn)
        order = _np.argsort(cn, kind="stable")
        return W[_np.ix_(order, order)], order, _np.where(_np.diff(cn[order]) != 0)[0] + 1

    def _q2_global(Y, Yhat):
        return 1.0 - ((Y - Yhat) ** 2).sum() / ((Y - Y.mean(axis=0)) ** 2).sum()

    def _rmse(Y, Yhat):
        return float(_np.sqrt(((Y - Yhat) ** 2).mean()))

    for name, attrs in (("pls_explorer.block_structure", {"reorder_by_clusters": _reorder_by_clusters}),
                        ("pls_explorer.metrics", {"q2_global": _q2_global, "rmse": _rmse})):
        mod = types.ModuleType(name)
        mod.__dict__.update(attrs)
        sys.modules[name] = mod

text = Path(FINAL, "notebooks", "figures.py").read_bytes().decode("utf-8").replace("\r\n", "\n")
cells = [c for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]
want = [int(x) for x in sys.argv[1].split(",")]
ns = {"__name__": "__main__", "__file__": str(Path(FINAL, "notebooks", "figures.py"))}


def make_fake_m(d):
    import numpy as np
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.model_selection import KFold
    from highconc import pls as _pls
    low = np.clip(d.low.to_numpy(float), 0, None)
    high = np.clip(d.high.to_numpy(float), 0, None)
    Xlog, Ylog = np.log1p(low), np.log1p(high)
    Yhat = np.zeros_like(Ylog)
    folds = []
    for tr, te in KFold(4, shuffle=True, random_state=0).split(Xlog):
        pl = PLSRegression(n_components=8, scale=False).fit(Xlog[tr], Ylog[tr])
        Yhat[te] = np.clip(pl.predict(Xlog[te]), 0, None)
        folds.append((tr, te))
    return _pls.PLSModel(low=low, high=high, Xlog=Xlog, Ylog=Ylog, Yhat=Yhat, folds=folds,
                         chem=np.asarray(d.chem_class), chem_pooled=np.asarray(d.chem_class),
                         names=d.odorants["Name"].to_numpy(), cluster_index=d.cluster_index,
                         cluster_number=d.cluster_number, xy=d.xy, subject=d.subject,
                         odorants=d.odorants, variant=d.variant)


for i in want:
    code = cells[i]
    if code.startswith("# %% [markdown]"):
        continue
    for pair in os.environ.get("PATCH", "").split("||"):
        if "=>" in pair:
            a, b = pair.split("=>"); code = code.replace(a, b)
    if os.environ.get("FAKE_M") and "m" not in ns and "d" in ns:
        ns["m"] = make_fake_m(ns["d"])
        print("[test stand-in model m created: NOT the paper's PLS]", flush=True)
    print(f"--- cell {i}", flush=True)
    exec(compile(code, f"<cell {i}>", "exec"), ns)
print("done")
