"""Read-only exploration: fold nonresponders into low, recluster (seeded), compare to the dropped clustering."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.dont_write_bytecode = True
sys.path.insert(0, r"C:\Users\u6068031\Documents\highconc-merge\analysis\GlomExplorer")
import numpy as np, pandas as pd
import glom_explorer as gx

odorants, low_raw, high_raw = gx.load_processed_data()
print("raw low/high:", low_raw.shape, high_raw.shape, "| column levels:", list(low_raw.columns.names))

high_recruited = pd.concat([(low_raw == 0).all(0), (high_raw != 0).any()], axis=1).all(1)
print("glomeruli silent at low but responding at high ('recruited nonresponders'):", int(high_recruited.sum()), "of", low_raw.shape[1])

new_low_raw = low_raw.copy()
new_low_raw.loc[:, high_recruited.values] = high_raw.loc[:, high_recruited.values]

mask = (new_low_raw.std() == 0)
cluster_on = new_low_raw.loc[:, ~mask]
print("glomeruli kept for clustering (std>0):", cluster_on.shape[1], "| dropped (flat):", int(mask.sum()))
row_labels, col_labels = gx.assign_cluster_indices(cluster_on, method="nmf")   # seeded defaults (706231266, 1000 iter)
print("new cluster sizes (raw NMF labels):", dict(zip(*np.unique(col_labels, return_counts=True))))

# the original clustering
low_d, high_d = gx.load_clustered_variant("nonresponders_dropped")
print("dropped variant low:", low_d.shape, "| column levels:", list(low_d.columns.names))
key_old = list(zip(low_d.columns.get_level_values("subject").astype(str), low_d.columns.get_level_values("gloms").astype(str)))
key_new = list(zip(cluster_on.columns.get_level_values("subject").astype(str), cluster_on.columns.get_level_values("gloms").astype(str)))
print("unique keys: old", len(set(key_old)), "/", len(key_old), " new", len(set(key_new)), "/", len(key_new))
old_map = dict(zip(key_old, low_d.columns.get_level_values("parent_cluster")))
src = np.array([old_map.get(k, "NA") for k in key_new], dtype=object)
print("new glomeruli found in the dropped clustering:", int((src != "NA").sum()), "| not found:", int((src == "NA").sum()))
ct = pd.crosstab(pd.Series(src, name="old(raw NMF idx or NA)"), pd.Series(col_labels, name="new(raw NMF idx)"))
print(ct.to_string())
