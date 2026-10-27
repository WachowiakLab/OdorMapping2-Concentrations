import json, sys, re
p = sys.argv[1]
nb = json.load(open(p, encoding="utf-8"))
cells = nb["cells"]
print("nbformat", nb.get("nbformat"), "| cells:", len(cells), "| kernel:", nb["metadata"].get("kernelspec", {}).get("name"))
for i, c in enumerate(cells):
    src = "".join(c["source"])
    first = next((l for l in src.split("\n") if l.strip()), "")[:95]
    flag = "  <<< SANKEY" if re.search(r"sankey", src, re.I) else ""
    n_out = len(c.get("outputs", [])) if c["cell_type"] == "code" else 0
    print(f"{i:3d} {c['cell_type'][:4]:4s} {len(src.splitlines()):4d}L out={n_out:2d} {first}{flag}")
