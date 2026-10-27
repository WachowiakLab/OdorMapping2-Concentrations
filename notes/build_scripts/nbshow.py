import json, sys
p = sys.argv[1]; idx = [int(x) for x in sys.argv[2].split(",")]
nb = json.load(open(p, encoding="utf-8"))
for i in idx:
    c = nb["cells"][i]; src = "".join(c["source"])
    live = sum(1 for l in src.split("\n") if l.strip() and not l.lstrip().startswith("#"))
    tot = sum(1 for l in src.split("\n") if l.strip())
    print(f"\n######## cell {i} [{c['cell_type']}] non-comment lines: {live}/{tot}")
    print("\n".join(l[:150] for l in src.split("\n")))
