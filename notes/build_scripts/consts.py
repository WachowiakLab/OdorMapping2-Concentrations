import re, sys
t = open(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py", encoding="utf-8").read().replace("\r\n", "\n")
cells = [c for c in re.split(r"(?m)^(?=# %%)", t) if c.startswith("# %%")]
sec, out = 0, {}
for c in cells:
    if c.startswith("# %% [markdown]"):
        m = re.search(r"(?m)^# ## (\d+)\.", c)
        if m:
            sec = int(m.group(1))
        continue
    for l in c.split("\n"):
        m = re.match(r"^([A-Z][A-Z0-9_]*(?:, *[A-Z][A-Z0-9_]*)*) *= *(.*)$", l)
        if m:
            val = re.sub(r"\s+#.*$", "", m.group(2))[:55]
            out.setdefault(sec, []).append(f"{m.group(1)} = {val}")
for s in sorted(out):
    print(f"--- §{s}")
    print("\n".join("   " + x for x in out[s]))
