import re,sys
from pathlib import Path
t=Path(sys.argv[1]).read_bytes().decode("utf-8").replace("\r\n","\n")
parts=[c for c in re.split(r"(?m)^(?=# %%)",t) if c.startswith("# %%")]
sec="0"
for i,c in enumerate(parts):
    md=c.startswith("# %% [markdown]")
    tag=re.search(r"tags=\[([^\]]*)\]",c.split("\n")[0]); tag=tag.group(1) if tag else ""
    m=re.search(r"(?m)^# ## (\d+)\.",c) if md else None
    if m: sec=m.group(1)
    body=[l for l in c.split("\n")[1:] if l.strip()]
    first=(body[0] if body else "")[:70]
    print(f"{i:3d} §{sec:>2} {'MD  ' if md else 'CODE'} {tag:6s} {len(c.split(chr(10))):4d}L  {first}")
