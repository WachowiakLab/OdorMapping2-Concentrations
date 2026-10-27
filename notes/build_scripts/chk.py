import re,sys,difflib
from pathlib import Path
def cells(p):
    b=Path(p).read_bytes()
    try: t=b.decode("utf-8")
    except UnicodeDecodeError:
        print("   !! not UTF-8, decoding as cp1252:",Path(p).parent.parent.name,Path(p).name); t=b.decode("cp1252")
    t=t.replace("\r\n","\n")
    parts=re.split(r"(?m)^(?=# %%)",t)
    out=[]
    for c in parts:
        if not c.startswith("# %%"): continue
        md=c.startswith("# %% [markdown]")
        tag=re.search(r"tags=\[([^\]]*)\]",c.split("\n")[0])
        out.append((md,tag.group(1) if tag else "",c.rstrip("\n")))
    return out
def sections(p):
    t=Path(p).read_text(encoding="utf-8").replace("\r\n","\n")
    return re.findall(r"(?m)^# ## (\d+\. .*)$",t)
for name,d in [("baseline","github/original_codebase"),("lab","highconc-figures-lab"),("main","highconc-figures-lab-main"),("copy","highconc-figures-lab-main - Copy")]:
    root=Path(sys.argv[1])/d/"notebooks"
    src=cells(root/"figures.py")
    code_src=[c[2] for c in src if not c[0]]
    print(f"== {name}: {len(src)} cells, {len(code_src)} code")
    seen=[]
    for s in sections(root/"figures.py"):
        if s not in seen: seen.append(s)
    print("   sections:", " | ".join(x[:48] for x in seen))
    for ed,fn in [("human","figures_human.py"),("llm","figures_LLM.py")]:
        p=root/fn
        if not p.exists(): print("   ",fn,"missing"); continue
        code_ed=[c[2] for c in cells(p) if not c[0]]
        if code_ed==code_src: print(f"   {fn}: code identical to figures.py")
        else:
            n=sum(1 for l in difflib.unified_diff("\n".join(code_src).split("\n"),"\n".join(code_ed).split("\n"),lineterm="",n=0) if l[:1] in "+-" and l[:3] not in("+++","---"))
            print(f"   {fn}: code DIFFERS from figures.py ({len(code_src)} vs {len(code_ed)} code cells, {n} changed lines)")
