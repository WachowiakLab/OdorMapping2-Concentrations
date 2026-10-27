"""Compile figures.py as one script and report names that are read but never bound anywhere."""
import ast, builtins, sys, re
from pathlib import Path
t = Path(sys.argv[1]).read_text(encoding="utf-8")
tree = ast.parse(t)                       # SyntaxError if any cell is broken
bound = set(dir(builtins)) | {"__file__", "get_ipython", "display"}
for n in ast.walk(tree):
    if isinstance(n, (ast.Import, ast.ImportFrom)):
        for a in n.names: bound.add((a.asname or a.name).split(".")[0])
    elif isinstance(n, (ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)):
        bound.add(n.name)
        if not isinstance(n, ast.ClassDef):
            for a in n.args.args + n.args.kwonlyargs + n.args.posonlyargs: bound.add(a.arg)
            if n.args.vararg: bound.add(n.args.vararg.arg)
            if n.args.kwarg: bound.add(n.args.kwarg.arg)
    elif isinstance(n, ast.Lambda):
        for a in n.args.args: bound.add(a.arg)
    elif isinstance(n, ast.Name) and isinstance(n.ctx, (ast.Store, ast.Del)):
        bound.add(n.id)
    elif isinstance(n, ast.ExceptHandler) and n.name: bound.add(n.name)
    elif isinstance(n, ast.arg): bound.add(n.arg)
missing = {}
for n in ast.walk(tree):
    if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load) and n.id not in bound:
        missing.setdefault(n.id, n.lineno)
print("compiles: yes | names read but never bound:", missing or "none")
