"""Modul asılılıq qrafiki — `imex2d/` paketinin idxal əlaqələri (Seans 44).

    python tools/module_graph.py                       graph.json cari qovluğa
    python tools/module_graph.py --out docs/x.json

Hər `.py` faylı bir node-dur (boş `__init__.py`-lər atlanır); kənar —
paket daxilindəki `import` / `from … import`. Node-da modulun docstring-inin
ilk sətri, sətir sayı, açıq siniflər və funksiyalar saxlanılır.

Təqdimat üçün interaktiv xəritənin məlumat mənbəyidir (bax
`docs/teqdimat/README.md`). Kodu İŞLƏTMİR — yalnız AST ilə oxuyur.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE = "imex2d"


def _modules(package_dir: str) -> dict:
    """`{modul.adı: fayl_yolu}` — `__pycache__` xaricində bütün .py fayllar."""
    modules = {}
    for folder, _, files in os.walk(package_dir):
        if "__pycache__" in folder:
            continue
        for name in files:
            if not name.endswith(".py"):
                continue
            path = os.path.relpath(os.path.join(folder, name), ROOT).replace(os.sep, "/")
            module = path[:-3].replace("/", ".")
            if module.endswith(".__init__"):
                module = module[:-len(".__init__")]
            modules[module] = path
    return modules


def _imports(tree: ast.AST, module: str, is_package: bool, known: dict) -> set:
    """Modulun paket daxilində idxal etdiyi modullar."""
    targets = set()
    parts = module.split(".")
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.level:
                base = parts if is_package else parts[:-1]
                base = base[:len(base) - (node.level - 1)]
                target = ".".join(base + ([node.module] if node.module else []))
            else:
                target = node.module or ""
            for candidate in [target] + [f"{target}.{a.name}" for a in node.names]:
                if candidate in known and candidate != module:
                    targets.add(candidate)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name in known and alias.name != module:
                    targets.add(alias.name)
    return targets


def build_graph() -> dict:
    known = _modules(os.path.join(ROOT, PACKAGE))
    nodes, links = [], set()
    for module, path in sorted(known.items()):
        source = open(os.path.join(ROOT, path), encoding="utf-8").read()
        tree = ast.parse(source)
        is_package = path.endswith("__init__.py")
        for target in _imports(tree, module, is_package, known):
            links.add((module, target))
        lines = source.count("\n")
        if is_package and lines < 5:
            continue
        doc = (ast.get_docstring(tree) or "").strip().split("\n")[0][:170]
        parts = module.split(".")
        nodes.append(dict(
            id=module, file=path, layer=parts[1] if len(parts) > 1 else "root",
            loc=lines, doc=doc,
            classes=[n.name for n in tree.body
                     if isinstance(n, ast.ClassDef) and not n.name.startswith("_")][:10],
            funcs=[n.name for n in tree.body
                   if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")][:10]))
    ids = {n["id"] for n in nodes}
    return dict(nodes=nodes,
                links=[dict(source=a, target=b) for a, b in sorted(links)
                       if a in ids and b in ids])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", default="graph.json")
    args = parser.parse_args()
    graph = build_graph()
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(graph, handle, ensure_ascii=False)
    print(f"{len(graph['nodes'])} modul, {len(graph['links'])} idxal əlaqəsi, "
          f"{sum(n['loc'] for n in graph['nodes'])} sətir -> {args.out}")
    print(dict(Counter(n["layer"] for n in graph["nodes"])))


if __name__ == "__main__":
    main()
