"""UI qatının statik yoxlanışı — PyQt quraşdırılmadan işləyir.

Səbəb: UI kodu avtomatik testlərlə əhatə olunmurdu və iki dəfə eyni
tipli səhv buraxıldı — istifadə olunan, lakin import edilməyən ad
(`QSpinBox`) və yaradılmadan çağırılan atribut (`pvt_axes`). Bu testlər
həmin sinif səhvləri tutur: Qt işə salınmır, yalnız AST təhlili aparılır.
"""

import ast
import os
import re

UI_DIRECTORY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "imex2d", "ui")

QT_NAME = re.compile(r"^Q[A-Z]")

# Qt baza siniflərindən (QWidget, QMainWindow, QThread) gələn üzvlər —
# metod kimi çağırılmadan ötürülə bilər (məsələn siqnala qoşulanda).
QT_INHERITED_MEMBERS = {
    "close", "show", "hide", "update", "repaint", "raise_", "deleteLater",
    "start", "quit", "exit", "accept", "reject", "isVisible", "parent",
}


def _modules():
    return [os.path.join(UI_DIRECTORY, name)
            for name in sorted(os.listdir(UI_DIRECTORY))
            if name.endswith(".py")]


def _tree(path):
    with open(path, encoding="utf-8") as handle:
        return ast.parse(handle.read(), filename=path)


def _bound_names(tree):
    """Modulda mövcud olan bütün adlar: import, təyinat, funksiya, sinif."""
    names = set(dir(__builtins__)) if isinstance(__builtins__, dict) else set(dir(__builtins__))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            names.update(alias.asname or alias.name for alias in node.names)
        elif isinstance(node, ast.Import):
            names.update((alias.asname or alias.name).split(".")[0]
                         for alias in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            names.add(node.id)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
    return names


def test_every_qt_name_used_in_ui_is_imported():
    problems = []
    for path in _modules():
        tree = _tree(path)
        available = _bound_names(tree)
        used = {node.id for node in ast.walk(tree)
                if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
                and QT_NAME.match(node.id)}
        missing = sorted(used - available)
        if missing:
            problems.append(f"{os.path.basename(path)}: {', '.join(missing)}")
    assert not problems, "Import edilməyən Qt adları: " + "; ".join(problems)


def test_self_attributes_used_in_ui_are_assigned_somewhere():
    """`self.x` oxunursa, həmin sinifdə `self.x = ...` təyinatı olmalıdır."""
    problems = []
    for path in _modules():
        for node in ast.walk(_tree(path)):
            if not isinstance(node, ast.ClassDef):
                continue
            assigned, used = set(), set()
            # sinif səviyyəsində təyin olunan adlar (siqnallar, sabitlər)
            for statement in node.body:
                if isinstance(statement, ast.Assign):
                    assigned.update(t.id for t in statement.targets
                                    if isinstance(t, ast.Name))
                elif isinstance(statement, ast.AnnAssign) and isinstance(
                        statement.target, ast.Name):
                    assigned.add(statement.target.id)

            # `self.x(...)` çağırışları — baza sinifdən (QWidget) gələ bilər,
            # ona görə yalnız MƏLUMAT kimi oxunan atributlar yoxlanılır
            called = {inner.func.attr for inner in ast.walk(node)
                      if isinstance(inner, ast.Call)
                      and isinstance(inner.func, ast.Attribute)
                      and isinstance(inner.func.value, ast.Name)
                      and inner.func.value.id == "self"}

            for inner in ast.walk(node):
                if not (isinstance(inner, ast.Attribute)
                        and isinstance(inner.value, ast.Name)
                        and inner.value.id == "self"):
                    continue
                (assigned if isinstance(inner.ctx, ast.Store) else used).add(inner.attr)
            methods = {m.name for m in node.body
                       if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef))}
            missing = sorted(used - assigned - methods - called
                             - QT_INHERITED_MEMBERS)
            if missing:
                problems.append(f"{os.path.basename(path)}.{node.name}: "
                                f"{', '.join(missing)}")
    assert not problems, ("Təyin edilməmiş atributlar: " + "; ".join(problems))


def test_builder_methods_are_actually_called():
    """`_build_*` metodu yazılıb, amma çağırılmayıbsa — tab görünməz.

    Məhz bu baş verdi: `_build_volume_tab()` yaradıldı, lakin
    `addTab` sətri tətbiq olunmadı. Metod daxilindəki `self.volume_time`
    təyinatı mövcud olduğu üçün atribut testi bunu TUTMADI — həmin
    kod heç vaxt icra olunmurdu.
    """
    problems = []
    for path in _modules():
        tree = _tree(path)
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            builders = {
                method.name for method in node.body
                if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef))
                and method.name.startswith("_build")
            }
            called = {
                inner.func.attr for inner in ast.walk(node)
                if isinstance(inner, ast.Call)
                and isinstance(inner.func, ast.Attribute)
                and isinstance(inner.func.value, ast.Name)
                and inner.func.value.id == "self"
            }
            missing = sorted(builders - called)
            if missing:
                problems.append(f"{os.path.basename(path)}.{node.name}: "
                                f"{', '.join(missing)}")
    assert not problems, ("Çağırılmayan qurucu metodlar: " + "; ".join(problems))


def test_ui_modules_compile():
    import py_compile
    for path in _modules():
        py_compile.compile(path, doraise=True)


def test_tabs_are_selected_by_name_not_index():
    """Sərt tab indeksləri qadağandır.

    3D tab əlavə olunanda bütün indekslər sürüşdü və proqram səhv
    səhifəyə keçdi. `show_tab("Ad")` bu problemi aradan qaldırır.
    """
    problems = []
    for path in _modules():
        source = open(path, encoding="utf-8").read()
        inside_helper = False
        for number, line in enumerate(source.split("\n"), start=1):
            stripped = line.strip()
            if stripped.startswith("def "):
                inside_helper = stripped.startswith("def show_tab")
            if inside_helper:
                continue          # köməkçinin öz daxili işi
            if stripped.startswith("self.tabs.setCurrentIndex("):
                problems.append(f"{os.path.basename(path)}:{number}")
    assert not problems, ("Sərt tab indeksi işlədilib (show_tab işlət): "
                          + ", ".join(problems))


def test_expected_tabs_match_the_tabs_the_window_registers():
    """`version.EXPECTED_TABS` faktiki `addTab` çağırışları ilə uyğun olmalıdır.

    Bu siyahı işə salınanda yoxlanış üçün işlədilir; kod ilə ayrılsa,
    yoxlama yalançı xəbərdarlıq verər.
    """
    import re
    import sys

    sys.path.insert(0, os.path.dirname(UI_DIRECTORY))
    from imex2d.version import EXPECTED_TABS

    source = open(os.path.join(UI_DIRECTORY, "main_window.py"),
                  encoding="utf-8").read()
    registered = re.findall(r'addTab\([^,]+,\s*"([^"]+)"\)', source)
    assert registered == EXPECTED_TABS, \
        f"kod: {registered}\nversion.py: {EXPECTED_TABS}"


def _self_calls(node):
    """`self.metod()` çağırışlarının adları."""
    return {inner.func.attr for inner in ast.walk(node)
            if isinstance(inner, ast.Call)
            and isinstance(inner.func, ast.Attribute)
            and isinstance(inner.func.value, ast.Name)
            and inner.func.value.id == "self"}


def _self_reads(node):
    """Oxunan `self.atribut` adları (çağırışlar istisna)."""
    called = _self_calls(node)
    return {inner.attr for inner in ast.walk(node)
            if isinstance(inner, ast.Attribute)
            and isinstance(inner.value, ast.Name)
            and inner.value.id == "self"
            and isinstance(inner.ctx, ast.Load)
            and inner.attr not in called}


def _reachable_reads(method_name, methods, seen=None):
    """Metodun ÖZÜ və çağırdığı bütün metodların oxuduğu atributlar.

    Zənciri izləmək vacibdir: `_on_mode_changed()` özü `self.info`-nu
    oxumur, lakin `_refresh_info()`-nu çağırır, o isə oxuyur.
    """
    seen = seen or set()
    if method_name in seen or method_name not in methods:
        return set()
    seen.add(method_name)
    target = methods[method_name]
    reads = set(_self_reads(target))
    for nested in _self_calls(target):
        reads |= _reachable_reads(nested, methods, seen)
    return reads


def test_init_does_not_call_methods_before_their_attributes_exist():
    """`__init__` daxilində metod çağırışı erkən ola bilər.

    Nümunə: `GridGeometryPanel.__init__` `_on_mode_changed()`-i
    çağırırdı; o özü `self.info`-nu oxumurdu, lakin `_refresh_info()`-nu
    çağırırdı, o isə oxuyurdu — və `self.info` hələ yaradılmamışdı.
    Proqram açılışda çökürdü.

    Ona görə çağırış ZƏNCİRİ izlənir, yalnız birinci səviyyə yox.
    """
    problems = []
    for path in _modules():
        tree = _tree(path)
        for klass in ast.walk(tree):
            if not isinstance(klass, ast.ClassDef):
                continue
            methods = {item.name: item for item in klass.body
                       if isinstance(item, ast.FunctionDef)}
            init = methods.get("__init__")
            if init is None:
                continue

            assigned_in_init = {
                node.attr for node in ast.walk(init)
                if isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "self"
                and isinstance(node.ctx, ast.Store)
            }

            defined_so_far = set()
            for statement in init.body:
                for call in _self_calls(statement):
                    required = _reachable_reads(call, methods)
                    missing = sorted((required & assigned_in_init)
                                     - defined_so_far - set(methods))
                    if missing:
                        problems.append(
                            f"{os.path.basename(path)}.{klass.name}: "
                            f"{call}() -> {', '.join(missing)}")
                for node in ast.walk(statement):
                    if (isinstance(node, ast.Attribute)
                            and isinstance(node.value, ast.Name)
                            and node.value.id == "self"
                            and isinstance(node.ctx, ast.Store)):
                        defined_so_far.add(node.attr)

    assert not problems, ("__init__-də erkən çağırış: " + "; ".join(problems))


def test_three_dimensional_tab_has_no_navigation_toolbar():
    """3D-də matplotlib alət paneli fırlatmanı bloklayır.

    `Axes3D._on_move` başlanğıcda yoxlayır:
        get_navigate_mode() is not None -> return
    Yəni panelin "pan" və ya "zoom" rejimi aktiv olan kimi siçanla
    fırlatmaq tamamilə dayanır. 2D xəritədə panel faydalıdır, 3D-də yox.
    """
    source = open(os.path.join(UI_DIRECTORY, "main_window.py"),
                  encoding="utf-8").read()
    builder = source.split("def _build_volume_tab")[1].split("\n    def ")[0]
    assert "NavToolbar" not in builder, \
        "3D tabına naviqasiya paneli qaytarılıb — fırlatma bloklanacaq"


def test_scroll_wheel_zoom_is_connected():
    source = open(os.path.join(UI_DIRECTORY, "main_window.py"),
                  encoding="utf-8").read()
    assert 'mpl_connect("scroll_event"' in source
    assert "_on_volume_scroll" in source
