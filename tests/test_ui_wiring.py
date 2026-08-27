"""İnterfeys bağlantılarının statik yoxlanışı — PyQt5 tələb etmir.

NİYƏ LAZIMDIR: A1-də `update_pvt_plot()` metodu yazıldı, amma `pvt_axes`
widget-i yaradılmadı. 57 testin hamısı keçdi, proqram isə açılan kimi
AttributeError verdi — çünki heç bir test interfeysə toxunmurdu.

Bu test main_window.py-ni AST kimi oxuyur və oxunan hər `self.X`
atributunun harasa mənimsədildiyini yoxlayır. Qt qurulmadan işləyir.
"""

import ast
import os

UI_FILE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "imex2d", "ui", "main_window.py")

# QMainWindow-dan miras qalan üzvlər (bu faylda mənimsədilmir)
INHERITED = {"tabs", "close", "statusBar", "menuBar", "show"}


def _class_node(tree, name):
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == name:
            return node
    raise AssertionError(f"{name} sinfi tapılmadı")


def _collect(class_node):
    assigned, read, called = set(), set(), set()
    for node in ast.walk(class_node):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) \
                and node.value.id == "self":
            if isinstance(node.ctx, ast.Store):
                assigned.add(node.attr)
            else:
                read.add(node.attr)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "self":
            called.add(node.func.attr)
    methods = {n.name for n in class_node.body if isinstance(n, ast.FunctionDef)}
    return assigned, read, called, methods


def test_every_read_attribute_is_assigned_somewhere():
    tree = ast.parse(open(UI_FILE, encoding="utf-8").read())
    assigned, read, called, methods = _collect(_class_node(tree, "MainWindow"))
    missing = read - assigned - called - methods - INHERITED
    assert not missing, ("Bu atributlar oxunur, amma heç yerdə yaradılmır: "
                         + ", ".join(sorted(missing)))


def test_expected_tabs_are_created():
    """Hər tab yaradılmalıdır — sayı və adları qorunur."""
    source = open(UI_FILE, encoding="utf-8").read()
    for title in ("Layihə", "Model", "Nəticələr", "Nisbi keçiricilik",
                  "PVT", "Validasiya (B-L)", "Jurnal"):
        assert f'"{title}"' in source, f"'{title}' tabı yaradılmır"


def test_panel_widgets_are_instantiated():
    source = open(UI_FILE, encoding="utf-8").read()
    for panel in ("GridGeometryPanel", "RockFluidPanel", "ScalPanel",
                  "PvtPanel", "WellPanel", "NumericalPanel"):
        assert f"{panel}()" in source, f"{panel} yaradılmır"


def test_every_panel_reaches_the_toolbox():
    """Hər panel ya birbaşa, ya da ara qutu vasitəsilə sol panelə düşür.

    `scal_panel` və `scal_source_panel` bir səhifədə birləşdirilib
    (SCAL mənbəyi + Corey parametrləri), ona görə birbaşa `addItem`
    yoxdur — yoxlama layout-a əlavə olunmanı da qəbul edir.
    """
    source = open(UI_FILE, encoding="utf-8").read()
    for attribute in ("grid_panel", "geology_panel", "rock_panel",
                      "scal_panel", "scal_source_panel", "pvt_panel",
                      "well_panel", "numerical_panel"):
        direct = f"self.toolbox.addItem(self.{attribute}" in source
        wrapped = f"addWidget(self.{attribute}" in source
        assert direct or wrapped, f"{attribute} sol panelə əlavə edilmir"


def test_toolbox_sections_are_numbered_in_order():
    """Bölmə nömrələri istifadəçi üçün naviqasiyadır — sıra pozulmamalıdır."""
    import re

    source = open(UI_FILE, encoding="utf-8").read()
    titles = re.findall(r'self\.toolbox\.addItem\([^,]+,\s*"(\d+)\s*·',
                        source)
    numbers = [int(value) for value in titles]
    assert numbers == sorted(numbers), numbers
    assert numbers == list(range(1, len(numbers) + 1)), numbers


def test_renderers_are_instantiated():
    source = open(UI_FILE, encoding="utf-8").read()
    for renderer in ("MapRenderer", "ProductionCurveRenderer", "ScalRenderer",
                     "PvtRenderer", "ValidationRenderer"):
        assert f"R.{renderer}()" in source, f"{renderer} yaradılmır"


# ── grid panelinin qalınlıq rejimi ────────────────────────────────────
def test_layer_thickness_modes_are_consistent():
    """Baza dərinliyi ilə DZ eyni kəmiyyəti təyin edir.

    Üç kəmiyyətdən (tavan, qalınlıq, baza) yalnız ikisi müstəqildir.
    Panel ikisini qəbul edir, üçüncünü hesablayır — istifadəçi hansı
    cütü verəcəyini seçir.
    """
    import ast
    import os

    source = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "imex2d", "ui", "panels.py")
    tree = ast.parse(open(source, encoding="utf-8").read())

    panel = next(node for node in ast.walk(tree)
                 if isinstance(node, ast.ClassDef)
                 and node.name == "GridGeometryPanel")
    methods = {item.name for item in panel.body
               if isinstance(item, ast.FunctionDef)}
    assert "layer_thickness" in methods
    assert "_on_mode_changed" in methods

    # `values()` xam `dz`-i deyil, hesablanmış qalınlığı verməlidir
    values = next(item for item in panel.body
                  if isinstance(item, ast.FunctionDef) and item.name == "values")
    text = ast.dump(values)
    assert "layer_thickness" in text, "values() xam dz işlədir"
