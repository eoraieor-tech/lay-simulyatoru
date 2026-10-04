"""Tutulmamış istisnalar üçün qlobal tutucu — Seans 53 (X-1).

NİYƏ. PyQt5-də Qt slotunun içində tutulmamış Python istisnası qalxanda,
`sys.excepthook` standart olduqda PyQt `qFatal()` çağırır və proses
dərhal dayanır (Windows-da çıxış kodu `0xC0000409`). Seans 52-də ölçüldü:
istifadəçi heç bir mesaj görmürdü, `logs/imex2d.log`-a heç nə yazılmırdı,
iz yalnız stderr-də qalırdı. `sys.excepthook` dəyişdirildikdə PyQt onu
çağırır və prosesi DAYANDIRMIR.

Bu modul Qt-dən asılı deyil (testlər `MainWindow` qurmadan yoxlayır):
mesajı göstərən funksiya `notify` kimi kənardan verilir (`app.py`).
"""

from __future__ import annotations

import logging
import sys
import traceback
from typing import Callable, Optional

#: istifadəçiyə göstərilən mətnin sonu — iz harada axtarılsın
LOG_HINT = "Ətraflı iz jurnal tabında və logs/imex2d.log faylındadır."


def user_message(exc_type, exc) -> str:
    """İstifadəçi üçün qısa mətn: xətanın növü və mətni + log ipucu."""
    detail = "".join(traceback.format_exception_only(exc_type, exc)).strip()
    return ("Gözlənilməz xəta baş verdi, əməliyyat yarımçıq qaldı. "
            "Proqram işləməyə davam edir — işinizi saxlamağınız tövsiyə "
            f"olunur.\n\n{detail}\n\n{LOG_HINT}")


def make_hook(logger: logging.Logger,
              notify: Optional[Callable[[str], None]] = None):
    """`sys.excepthook` imzalı funksiya qaytarır.

    * `KeyboardInterrupt` (Ctrl+C) olduğu kimi standart tutucuya ötürülür;
    * qalanları TAM izi ilə loga yazılır (`CRITICAL`) və `notify`-a
      qısa mətn verilir;
    * `notify`-ın özü xəta versə və ya mesaj göstərilərkən yeni istisna
      qalxsa (təkrar giriş), yalnız loga yazılır — sonsuz dövrə olmur.
    """
    state = {"busy": False}

    def hook(exc_type, exc, tb):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc, tb)
            return
        logger.critical("Tutulmamış istisna", exc_info=(exc_type, exc, tb))
        if notify is None or state["busy"]:
            return
        state["busy"] = True
        try:
            notify(user_message(exc_type, exc))
        except Exception:
            logger.exception("Xəta mesajı göstərilə bilmədi")
        finally:
            state["busy"] = False

    return hook


def install(logger: logging.Logger,
            notify: Optional[Callable[[str], None]] = None):
    """Tutucunu `sys.excepthook`-a qoyur; əvvəlkini qaytarır."""
    previous = sys.excepthook
    sys.excepthook = make_hook(logger, notify)
    return previous
