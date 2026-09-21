"""İş seansı: son layihələr, dəyişiklik izi, bərpa faylı — Seans 46.

**Qt-dən ASILI DEYİL** — `MainWindow` pytest-də qurula bilmir (VTK
səhnəsi qaçırıcını çökdürür, bax `ui/playback.py`), ona görə qərar
məntiqi burada saxlanılır və ayrıca test olunur. UI yalnız saxlama
yerini (QSettings) və dialoqları idarə edir.
"""

from __future__ import annotations

import json
import os
from typing import Iterable, List, Optional

#: «Son layihələr» menyusunda neçə fayl göstərilir.
RECENT_LIMIT = 5

#: Bərpa faylının adı — hər hesablama bitəndə layihənin SURƏTİ bura
#: yazılır. İstifadəçinin öz `.imx` faylına HEÇ VAXT toxunulmur.
RECOVERY_FILE_NAME = "berpa.imx"

#: Test və xüsusi quraşdırma üçün məlumat qovluğunu dəyişmək imkanı.
DATA_DIR_ENV = "IMEX2D_DATA_DIR"


def _key(path: str) -> str:
    """Eyni faylın iki yazılışı (`C:\\a\\b.imx` / `c:/a/b.imx`) bir sayılır."""
    return os.path.normcase(os.path.abspath(path))


def push_recent(paths: Iterable[str], path: str,
                limit: int = RECENT_LIMIT) -> List[str]:
    """`path`-i siyahının BAŞINA qoyur; təkrarı silir, uzunluğu kəsir."""
    result = [os.path.abspath(path)]
    seen = {_key(path)}
    for item in paths or []:
        if not item or _key(item) in seen:
            continue
        seen.add(_key(item))
        result.append(item)
    return result[:limit]


def remove_recent(paths: Iterable[str], path: str) -> List[str]:
    return [item for item in paths or [] if _key(item) != _key(path)]


def data_dir() -> str:
    """Proqramın öz məlumat qovluğu (bərpa faylı üçün)."""
    override = os.environ.get(DATA_DIR_ENV)
    if override:
        return override
    base = (os.environ.get("LOCALAPPDATA")
            or os.path.join(os.path.expanduser("~"), ".local", "share"))
    return os.path.join(base, "IMEX2D")


def recovery_path() -> str:
    return os.path.join(data_dir(), RECOVERY_FILE_NAME)


def state_signature(ui_state: dict, run_ids: Iterable[str],
                    geology_rows: Iterable[dict]) -> str:
    """İşin «barmaq izi» — saxlanandan sonra dəyişib-dəyişmədiyini bilmək üçün.

    Panellər, işə salınmaların siyahısı və geologiya cədvəli daxildir.
    Bağlayanda iz son saxlanmadakı ilə eyni deyilsə, «Saxlanılsın?» soruşulur.
    """
    return json.dumps({"ui": ui_state, "runs": sorted(run_ids),
                       "geology": list(geology_rows)},
                      sort_keys=True, ensure_ascii=False, default=str)


def is_recovery_file(path: Optional[str]) -> bool:
    return bool(path) and _key(path) == _key(recovery_path())
