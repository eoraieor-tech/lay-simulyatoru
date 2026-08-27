"""Versiya məlumatı.

Səbəb: eyni adlı fayllar müxtəlif qovluqlarda qarışa bilir və hansı
buraxılışın işlədiyini görmək mümkün olmurdu. İndi versiya pəncərənin
başlığında, jurnal faylında və "Kömək → Versiya" pəncərəsindədir.

Yeni buraxılışda YALNIZ bu fayl dəyişir.
"""

from __future__ import annotations

VERSION = "67"
RELEASE_DATE = "2026-08-21"

FEATURES = [
    ("Qatlı arxitektura", "domain / interfaces / application / simulation / rendering / ui"),
    ("A1 · PVT", "təzyiqdən asılı Bo, μo, Rs, Bw; korrelyasiyalar"),
    ("A3 · Equilibration", "hidrostatik təzyiq, OWC, maili lay"),
    ("A4 · Kapilyar + cazibə", "Brooks-Corey Pc, keçid zonası, faza potensialı"),
    ("A5 · 3D grid", "nz > 1, şaquli axın, Kv/Kh, kəsik görüntüsü"),
    ("A6 · Fully implicit", "Nyuton, analitik Jakobian, adaptiv Δt, CPR"),
    ("B1 · Layihə faylı", ".imx saxlama/açma, işə salınmaların müqayisəsi"),
    ("B2 · Geoloji import", "quyu CSV, IDW / Kriging / ən yaxın qonşu"),
    ("3D həcm görüntüsü", "siçanla fırlatma/sürüşdürmə/yaxınlaşdırma, kəsim həddi"),
    ("B4 · Region əsaslı SCAL", "laboratoriya cədvəlləri, SATNUM regionları"),
    ("B5 · Eclipse mübadiləsi", "GRDECL oxuma, .DATA yazma, SWOF importu"),
    ("C6 · Həssaslıq analizi", "Tornado diaqramı, yerli elastiklik"),
    ("C5 · Tarixçə uyğunluğu", "müşahidə CSV, NRMSE ölçüsü, çarpaz qrafik"),
    ("B3 · Fault transmissivliyi", "I/J/K müstəvisi, çarpan, sealing, CSV/Eclipse FAULTS, 3D-də görünən"),
    ("B6 · PDF hesabat", "çoxsəhifəli hesabat: model, xəritə, SCAL, PVT, nəticələr"),
    ("A7 · Qaz fazası (sınaq statusunda)", "PVT tabında aktivləşdirilir — GOC, qaz-neft SCAL, təhlükəsiz uğursuzluqla"),
    ("OPM Flow idxalı", "OPM Flow nəticələrini öz 3D görüntümüzlə göstər — .EGRID+.UNRST"),
    ("VTK 3D motoru", "ResInsight tipli — quyular, faultlar, ölçü oxları, istiqamət oxu"),
    ("Diaqnostika", "xəta/xəbərdarlıq səviyyələri, quyu rejimi yoxlamaları"),
]

EXPECTED_TABS = [
    "Layihə", "Model", "Nəticələr", "Nisbi keçiricilik", "3D görüntü",
    "PVT", "Validasiya (B-L)", "Müqayisə", "Tarixçə", "Uyğunlaşdırma",
    "Həssaslıq", "Jurnal",
]
"""Bu buraxılışda olmalı olan tablar.

Proqram açılanda faktiki tablarla müqayisə olunur — köhnə fayl
qarışıbsa, jurnalda dərhal xəbərdarlıq görünür.
"""


def title() -> str:
    return f"IMEX-2D v{VERSION}  ·  Geoloji modelləşdirmə və rezervuar simulyasiyası"


def summary() -> str:
    lines = [f"IMEX-2D versiya {VERSION}   ({RELEASE_DATE})", ""]
    for name, description in FEATURES:
        lines.append(f"  {name}")
        lines.append(f"      {description}")
    return "\n".join(lines)
