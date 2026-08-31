# İş hesabatı — Geologiya bölməsinin quyu cədvəli ilə əvəzlənməsi

Bu fayl işin gedişini, tapşırıqla kod arasındakı fərqləri və öz təşəbbüsü ilə
verilmiş qərarları qeyd edir. Hər mərhələdən sonra yenilənir.

## Mərhələ 1 — Mövcud kodun təhlili

### Əsas tapıntı: `.imx` sxemi tapşırıqdakı təsvirdən fundamental fərqlidir

Tapşırıqda təsvir olunan sxem (bölmə 8) tək-modelli, düz strukturdur:
`{"schema_version": 3, "geology": {...}, "wells": [...]}`.

Faktiki kod (`imex2d/application/serialization.py`, `imex2d/application/project.py`)
**çox-modelli** `Project` kökünə əsaslanır:

```
Project
 ├─ geological_models: Dict[str, GeologicalModel]   (bir neçə ola bilər)
 ├─ reservoir_models:  Dict[str, ReservoirModel]     (bir neçə ola bilər)
 └─ runs:              Dict[str, SimulationRun]
```

Fayl formatı `"version": 1` (tam ədəd) açarı ilə işləyir, `"schema_version"` yox,
və `ProjectSerializer.load()` **dəqiq bərabərliklə** (`!=`) rədd edir — versiyanı
sadəcə 3-ə qaldırmaq bütün mövcud `.imx` fayllarını dərhal sındırardı.

**Qərar:** `Project` kökünə yeni sahələr əlavə edildi (aşağıda), `FORMAT_VERSION`
1 → 2 qaldırıldı, `load()`-dakı yoxlama `version > FORMAT_VERSION` şərtinə
dəyişdirildi (yalnız GƏLƏCƏK versiyanı rədd edir, köhnəni yox). Bu, tapşırığın
"geriyə uyğunluq" tələbini yerinə yetirir, sadəcə hərfi "schema_version": 3
adı ilə deyil.

### `Well`/`Perforation` — koordinat heç vaxt olmayıb

`imex2d/domain/wells.py`-də `Perforation(i, j, k, open, skin)` YALNIZ indekslə
işləyir, heç vaxt X/Y saxlamayıb. Deməli tapşırığın 8-ci bölməsindəki "əgər
`wells` blokunda i/j və ya koordinat varsa, miqrasiya et" şərti əslində baş
vermir (heç vaxt olmayıb) — amma köhnə layihələrdə geologiya CƏDVƏLİ
ÜMUMİYYƏTLƏ yox idi (yalnız CSV və ya sintetik generator var idi). Ona görə
miqrasiya MƏNASI DƏYİŞDİRİLƏRƏK saxlanıldı: köhnə `.imx` açılanda (yeni
`project.geology_wells` açarı yoxdursa) son rezervuar modelinin quyularından
i/j → təxmini X/Y (`(i+0.5)*dx`, `(j+0.5)*dy`) hesablanıb geologiya cədvəlinə
bir dəfə yazılır (bax `ProjectSerializer._migrate_geology_wells`).

### `geology_service.py` toxunulmazdır, amma tapşırıq onun `top`/`bottom`-u
### səth kimi işlətməsini gözləyir — ziddiyyət

Tapşırıq (bölmə 1): "`top`/`bottom`... üst səthin interpolyasiyası üçün
istifadə olunur". Amma bölmə 9: "`geology_service.py` TOXUNULMUR". Faktiki
kodda `WellBasedGeologicalModelBuilder._surface()` yalnız `spec.dip_x/dip_y`-
dən sabit maili səth qurur, quyu nöqtələrindən İNTERPOLYASİYA OLUNMUŞ səthi
qəbul etmir. Bunu düzgün etmək `_surface()`-i dəyişməyi tələb edərdi.

**Qərar:** `geology_service.py` toxunulmadı (bölmə 9 üstünlük təşkil edir).
`top`/`bottom` adi xassə kimi (`TOP`, `BOTTOM` açarları ilə) `WellDataset`-ə
ötürülür və interpolyasiya olunub hesabatda göstərilir, amma nəticə grid-in
`top_depth_map`-ına AVTOMATİK yazılmır. Bu, "6-cı bənd: buraxılmış iş" kimi
qeyd olunur. Gələcəkdə `geology_service.py`-ə `top_depth_map` parametri əlavə
edilə bilər (ayrıca tapşırıq kimi).

### `GeologicalModel` özü quyu saxlamır (qəsdən)

`imex2d/domain/geological_model.py`-nin başlığı: "Nə YOXDUR (qəsdən):
quyular". `GeologicalWell` siyahısını ora əlavə etmək bu invariantı pozardı.

**Qərar:** redaktə olunan `List[GeologicalWell]` cədvəli, interpolyasiya
üsulu/parametrləri və defoltlar `Project` səviyyəsində saxlanılır (`geology_wells`,
`geology_method`, `geology_params`, `geology_defaults`) — `GeologicalModel`
təmiz qalır, yalnız NƏTİCƏ (grid xassələri) daşıyır. `source` ("wells"/
"synthetic") ayrıca saxlanmır — `geology_wells` boşdursa avtomatik "synthetic"
kimi işlənir (təkrarlanan vəziyyəti önləyir).

### `Well` dataclass-a əlavə (qeyri-kəsici)

7-ci bölmə "Perf üst, m / Perf alt, m" metrlə giriş istəyir, amma indeksləri
(k) saxlanılan yeganə şey Peaceman hesablamasına gedir. `Well`-ə iki YENİ,
defolt `None` sahə əlavə edildi: `perf_top: Optional[float]`,
`perf_bottom: Optional[float]`. Bunlar YALNIZ metr girişinin fayl arasında
itməməsi üçündür (indeks grid ölçüsü dəyişəndə köhnəlir, metr qalır);
`perforations` (i,j,k siyahısı) bunlardan HESABLANIR və mühərrikə əvvəlki
kimi gedir — `ReservoirModel`/simulyasiya kodu toxunulmayıb.

### UI: `self.geology_panel` adı sabit qalmalıdır

`tests/test_ui_wiring.py` `self.geology_panel` atributunu və
`self.toolbox.addItem(self.geology_panel, ...)` sətrini axtarır (sinif adını
YOX). Ona görə panel sinfi `WellDataPanel` → `GeologyPanel` adlandırıldı,
amma `main_window.py`-də dəyişən adı `self.geology_panel` olaraq saxlanıldı.

## Mərhələ 2 — `GeologicalWell` + serialization + miqrasiya

`imex2d/domain/geology.py` yaradıldı, `Project`-ə `geology_wells` /
`geology_method` / `geology_params` / `geology_defaults` əlavə olundu
(`GeologicalModel`-ə YOX — bax Mərhələ 1-dəki qərar). `Well`-ə
`perf_top`/`perf_bottom` (defolt `None`) əlavə olundu. `FORMAT_VERSION`
1 → 2, yoxlama `version > FORMAT_VERSION`. Miqrasiya funksiyası köhnə
faylların son rezervuar modelinin quyularından təxmini X/Y qurur.
28 test (`tests/test_geology_wells.py` + mövcud `test_serialization.py`)
keçdi, tam dəst (589 test) də təsirlənmədi.

## Mərhələ 3 — `xy_to_ij` / `depth_to_k`

`imex2d/domain/geometry.py`-ə əlavə edildi. İmza tapşırıqdakı kimi YOX
(`xy_to_ij(x, y, grid)`) — `grid` parametri əvəzinə mövcud `CellGeometry`
işlədildi, çünki bu sinif artıq `grid` + `dx`/`dy`/`dz`/`top_depth`/
`top_depth_map`-ı birləşdirir (məhz "sütun tapıb sonra Z sərhədlərini
oxumaq" üçün lazım olan hər şey). Origin (0,0) qəbul edilir — bu, artıq
`GeologicalGridSpec._cell_centres`-də işlədilən konvensiyadır, ayrıca
origin sahəsi yoxdur. 12 test keçdi.

## Mərhələ 4 — `validate_wells`

`imex2d/domain/geology.py`-ə `ValidationIssue` + `validate_wells()`
əlavə olundu. **Bir qərar tapşırıqdan kənara çıxır:** "Seçilmiş üsul
üçün quyu sayı azdır" bəndi tapşırıqda `error` kimi yazılıb, amma bölmə
3-ün mətni açıq deyir ki, bir xassənin çatışmazlığı DİGƏR xassələrin
hesablanmasını BLOKLAMAMALIDIR ("proqram çökmür, digər xassələr
hesablanır"). Əgər bunu `error` edib "İnterpolyasiya et" düyməsini
söndürsək, bir xassə (məs. Sw) azlıq edən kimi φ/k də hesablana
bilməyəcək — bu, bölmə 3-ün tələbini pozar. Ona görə: xassə-üzrə
çatışmazlıq `warning` (bloklamır, yalnız məlumatlandırır); YALNIZ heç
bir xassə üçün kifayət qədər quyu yoxdursa (interpolyasiyanın nəticəsi
tam boş olacaqsa) `error`. Bu, `wells_to_dataset` adapteri (Mərhələ 5)
ilə eyni məntiqi paylaşır: yetərsiz xassə sadəcə göndərilən dəstdən
çıxarılır, `geology_service.py` onu heç görmür.

İkinci kiçik fərq: "IDW ≥ 1" tapşırıqda yazılıb, amma paylaşılan
`WellDataset.validate()` (CSV importunda da işlədilir) hər xassə üçün
minimum 2 nöqtə tələb edir (`test_single_point_dataset_is_rejected`).
Bunu boşaltmaq CSV importunun mövcud sınağını sındırardı, ona görə
faktiki minimum bütün üsullar üçün 2-yə qaldırıldı (`method_minimum()`,
bax kodun içindəki şərh). 17 test yazıldı, hər yoxlama üçün ən azı biri.

## Mərhələ planı (bu fayl hər mərhələdən sonra yenilənəcək)

- [x] 1 — təhlil
- [x] 2 — `GeologicalWell` + serialization + miqrasiya + testlər
- [x] 3 — `xy_to_ij` / `depth_to_k` + testlər
- [x] 4 — `validate_wells` + testlər
- [ ] 5 — 2-ci bölmənin UI-si (+ `wells_to_dataset` adapteri)
- [ ] 6 — 7-ci bölmə ilə birləşmə + ssenari generatoru
- [ ] 7 — yekun: golden + tam test dəsti + `run.bat`
