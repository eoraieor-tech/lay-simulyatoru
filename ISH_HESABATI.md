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

## Mərhələ 5 — 2-ci bölmənin UI-si

`imex2d/application/geology_adapter.py` (`wells_to_dataset`) + yeni
`imex2d/ui/geology_map.py` (`GeologyMapWidget`, sadə `QPainter` xəritəsi:
grid düzbucaqlısı + quyu nöqtələri, rəng modeldə/məlumat-yalnız/kənar üzrə)
+ `imex2d/ui/panels.py`-də `WellDataPanel` → `GeologyPanel` (sinif adı
dəyişdi, `self.geology_panel` atribut adı EYNİ qaldı — `test_ui_wiring.py`
bunu tələb edir). Cədvəl: `Ad | Modeldə | X,m | Y,m | (i,j) | Lay üstü,m |
Lay altı,m | φ | k,mD | Sw | Qeyd`. Alət paneli (Quyu əlavə et / Dublikat /
Sil / Grid mərkəzinə at) + yoxlama paneli (`validate_wells` nəticəsi) +
"İnterpolyasiya et" düyməsi + "Nəticə köhnəlib" göstəricisi.

`main_window.py`-də `geology_panel.changed` BİLƏRƏKDƏN `rebuild_model`-ə
qoşulmayıb (yalnız `interpolate_requested` → `_interpolate_geology`) —
cədvəl redaktəsi böyük gridi avtomatik yenidən interpolyasiya etmir.
Hesablanan model `_geology_model_from_wells`-də keşlənir; cədvəl boşdursa
və ya hələ interpolyasiya edilməyibsə `_build_geological_model()` sintetik
modelə düşür (izahlı mesajla, çökmür).

Pəncərə başlığının `*`-si `geology_panel.changed` → `_mark_dirty()` ilə
idarə olunur (yalnız 2-ci bölmə redaktəsi, tapşırıqda göstərildiyi kimi).

**Diqqətəlayiq addım:** startup zamanı (`__init__`) əvvəllər YALNIZ 7-ci
bölmə default five-spot alırdı, geologiya cədvəli boş qalırdı. Bu artıq
mümkün DEYİL — 7-ci bölmənin i/j/k-sı indi TAM geologiya bağlantısından
hesablanır, boş geologiya cədvəli ilə hər iki default quyu (0,0) hüceyrəsinə
düşərdi (bax aşağı, Mərhələ 6). Ona görə startup-da default five-spot HƏM
`well_panel`-ə, HƏM `geology_panel`-ə (yalnız X/Y, petrofizikasız) yazılır.
`geology_source` nəticədə startup-da "wells" olur, AMMA `_geology_model_from_wells`
hələ `None` olduğu üçün grid xassələri yenə SİNTETİKDİR — status mesajı
bunu düzgün əks etdirir. Defolt davranış (sintetik xassələr) qorunur,
yalnız daxili mexanizm dəyişib.

## Mərhələ 6 — 7-ci bölmə ilə birləşmə + ssenari generatoru

`WellPanel` tam yenidən yazıldı: `Ad | i | j | Perf üst,m | Perf alt,m | k
| Tip | İdarə | Qiymət | rw`. `Ad`/`i`/`j`/`k` REDAKTƏSİZ (boz) — `Ad`
geologiya adına bağlıdır, `i`/`j` `xy_to_ij`, `k` `depth_to_k` ilə
hesablanır. `set_geology_context(wells, geometry)` sətirləri `in_model`
dəstinə görə avtomatik əlavə/silir; silinən sətirin rejimi `_retained`
lüğətində qalır — yenidən işarələnəndə geri qayıdır (tapşırığın tələbi).

`clamp_to_grid`/`set_layer_count`/`add_button`/`remove_button` silindi —
mənasız qaldı, çünki indeks artıq HEÇ VAXT əl ilə redaktə olunmur (yalnız
metrdən hesablanır), grid ölçüsü dəyişəndə köhnəlmə problemi öz-özünə
aradan qalxır.

Ssenari generatoru (`_apply_pattern`): daxili `five_spot()` və s. İNDEKSLƏ
işləməyə davam edir (dəyişmədi). Tətbiq ediləndə nəticə metrə çevrilib
(`_wells_to_geology_rows`) HƏM `well_panel`, HƏM `geology_panel`-ə yazılır.
Mövcud geologiya sətirləri varsa əvvəlcə təsdiq soruşulur (`QMessageBox`).

Perforasiya lay qalınlığından kənardadırsa (`depth_to_k` → `None`)
xəbərdarlıq `WellPanel.warning_label`-də göstərilir (bloklamır).

**Uçdan-uca sınaq** (`QT_QPA_PLATFORM` OFFSCREEN YOX — VTK offscreen-də
segfault verir, real Windows platforması ilə skript sınandı): startup →
default five-spot iki bölmədə də düzgün göründü, INJ-1/PROD-1 küncləri
(0,0)/(40,40) dəqiq tutdu; petrofizika əlavəsi + 3-cü quyu + interpolyasiya
işlədi; perforasiya metr redaktəsi k-nı yenilədi; `in_model` söndürülüb-
yandırılanda rejim qorundu; boş ad + kənar koordinat "İnterpolyasiya et"-i
blokladı; `.imx` saxla/aç dövrəsi bütün quyuları və `geology_source`-u
qorudu.

**Kənara çıxma (commit strukturu):** `GeologyPanel` və `WellPanel` eyni
faylda (`panels.py`) sıx bağlı şəkildə yazıldı və birlikdə sınandı — ona
görə Mərhələ 5 və 6 AYRI-AYRI commit əvəzinə BİR birləşmiş UI commiti kimi
yazılır ("ən mühafizəkar" seçim: yarımçıq/qeyri-ardıcıl aralıq vəziyyəti
commit etməkdənsə, tam sınanmış vahid dəyişiklik).

## Mərhələ 7 — yekun

**Tam test dəsti:** `pytest -q` → **624 keçdi, 1 keçildi (skip)**, 0 uğursuz.
Skip olunan test bu işdən əvvəl də mövcud idi və bu işə aid deyil (yavaş
test, `IMEX_SKIP_SLOW` ilə əlaqəli — toxunulmayıb). Bu işlə əlaqəli yeni
testlər: `test_geology_wells.py` (14), `test_geology_geometry.py` (12),
`test_validate_wells.py` (17), `test_geology_adapter.py` (6) — cəmi 49.

**Golden:** `python tools/golden.py` (yazma REJİMİ İŞLƏDİLMƏDİ) → üç keys
də (`five_spot`, `bl_1d`, `five_spot_small`) **UYĞUNDUR**. Gözlənilən
nəticə — bu iş yalnız daxiletmə yolunu dəyişir, `geology_service.py` və
simulyasiya mühərriki toxunulmayıb.

**`run.bat` sınağı — qismən, alət mühiti məhdudiyyəti ilə:** Bu mühitdə
(sandboxlanmış alət icrası) uzunmüddətli GUI prosesini bir neçə alət
çağırışı arasında canlı saxlamaq mümkün olmadı (`cmd.exe /c run.bat`
Bash alətində sükutla heç nə icra etmir — köhnə davranış, mənim işimlə
əlaqəli deyil; PowerShell `Start-Process` isə prosesi düzgün başladır,
amma alət çağırışları arasında iş meneceri onu dayandırır). Ona görə
əvəzinə üç səviyyəli yoxlama aparıldı:
  1. `venv\Scripts\python.exe app.py` (run.bat-ın aktivləşdirdiyi EYNİ venv)
     — İKİ DƏFƏ təmiz başladı, jurnalda "12 tab" doğrulaması keçdi, xəta yoxdur.
  2. PowerShell `Start-Process` ilə `run.bat`-ın ÖZÜ işə salındı — jurnalda
     "IMEX-2D v69 başladıldı" sətri göründü (yəni `cd`/venv aktivləşdirmə/
     `python app.py` zənciri düzgün işləyir), sonra alət mühiti prosesi
     kəsdi (`^C`) — bu, `run.bat`-ın özündə DEYİL, mənim test roll-umda idi.
  3. Real Qt platforması ilə (OFFSCREEN YOX) ssenari üzrə skriptli sınaq:
     `MainWindow` qurulması, default five-spot-un HƏM 2-ci, HƏM 7-ci
     bölmədə düzgün göründüyü (INJ-1/PROD-1 küncləri (0,0)/(40,40) dəqiq),
     petrofizika daxil edilib interpolyasiya, perforasiya metr redaktəsi,
     `in_model` söndürülüb-yandırılanda rejimin qorunması, boş ad/kənar
     koordinatın "İnterpolyasiya et"-i bloklaması, `.imx` saxla/aç dövrəsi.
     Hamısı gözlənilən nəticəni verdi (bax Mərhələ 6).
  4. `QT_QPA_PLATFORM=offscreen` ilə TAM `MainWindow` qurulması SEGFAULT
     verir — səbəb VTK-nin (3D görüntü tabı) offscreen rejimdə OpenGL
     kontekst ala bilməməsidir. Bu, MƏNİM DƏYİŞİKLİYİMLƏ ƏLAQƏLİ DEYİL —
     VTK-nin real ekran/GPU tələb etməsi əvvəldən mövcud məhdudiyyətdir.

  **Tövsiyə:** qayıdanda `run.bat`-ı əl ilə açıb ən azı bunlara bax: (a)
  2-ci bölmənin kiçik xəritəsi düzgün göstərilir, (b) "Grid mərkəzinə at"
  düyməsi, (c) 7-ci bölmədə Perf üst/alt redaktəsi ilə k sütununun canlı
  yenilənməsi. Bunlar VİZUAL layout məsələləridir, skriptlə tam yoxlanıla
  bilmədi (yalnız MƏNTİQ yoxlanıldı, PİKSEL-səviyyəli görünüş yox).

### Bilinən məhdudiyyətlər / sonraya buraxılan işlər

- **`top`/`bottom` grid səthinə köçmür:** quyu `top`/`bottom` dəyərləri
  `TOP`/`BOTTOM` adı ilə interpolyasiya olunur və hesabatda göstərilir,
  AMMA nəticə `CellGeometry.top_depth_map`-a avtomatik yazılmır — bunun
  üçün `geology_service.py`-ə toxunmaq lazım olardı, bu isə bölmə 9-un
  qadağasına ziddir (bax Mərhələ 1). Gələcək iş: `_surface()`-ə real
  interpolyasiya olunmuş səthi qəbul etmək imkanı əlavə etmək.
- **Xəritə statikdir:** `GeologyMapWidget` nöqtəni klikləyib sətri seçmə
  imkanı vermir (yalnız cədvəldə seçilən sətir xəritədə fərqli rənglə
  göstərilir, əksi yox). Vaxt məhdudiyyətinə görə minimal saxlanıldı.
- **CSV idxalı** (`imex2d/geology/well_data_io.py`) toxunulmadan qalıb,
  UI-dən ayrılıb — tapşırığın 13-cü bəndinə uyğun olaraq bilərəkdən indi
  geri qoşulmayıb.
- **PyQt widget-lərinin avtomatlaşdırılmış testi yoxdur** — layihənin
  mövcud konvensiyası (`test_ui_wiring.py`/`test_ui_static.py`) Qt
  ekranı tələb etməyən AST-əsaslı statik yoxlamalardır; mən də bu
  qaydaya uyğunlaşdım. `GeologyPanel`/`WellPanel`-in davranışı yuxarıdakı
  skriptli sınaqla (commit mesajında təsvir olunan) doğrulandı, amma bu
  sınaq test dəstinə ƏLAVƏ OLUNMAYIB (təkrarlana bilən avtomatik test
  deyil, real Qt platforması və müvəqqəti fayl tələb edir).

### Öz təşəbbüsümlə verilmiş əsas qərarlar (xülasə)

Hamısı yuxarıda müvafiq mərhələdə ətraflı izah olunub, burada siyahı kimi:
1. `.imx` sxemi tapşırıqdakı kimi deyil, mövcud `Project`-based struktura
   uyğunlaşdırıldı (Mərhələ 1).
2. `FORMAT_VERSION` "3" yox, `1 → 2`; yoxlama `>` (Mərhələ 1/2).
3. Geologiya quyu cədvəli `GeologicalModel`-də deyil, `Project`-də saxlanılır
   (Mərhələ 1).
4. `xy_to_ij`/`depth_to_k` imzası `grid` əvəzinə mövcud `CellGeometry`
   qəbul edir (Mərhələ 3).
5. Metod-üzrə minimum quyu sayı "IDW ≥ 1" yox, faktiki "IDW ≥ 2" (paylaşılan
   `WellDataset.validate()`-in tələbinə görə) (Mərhələ 4).
6. Xassə-üzrə az-quyu yoxlaması `error` yox, `warning` (bölmə 3-ün "digər
   xassələr hesablanır" tələbini qorumaq üçün) (Mərhələ 4).
7. `top`/`bottom` interpolyasiya olunur, amma grid səthinə avtomatik
   köçmür (Mərhələ 1/7, yuxarı bax).
8. Startup-da default five-spot HƏM geologiya, HƏM quyu cədvəlinə yazılır
   (əvvəllər yalnız quyu cədvəlinə) — 7-ci bölmənin i/j-si artıq TAM
   geologiya bağlantısından asılı olduğu üçün zəruri oldu (Mərhələ 5).
9. Mərhələ 5 və 6 ayrı-ayrı yox, bir commit kimi yazıldı (Mərhələ 6).

## Mərhələ planı

- [x] 1 — təhlil
- [x] 2 — `GeologicalWell` + serialization + miqrasiya + testlər
- [x] 3 — `xy_to_ij` / `depth_to_k` + testlər
- [x] 4 — `validate_wells` + testlər
- [x] 5 — 2-ci bölmənin UI-si (+ `wells_to_dataset` adapteri)
- [x] 6 — 7-ci bölmə ilə birləşmə + ssenari generatoru
- [x] 7 — yekun: golden (UYĞUNDUR) + tam test dəsti (624 keçdi) + `run.bat` (qismən, yuxarı bax)
