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

### Sonradan tapılan həqiqi bug (istifadəçi bildirdi, düzəldildi)

İstifadəçi bildirdi: proqramda 2-ci bölməni AÇMAQ (toolbox səhifəsini
göstərmək) proqramı bağlayırdı. Səbəb: `GeologyMapWidget.paintEvent`-də
`painter.drawEllipse(sx - _RADIUS, ...)` və `painter.drawText(sx + ..., ...)`
XAM PYTHON `float` DƏYƏRLƏRİ ilə çağırılırdı — PyQt5-in bu overload-ları
(4 mövqeli float arqument) qəbul ETMİR, yalnız `QRectF`/`QPointF` və ya
tam `int` qəbul edir. Nəticə: hər dəfə xəritə RƏNGLƏNƏNDƏ (yəni panel
göstəriləndə) Python `TypeError` PyQt5-in virtual metod override-ı
daxilində tutulmadan atılır, PyQt5 bunu `qFatal`/abort kimi işləyir və
bütün proqram bağlanır.

**Niyə mən bunu tutmadım:** əvvəlki skriptli sınaqlarımda `set_data()`
çağırırdım (bu, `self.update()` ilə YALNIZ təzələməni PLANLAŞDIRIR), amma
heç vaxt `app.processEvents()` və ya real `.show()` ilə faktiki
`paintEvent`-i işə salmırdım — ona görə bug görünməz qaldı. Dərs: PyQt
widget sınağı `update()`-dən sonra HƏMİŞƏ `app.processEvents()` və ya
`.repaint()` ilə tamamlanmalıdır, əks halda paint kodu heç vaxt icra
olunmur.

**Düzəliş:** `imex2d/ui/geology_map.py` — `drawEllipse` üçün `QRectF`,
`drawText` üçün `QPointF` işlədildi. Sonra bütün toolbox (8) və görüntü
(12) tablarını proqramlı şəkildə açıb-bağlayan bir skriptlə YENİDƏN
sınandı (`window.toolbox.setCurrentIndex(i)` + `app.processEvents()`
hər i üçün) — heç bir tabda xəta qalmadı. Tam test dəsti (624) və golden
yenidən keçdi.

### İkinci həqiqi bug (istifadəçi bildirdi: "quyu əlavə etmək istəyəndə bağlandı")

`GeologyPanel.__init__`-də `self.add_button.clicked.connect(self.add_row)`
birbaşa qoşulmuşdu. `QPushButton.clicked` siqnalı `bool checked` arqumenti
göndərir; `add_row(self, well=None)` bir SEÇİMLİ parametr qəbul etdiyi
üçün PyQt5 bu bool-u `well`-ə ötürür → `well.name` sətrində
`AttributeError: 'bool' object has no attribute 'name'` — bu, `clicked`
callback-i daxilində tutulmadan atılır və proqram bağlanır (dəqiq
istifadəçinin təsvir etdiyi kimi). Digər düymələr (`duplicate_button`,
`delete_button`, `centre_button`) təhlükəsiz idi, çünki bağlı olduqları
metodlar HEÇ bir əlavə parametr qəbul etmir (PyQt yalnız slotun qəbul
etdiyi qədər arqument ötürür).

**Düzəliş:** `_on_add_clicked(self): self.add_row()` — sıfır-arqumentli
"wrapper" metod əlavə edilib, düymə ona qoşulub (`lambda: self.add_row()`
DEYİL — bu, `test_ui_static.py`-dəki statik "erkən çağırış" yoxlamasını
yalançı-müsbətlə pozurdu, çünki `ast.walk` lambda daxilindəki çağırışı da
`__init__`-in "birbaşa" çağırışı kimi görür). Əlavə olaraq `add_row`-un
`blockSignals`/`insertRow`/`setItem` ardıcıllığı `try/finally` ilə
bağlandı ki, gələcəkdə bənzər bir xəta cədvəli daimi "siqnalsız" vəziyyətdə
qoymasın.

**Necə tapıldı:** proqramlı şəkildə HƏR düymənin `.click()`-i çağırılıb
(modal `QMessageBox`-lar aftomatik bağlanaraq), nəinki yalnız tab
keçidləri. Əvvəlki sınaq YALNIZ tab açılışlarını yoxlamışdı (bax yuxarı,
paintEvent bugu), buna görə bu bug ötürülmüşdü. Golden və tam test dəsti
(624) yenidən keçdi.

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

---

# A1 — Geologiya cədvəlindən gələn `SW` xəritəsinin simulyasiyaya qoşulması

## Problem

`SW` xəritəsi geologiya cədvəlindən modelə çatırdı (3D-də görünürdü), amma
HEÇ BİR mühərrik onu oxumurdu: ilkin doyumluluq həmişə skalyar
`InitialConditions.water_saturation`-dan qurulurdu, ona görə OOIP və RF
istifadəçinin verdiyi Sw məlumatını əks etdirmirdi.

## Həll

Mühərriklərə TOXUNULMADI — hər ikisi onsuz da `IInitializationProvider`-ə
hörmət edir, ona görə dəyişiklik yalnız provider qatında və onu quran
servisdədir. Yeni bayraq: `InitialConditions.use_saturation_map`
(defolt `False` → köhnə davranış). Prioritet cədvəli `domain/initial.py`
docstring-indədir; hər iki bayraq açıq olanda təzyiq hidrostatik qalır,
doyumluluq isə xəritədən gəlir (ölçmə modelin üstündədir) və keçid
zonasının əvəz olunduğu jurnala INFO kimi yazılır.

Fayllar: `domain/initial.py`, YENİ
`simulation/initialization/saturation_map.py` (saf oxuma/yoxlama funksiyası
+ iki provider), `application/simulation_service.py`
(`_initialization` — yeganə qoşulma nöqtəsi), `domain/reservoir_model.py`
(diaqnostika), `application/serialization.py` (bayraq `.imx`-ə yazılır),
`history/parameters.py`, `ui/panels.py`, `ui/main_window.py`.

## Öz təşəbbüsümlə verilmiş qərarlar

1. **Xəritə statistikası üçün ayrıca `count_outside_scal_limits`** — kəsmənin
   ÖZÜ mühərrikdə qalır (`relperm.saturation_limits()`), modul yalnız neçə
   hüceyrənin kəsiləcəyini sayır. Saf funksiyanın imzası tapşırıqdakı kimi
   `-> np.ndarray` qaldı, say ayrı funksiyadadır (diaqnostika da onu işlədir).
2. **§3.8 (history matching):** kodda ayrıca "ilkin Sw" `ParameterDefinition`
   YOXDUR — skalyar yalnız `_reconcile_initial_saturation` vasitəsilə (SOR/SWC
   dəyişəndə) toxunulurdu. Ona görə variant (a) funksiya səviyyəsində tətbiq
   olundu: xəritə rejimində reconcile heç nə etmir, xəbərdarlıq isə
   `standard_parameters()` qurulanda BİR DƏFƏ yazılır (hər `apply()`-da yox —
   optimallaşdırma yüzlərlə dəfə çağırır, jurnal dolardı).
3. **UI-də checkbox söndürülmə qaydası:** xəritə yoxdursa seçim SEÇİLƏ BİLMİR,
   AMMA artıq seçilibsə (məs. `.imx` faylından belə açılıb) söndürülmür —
   əks halda istifadəçi qutunu geri qaldıra bilməz və xəta vəziyyətində
   ilişib qalardı. Səssiz "geri qayıtma" yoxdur: vəziyyət etiketdə yazılır,
   model diaqnostikası isə ERROR verir.
4. **Diaqnostikada `SW` xəritəsi var, bayraq bağlı → WARNING** (tapşırıqda
   tələb olunduğu kimi) — məhz bu səhvin təkrarlanmasının qarşısını alır.

## Yoxlama

- `python -m pytest -q` → **1976 keçdi**, 1 uğursuz:
  `test_phase_b_production_integration.py::test_differential_phase_a_vs_phase_b_permx_differ_and_why`.
  Bu, bu işdən ƏVVƏL də uğursuzdur (təmiz `git stash` ağacında yoxlanılıb) —
  test `_cell_centres(grid, geometry)`-ə `geometry` yerinə `spec` ötürür,
  A1 ilə əlaqəsi yoxdur.
- `python tools/golden.py` → **ÜÇÜ DƏ UYĞUNDUR** (RF rəqəmləri dəyişməyib).
- Yeni testlər: `tests/test_initial_saturation_map.py` (38 test) — saf
  funksiya, 4 provider halı, hər iki mühərrik, OOIP/RF, kəsmə (+ kəsilən
  hüceyrə sayının JURNALA düşməsi), diaqnostika, `.imx` round-trip (köhnə
  açarsız sözlük daxil), history mühafizəsi, UI.
- §2.1 iddiası yoxlanıldı: `SyntheticGeologicalModelBuilder` yalnız
  `PORO`/`PERMX`/`PERMY`/`PERMZ` yaradır (`scenarios.py:127-130`), `SW`
  YARATMIR — yəni golden keysləri bu dəyişiklikdən təsirlənə bilməz.
- **GUI əl ilə yoxlanmayıb:** bu mühitdə `python app.py` açıla bilmir —
  offscreen rejimdə VTK `failed to get valid pixel format` verib
  segmentation fault-la dayanır. Bu, TƏMİZ ağacda da eynidir (yoxlanılıb),
  yəni mühit məhdudiyyətidir. Panel məntiqi (checkbox → `initial_conditions()`,
  skalyar sahənin sönməsi, etiket/tooltip) offscreen Qt widget testləri ilə
  örtülüb, amma real klikləmə axını (cədvələ Sw yaz → interpolyasiya → RF
  dəyişir) istifadəçi tərəfindən təsdiqlənməlidir.


---

# A2 — Quyulardan gələn lay üstü/altı (TOP/BOTTOM) səthlərinin grid həndəsəsinə qoşulması

## Problem

Geologiya cədvəlindəki «lay üstü»/«lay altı» dərinlikləri toplanır
(`geology_adapter.AREAL_TARGETS`), interpolyasiya olunur
(`property_maps["TOP"]/["BOTTOM"]`), hesabatda və 3D-də göstərilirdi — amma
GRID HƏNDƏSƏSİ onlardan QURULMURDU. `_surface()` yalnız `top_depth +
i·dip_x + j·dip_y` düz maili müstəvi verirdi, `CellGeometry.dz` isə LAY
ÜZRƏ (`(nz,)`) olduğu üçün sütundan sütuna dəyişən qalınlığı ifadə edə
bilmirdi. Nəticədə istifadəçinin doldurduğu struktur dərinlikləri məsamə
həcminə, cazibəyə, equilibration-a və perforasiya→K uyğunlaşdırmasına
SƏSSİZCƏ heç cür təsir etmirdi.

## Həll

Yeni həndəsə motoru YAZILMADI — `CornerPointGeometry` (770 sətir, hazır və
test edilmiş) İSTİFADƏ olundu. Boşluq yalnız konstruksiyada idi: yeni saf
modul `domain/structural_grid.py` interpolyasiya edilmiş TOP/BOTTOM
səthlərindən `(ncell, 8, 3)` təpə massivi qurur (`structural_nodes`,
`proportional` bölgü) və onu yoxlayır (`validate_structure`,
`structure_statistics`); qalanını `CornerPointGeometry.from_nodes()` edir.
Serializasiya (`nodes`) və VTK CPG dəstəyi onsuz da hazır idi.

`GeologicalGridSpec`-ə beş sahə: `structure_from_wells` (defolt `False` →
HEÇ NƏ dəyişmir), `thickness_source` (`"wells"`/`"constant"`), `layering`
(`"proportional"`), `on_zero_thickness` (`"error"`/`"clamp"`/`"deactivate"`),
`min_thickness`. `build()` İKİ MƏRHƏLƏLİ oldu: müvəqqəti düz həndəsə →
YALNIZ areal TOP/BOTTOM (`_interpolate_areal`) → həqiqi həndəsə
(`_build_structural_geometry`) → `model.geometry` → QALAN xassələr (artıq
düzgün dərinliklərlə).

Fayllar: YENİ `domain/structural_grid.py`, `application/geology_service.py`,
`domain/geological_model.py` (`structural_issues`),
`application/serialization.py` (həmin siyahı `.imx`-ə yazılır),
`ui/panels.py` (checkbox + qalınlıq mənbəyi comboxu + söndürmə),
`ui/main_window.py` (spec-ə ötürmə + quyu i/j/K yenidən hesablanması),
YENİ `tests/test_structural_grid.py` (45 test).

## Öz təşəbbüsümlə verilmiş qərarlar

1. **Areal interpolyasiya QƏSDƏN 2D-dir** (nümunə və hədəf Z-si sıfır).
   Tapşırıq §3.2 areal addımın "yalnız (x, y) işlətdiyini" fərz edirdi,
   amma `_estimate_layer` HƏMİŞƏ 3D-dir (Z sütunu əlavə edir) — yəni Z
   müvəqqəti həndəsədən gəlsəydi, `top_depth`/`dip_x`/`dip_y`/`dz` şaquli
   məsafə vasitəsilə kriginq çəkilərinə, oradan da QURULAN HƏNDƏSƏYƏ təsir
   edərdi və "bu parametrlər iştirak etmir" iddiası YALAN olardı. Bunu test
   aşkarladı (`test_structural_geometry_ignores_top_depth_and_dip_completely`
   əvvəlcə UĞURSUZ oldu), sonra düzəldildi. Nümunələr həm də `depth=None`
   nüsxə kimi ötürülür — CSV-dən gələn `depth` sütunu süni şaquli məsafə
   yaratmasın.
2. **`GeologicalModel.structural_issues`** — `on_zero_thickness="error"`
   halında model QAYTARILMALI (istifadəçi görsün), amma
   `ReservoirModelBuilder` QƏBUL ETMƏMƏLİDİR. Mövcud yeganə qapı
   `GeologicalModel.validate()`-dir; həndəsə validasiyası isə MÜSBƏT-amma-
   nazik sütunu tutmur (`min_thickness` istifadəçinin öz meyarıdır). Ona
   görə `completeness_issues()` ilə EYNİ məntiqli açıq sahə əlavə olundu və
   `.imx`-ə yazılır — əks halda yararsız layihəni saxlayıb açmaq qapını
   SÜKUTLA açardı (köhnə fayllarda açar yoxdur → boş siyahı).
3. **UI paneli `GridPanel` yox, `GridGeometryPanel` adlanır** (tapşırıqdakı
   ad köhnədir) — dəyişiklik həmin sinfə edildi. Yeni açarlar `values()`-ə
   ƏLAVƏ EDİLMƏDİ, ayrıca `structure_values()`-dədir: `values()` sintetik
   qurucuya `**` ilə birbaşa açılır (`geology_builder.build(**grid_values)`)
   və ona bu açarlar naməlumdur.
4. **`thickness_source="constant"` halında DZ söndürülmür** — orada
   qalınlıq MƏHZ DZ-dən gəlir; söndürülən yalnız tavan dərinliyi və
   mailliklərdir.

## Yoxlama

- `python -m pytest -q` → **2027 keçdi** (45 yeni test), 1 uğursuz:
  `test_phase_b_production_integration.py::test_differential_phase_a_vs_phase_b_
  permx_differ_and_why`. Bu, A2-dən ƏVVƏL də uğursuzdur (təmiz ağacda
  yoxlanılıb, eyni sətir, eyni səbəb: test `_cell_centres`-ə `geometry`
  yerinə `spec` ötürür).
- `python tools/golden.py` → **ÜÇÜ DƏ UYĞUNDUR** (RF: 3.179645 / 8.411893 /
  6.766724 — dəyişməyib).
- Fizika sübutları (§5) ayrıca testlərdədir: PV nisbəti = qalınlıq nisbəti
  (1:2:3:4), düz strukturda PV və TPFA transmissibillikləri Kartezian ilə
  MAŞIN DƏQİQLİYİNDƏ eyni, maili strukturda YANAL üzlərdə cazibə yaranır
  (düz halda sıfırdır), equilibration-da qalxıqda neft zonası qalın,
  `interval_layers` sütun-həssasdır, `.imx` round-trip bit-bə-bit eynidir.
- Uçdan-uca qəbul meyarı: eyni 4 quyu ilə struktur rejimi PV-ni
  1.22 mln m³ → 1.98 mln m³ dəyişir (OOIP və RF ilə birlikdə), checkbox
  söndürüləndə isə BİT-BƏ-BİT köhnə rəqəmə qayıdır.
- **GUI əl ilə yoxlanmayıb:** bu mühitdə `python app.py` açıla bilmir —
  offscreen rejimdə VTK `failed to get valid pixel format` verib
  segmentation fault-la dayanır. A1-də olduğu kimi, bu TƏMİZ ağacda da
  eynidir (`git stash` ilə yenidən yoxlanılıb), yəni mühit
  məhdudiyyətidir. Panel məntiqi (checkbox → söndürmə → tooltip →
  `structure_values()`) offscreen Qt widget testləri ilə, `_interpolate_
  geology` bağlantısı isə AST + saf köməkçi (`_describe_moved_wells`)
  testləri ilə örtülüb. 3D görüntüdə əyri layın GÖRÜNMƏSİ və
  interpolyasiyadan sonra quyu indeksləri dialoqu istifadəçi tərəfindən
  təsdiqlənməlidir.
