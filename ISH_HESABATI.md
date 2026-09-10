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

---

# A3 — Bakli-Leverett analitik istinad əyrisinin düzəldilməsi

Bu mərhələdə düzəldilən şey fizika mühərriki DEYİL, **onu ölçən etalondur**
(`imex2d/simulation/analytical.py`). Sınıq etalon fizika səhvindən daha
təhlükəlidir, çünki gələcək səhvləri maskalayır.

## Səhv nə idi

Profil yığılmasında (köhnə sətir 44-51) cəbhədən sonra cəmi İKİ nöqtə
qalırdı — `(x_front, sw_shock)` və `(1.6·x_front, swi)`. `np.interp` və
matplotlib onların arasını XƏTTİ birləşdirdiyi üçün şok sıçrayışı əvəzinə
uzun mailli enən pillə çəkilirdi. Riyaziyyat (fraksion axın, Welge
toxunanı, `x_front`) düzgün idi — səhv yalnız yığılmada idi.

## Diaqnoz və qəbul meyarı: kütlə eyniliyi

Bakli-Leverett həllində DƏQİQ ödənir:

    ∫₀^∞ (Sw(x) − Swi) dx  =  q·t / (φ·A)  =  velocity · t

| | ∫(Sw−Swi)dx | gözlənilən | xəta |
|---|---|---|---|
| düzəlişdən əvvəl | 93.5635 | 75.0 | **24.75 %** |
| düzəlişdən sonra | 75.0000267 | 75.0 | **0.0000356 %** |

## Nə DƏYİŞMƏDİ (ölçülüb təsdiqlənib)

- `sw_shock = 0.5386941`, `x_front = 182.61886` — düzəlişdən əvvəl və sonra
  **bit-bə-bit eyni**. Müstəqil yoxlama: testin içində `scipy.optimize.brentq`
  ilə, ANALİTİK `df/ds` (zəncir qaydası, `np.gradient`-siz) üzərindən tapılan
  Welge kökü ilə `rtol=1e-4` daxilində üst-üstə düşür.
- Rarefaksiya yelpiyi və `np.gradient(f, s)` toxunulmadı.
- Mühərrikin çıxışı dəyişmədi: `nx = 60/120/240` üçün ədədi profilin öz həcm
  balansı düzəlişdən əvvəl və sonra eyni qaldı (71.137 / 72.968 / 73.950 m,
  75.0-ə yaxınsayır — yəni mühərrik onsuz da düzgün idi).

## Düzəliş

1. Profil dörd hissədən, ARTIQ SIRALANMIŞ şəkildə qurulur: süpürülmüş
   zona → rarefaksiya yelpiyi → şok → toxunulmamış zona.
2. Şok `np.nextafter(x_front, np.inf)` ilə ŞAQULİdir — fiziki qalınlığı
   sıfır, massiv isə ciddi artan. `x_front · 1.001` kimi süni ofset
   İŞLƏDİLMƏDİ: o, qalınlığı olan saxta keçid zonası yaradıb problemi
   kiçildər, aradan qaldırmazdı.
3. `np.argsort` **silindi**. Defolt `quicksort` stabil deyil; düzəlişdən
   sonra eyni `x`-də iki fərqli `sw` (`sw_shock` və `swi`) olacağı üçün bu,
   gizli səhvdən aktiv səhvə çevrilərdi.
4. `1.6` sehrli sabiti getdi → `length` parametri (verilməzsə `2·x_front`,
   YALNIZ göstərmə üçün). Çağıran `length = nx · dx` ötürür.
5. `breakthrough` bayrağı: `x_front ≥ length` olanda UI "cəbhə modeldən
   çıxıb" deyir, səssizcə mənasız faiz göstərmir.
6. `np.maximum.accumulate` saxlanıldı, amma düzəltdiyi nöqtələr SAYILIR və
   `monotonicity_corrections` + `notes` kimi qaytarılır — səssiz düzəliş
   qadağandır. Corey defoltlarında say sıfırdır.
7. Giriş yoxlaması əlavə olundu (`time`, `total_rate`, `area`, `porosity`,
   lözlüklər > 0; `swc ≤ Swi < 1 − sor`) — əvvəl heç bir yoxlama yox idi.

## Ölçmə metrikası (`main_window.run_validation`)

Köhnə metrika (`Sw > Swc + 0.01` olan SON hüceyrə) yayılmış cəbhənin ÖN
KƏNARIDIR və sistematik olaraq şişirdir. İndi üç metrika göstərilir:

1. **Əsas:** bütün profil üzrə RMS xəta (hüceyrə mərkəzlərində).
2. **Cəbhə:** `Sw = (sw_shock + Swi)/2` səviyyəsinin keçidi — yayılma bu
   nöqtəyə görə simmetrik olduğundan qərəzsizdir. Ölçdüm: ön kənar 7.33 %,
   orta nöqtə 4.53 %.
3. **Mühərrikdən asılı olmayan:** ədədi profilin `∫(Sw−Swi)dx` həcm balansı.

Validasiya modeli artıq SIXILMAZ flüidlə qurulur (`base.fluids`-dən yalnız
lözlüklər götürülür, `c = 0`) — BL sıxılmaz nəzəriyyədir və etalonda "kiçik
naməlum fərq" qalmamalıdır.

## ƏSL SÜBUT — grid yaxınsaması

Eyni fiziki məsələ (`length = 960 m`, `t = 250 gün`), üç şəbəkə:

| nx | RMS, düzəlişdən ƏVVƏL | RMS, SONRA |
|---|---|---|
| 60 | 0.0498 | 0.0397 |
| 120 | 0.0564 | 0.0231 |
| 240 | 0.0606 | 0.0163 |

Əvvəl xəta şəbəkə xırdalandıqca **artırdı** — çünki ölçü cihazının özü səhv
idi. İndi monoton azalır. `tests/test_analytical_bl.py::
test_refining_the_grid_reduces_error_against_analytical` bunu qoruyur.

## Yoxlama

- `tests/test_analytical_bl.py` — 24 yeni test (kütlə eyniliyi, dəqiq
  quyruq, ciddi artan məsafə, monotonluq, müstəqil Welge kökü, zaman
  miqyası, breakthrough, giriş yoxlaması, grid yaxınsaması,
  kapilyar/cazibə izolyasiyası). §8-ə uyğun olaraq kütlə eyniliyi və grid
  yaxınsaması testləri düzəlişdən ƏVVƏL yazılıb və QIRMIZI olduqları
  görülüb (16 uğursuz / 6 keçdi).
- `tests/test_physics.py::test_buckley_leverett_front_position` silinmədi,
  tolerantlığı 12 % → **9 %** sıxlaşdırıldı (faktiki ön-kənar xətası
  7.33 %). Metrikanın qərəzli olduğu docstring-də açıq yazıldı.
- `python tools/golden.py` → **ÜÇÜ DƏ UYĞUNDUR** (3.179645 / 8.411893 /
  6.766724 — dəyişməyib). Gözlənildiyi kimi: bu modul simulyasiya
  zəncirinin hissəsi deyil. Məhz buna görə səhv golden ilə tutulmamışdı —
  bu qeyd `tools/golden.py` başlığına da yazıldı.
- `python -m pytest -q` → **2058 keçdi, 1 buraxıldı, 1 uğursuz**. Uğursuz
  olan `test_phase_b_production_integration.py::test_differential_phase_a_
  vs_phase_b_permx_differ_and_why`-dır; A2 hesabatında artıq qeyd olunub və
  `git stash` ilə TƏMİZ ağacda yenidən yoxlanılıb — eyni sətir, eyni səbəb
  (`_cell_centres`-ə `geometry` yerinə `spec` ötürülür). A3 ilə əlaqəsi
  yoxdur.
- **Vizual yoxlama edilib** (`ValidationRenderer` başsız/offscreen
  matplotlib ilə real ədədi nəticə üzərində PNG-yə çəkilərək): analitik
  əyri indi cəbhədə ŞAQULİ düşür və sonra düz üfüqi `Swi = 0.20` xətti
  kimi gedir; köhnə mailli enən xətt YOX OLUB. Ədədi əyri yelpik zonasında
  demək olar üst-üstə düşür, yalnız cəbhə ətrafında yayılma qədər
  fərqlənir.
- `python app.py` → açılır, 12 tab qurulur, aralarında **Validasiya
  (B-L)**. (A1/A2-dən fərqli olaraq bu dəfə mühit tətbiqi işə saldı.)
  Tabın DAXİLİNDƏ düyməyə basıb rəqəmləri gözlə yoxlamaq isə interaktiv
  seans tələb edir — istifadəçi tərəfindən təsdiqlənməlidir.


---
---

# ⬇️ Aşağıdakılar `LAY-SIMULYATIR-MODELI-` reposundan birləşdirildi

**Birləşdirmə tarixi:** 10 sentyabr 2026

10 sentyabr 2026-da ayrıca `LAY-SIMULYATIR-MODELI-` reposu açılmışdı —
məqsəd sıfırdan yeni simulyator qurmaq idi. Audit göstərdi ki, bu
layihə (`lay-simulyatoru`) artıq işin böyük hissəsini ehtiva edir
(bax `AUDIT_2026-09-10.md`), ona görə iki repo birləşdirildi.

Hər iki tarixçə qorunub (`git merge --allow-unrelated-histories`).
Aşağıdakı qeydlər həmin repodan olduğu kimi köçürülüb.

---

---

## 10 sentyabr 2026 — Seans 1: Layihənin qurulması

### Görülən iş

1. **GitHub repo-su bağlandı**
   `https://github.com/eoraieor-tech/LAY-SIMULYATIR-MODELI-` (private, boş idi)
   Lokal qovluq: `C:\Users\LOQ\LaySimulyator`
   Branch: `main` · Remote: `origin`

2. **Sənədləşdirmə strukturu quruldu**
   - `README.md` — giriş nöqtəsi və naviqasiya
   - `ARCHITECTURE.md` — struktur (boş şablon, doldurulacaq)
   - `ROADMAP.md` — mərhələlər
   - `QARARLAR.md` — qərar jurnalı
   - `ISH_HESABATI.md` — bu fayl
   - `CLAUDE.md` — iş qaydaları
   - `docs/` — dərin mövzular üçün
   - `.gitignore`

### Qərarlar

- **Sənədlər Azərbaycan dilində yazılır.** Səbəb: komandanın iş dili budur.
  Kod və dəyişən adları ingiliscə qalacaq (standart praktika).
- **Şablonlar boş yaradıldı, uydurulmadı.** Layihənin texniki detalları
  hələ məlum deyil — uydurma məzmun sonradan səhv istiqamət verə bilər.
  Doldurulmamış yerlər ⏳ işarəsi ilə açıq göstərilib.

### Bilinməyənlər (növbəti seansda aydınlaşmalıdır)

| # | Sual |
|---|---|
| 1 | Simulyatorun məqsədi nədir — hansı prosesləri modelləşdirir? |
| 2 | Hansı proqramlaşdırma dili / texnologiya? |
| 3 | İnterfeys olacaq, yoxsa yalnız hesablama nüvəsi? |
| 4 | 2D, yoxsa 3D? Neçə faza (su/neft/qaz)? |
| 5 | Hazır layihədən (IMEX2D) nəsə götürülür, yoxsa tam təmiz başlanğıc? |

### Buraxılan iş

Heç nə — bu seansda yalnız struktur quruldu, kod yazılmadı.
Bu qəsdəndir: tələblər məlum olmadan kod yazmaq boş işdir.

---

## 10 sentyabr 2026 — Seans 2: Tələblər, mühit və mövcud kodun auditi

### Görülən iş

#### 1. Tələblər alındı

Sahibkar tam texniki plan verdi: **3D 3-fazalı (neft-su-qaz) hidrodinamik
lay simulyatoru + geoloji interpolyasiya platforması**, Python 3.11+ /
VS Code, **MPFA-O** diskretizasiyası (TPFA əvəzinə), 5-spot quyu şəbəkəsi,
RF/THP/BHP hesablanması, PyVista ilə 3D vizualizasiya. Altı mərhələlik
xəritə ilə birlikdə.

Bu, Faza 0-ın 1-4-cü açıq suallarını bağladı.

#### 2. Mühit yoxlandı — üç problem tapıldı

| Tapıntı | Nəticə |
|---|---|
| Python 3.11 **yoxdur** (3.14.7 və 3.12 var) | 3.12 seçildi (Q-03) |
| C++ alət zənciri yoxdur (MSVC, CMake, MinGW, LLVM) | C++ variantı praktiki deyil |
| **Smart App Control işlək** (`VerifiedAndReputablePolicyState = 1`) | imzasız DLL-lər bloklanır |

Smart App Control ölçmələri:

- Yeni yaradılan venv-lərdə (həm OneDrive-da, həm `C:\Dev`-də)
  `scipy.sparse.linalg` (`_superlu`), `scipy.interpolate` (`_specfun`),
  `scipy.spatial` (`_special_ufuncs`), PyKrige və VTK **bloklanır**
- **Sistem Python 3.14-də** isə SciPy-ın **bütün** modulları işləyir;
  yalnız VTK bloklanır
- Səbəb: köhnə fayllar reputasiya qazanıb, yeni yazılanlar yox

Yol boyu **iki səhv fərziyyəm oldu və düzəldildi**:

1. Əvvəlcə blokun səbəbini OneDrive qovluğunda gördüm — səhv idi.
   Yalnız `scipy.linalg` sınanmışdı; `scipy.sparse.linalg` sınananda
   məlum oldu ki, yeni yerdə də bloklanır.
2. Sonra "yalnız vizualizasiya bloklanır" dedim — bu da səhv idi;
   həqiqətdə xətti həlledici də bloklanırdı (yeni venv-də).

#### 3. Mövcud layihə tapıldı və auditdən keçirildi

Sahibkar bildirdi ki, proqramın bir hissəsi artıq yazılıb. `Downloads`
qovluğunda **14 zip versiyası** tapıldı (31 avqust → 10 sentyabr).
Ən yenisi (`(13).zip`) auditə verildi.

Bu, Faza 0-ın 5-ci sualını bağladı: **"IMEX2D" məhz bu layihədir.**

**Ölçülənlər:** 65 175 sətir Python (28 509-u test, 96 fayl), 27 sənəd,
təmiz hexagonal arxitektura, 0 TODO / 0 FIXME / 0 HACK / 0 çılpaq
`except:`.

**Testlər real icra edildi** (fərziyyə ilə deyil): **1 993 keçdi,
1 uğursuz, 1 ötürüldü, 1 xfail** — 7 dəq 48 san.

**Tapılan səhvlər:**

- **§8.3** — `pytest.importorskip("vtk")` pytest 8.2+ ilə işləmir
  (yalnız `ModuleNotFoundError`-da ötürür, burada `ImportError` gəlir).
  Nəticə: toplama dayanır, **1 995 testin heç biri işləmir**.
- **§8.4** — uğursuz test diaqnoz edildi: səhv **işlək kodda deyil,
  testin özündədir** (`_cell_centres`-ə `geometry` yerinə `spec`).
  Düzəliş tətbiq olunub yoxlanıldı → **23/23 keçdi**.
- **§8.1** — **MPFA-O yazılıb (1 699 sətir), test edilib, amma
  mühərrikə qoşulmayıb.** Yalnız testlər idxal edir; `impes_engine.py`
  MPFA-O ötürülsə qəsdən xəta atır; `jacobian.py` hər üz üçün dəqiq
  2 hüceyrə fərz edir.
- **§8.5** — THP / VFP modulu yoxdur (planda tələb olunur).

#### 4. Qaz fazası bərpa edildi

`A7_PLAN.md` göstərdi ki, üç fazalı black-oil işi demək olar
bitirilmişdi, sonra **v69-da tamamilə silinib** — çünki quyu öz BHP
həddinə yaxınlaşanda Nyuton yığılmırdı və kök tapılmamışdı. Bunun
üzərinə strateji dönüş edilib: fizika OPM Flow-a həvalə olunsun.

Kod 14 zip-in **heç birində yox idi**. Lakin GitHub reposu
(`eoraieor-tech/lay-simulyatoru`) public çıxdı və reponun **ilk
commit-i** məhz bu məqsədlə saxlanılıb:

```
66deca8   v67 - son 3 fazali versiya (qaz fazasi + matplotlib 3D)
```

Oradan **12 fayl / 4 118 sətir** çıxarıldı (2 165 mənbə + 1 953 test)
və `berpa/A7_qaz_fazasi/` qovluğunda saxlanıldı — Stone II, üç fazalı
qalıq tənlikləri, analitik Jakobian, Nyuton döngəsi, dəyişən keçidi.

### Qərarlar

- **Q-03** — Python 3.12 (3.11 yoxdur, 3.14 çox yenidir)
- **Q-04** — mövcud kod bazası atılmır, üzərində davam edilir
- **Q-05** — bu layihə üçün yeni venv yaradılmır; sistem Python 3.14
  işlədilir
- **Q-06** — Smart App Control söndürülmür

### Buraxılan iş

- Yeni repoda kod yazılmadı — audit nəticəsi strateji qərarı
  dəyişdirdiyi üçün qəsdən dayandırıldı
- `src/` altındakı boş modul skeletləri qaldı; hansı repoda davam
  ediləcəyi həll olunandan sonra taleyi müəyyənləşəcək
- VTK bloku həll olunmadı (qərar gözlənilir)
- MPFA-O qoşulmadı, qaz fazası əsas koda qaytarılmadı

### Açıq suallar (növbəti seansda həll olunmalıdır)

| # | Sual |
|---|---|
| 1 | **Öz fizikamız, yoxsa OPM Flow?** v69 qərarı yeni planla ziddiyyətdədir |
| 2 | **Hansı repoda davam edirik?** `lay-simulyatoru` (65k sətir) yoxsa `LAY-SIMULYATIR-MODELI-` (boş) |
| 3 | Qovluq təkrarı: `OneDrive\Desktop\L.S.M\...` və `C:\Dev\LSM` — hansı qalır |
| 4 | VTK bloku: Plotly/matplotlib, WSL2, yoxsa SAC söndürülsün |

Tam audit: [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md)


---

## 10 sentyabr 2026 — Seans 2 (davamı): Birləşdirmə, düzəlişlər, bərpa

### Sahibkarın strateji qərarları

1. **Öz fizikamız** — A7 (qaz fazası) bərpa edilib bitiriləcək.
   v69-un OPM Flow-a keçid qərarı ləğv olunur.
2. **İki repo birləşdirilir** — `lay-simulyatoru`-nun tarixçəsi
   saxlanılır, `LAY-SIMULYATIR-MODELI-`-nin sənədləri ora köçürülür.

### Görülən iş

#### 1. Repolar birləşdirildi

`git merge --allow-unrelated-histories` ilə **hər iki tarixçə qorundu**.
İş qovluğu: `C:\Dev\LSM`. Commit sayı: 52 → 58.

Yoxlanıldı — hər iki əsas commit əlçatandır:
`66deca8` (A7 mənbəyi) və `6ae83d7` (sənəd reposunun ilki).

Toqquşmaların həlli:

| Fayl | Qərar |
|---|---|
| `ARCHITECTURE.md` | köhnə saxlanıldı (51 KB real sənəd) |
| `requirements.txt` | köhnə saxlanıldı (real pinned siyahı) |
| `CLAUDE.md` | **birləşdirildi** — iş qaydaları + `graphify` |
| `ISH_HESABATI.md` | köhnə 644 sətir **toxunulmadı**, yenisi sonuna əlavə |
| `.gitignore` | birləşdirildi |

`src/` və `tests/__init__.py` çıxarıldı: birincisi `imex2d/` ilə
üst-üstə düşürdü, ikincisi köhnə repoda yox idi və pytest-in idxal
semantikasını dəyişə bilərdi. Hər ikisi `0315ae7`-də qalır.

#### 2. Auditdə tapılan iki səhv düzəldildi

- **§8.3** — `pytest.importorskip("vtk")` → `try/except ImportError` +
  `pytest.skip(allow_module_level=True)`. İki faylda.
- **§8.4** — `_cell_centres(model_b.grid, spec)` → `model_b.geometry`.

#### 3. A7 bərpası tamamlandı — 13 fayl, 4 323 sətir

**Öz səhvimin düzəlişi:** ilk bərpada 12 fayl / 4 118 sətir demişdim.
`tests/test_variable_switching.py` (205 sətir) buraxılmışdı — axtarış
açar sözlərim (`three_phase`, `stone`, `gas`) onun adına uyğun gəlmirdi.
v69 commit-lərinin diffi ilə tutuşdurma bunu üzə çıxardı.

Tamlıq yoxlandı: `66deca8`→`HEAD` arasında silinən fayl siyahısı (13) ilə
`berpa/A7_qaz_fazasi/` (13) **hərfi-hərfinə** üst-üstə düşür.

#### 4. Bərpa yolu tapıldı

Silmə mərhələli və adlandırılmış commit-lərlə aparılıb — bu, geri
qaytarmanı xeyli asanlaşdırır:

`f1037a5` (UI) → `b4ea22a` (application) → `36cf0db` (rendering) →
`01e95a6` + `8cb9c86` (well/newton) → `001cc12` → `00833db` (modullar) →
`ba59eb9` (qalıqlar)

#### 5. VTK bloku öz-özünə həll oldu

VTK indi sistem Python 3.14-də problemsiz işləyir (`9.7.0`), SAC isə
hələ də işlək (`= 1`). Səbəb: reputasiya buludda yoxlanılır və sorğu
tamamlanandan sonra icazə verilir. Bax `QARARLAR.md` → Q-07.

### Testlərin nəticəsi — birləşdirmədən SONRA

```
2 042 keçdi, 2 ötürüldü, 1 xfail  —  14 dəq 42 san  (exit 0)
```

Əvvəl 1 993 keçirdi; fərq (49) məhz açılan VTK testləridir.
**Uğursuz test yoxdur.**

### Buraxılan iş

- A7 hələ əsas kod bazasına qaytarılmayıb — yalnız `berpa/`-da saxlanılır
- MPFA-O hələ mühərriyə qoşulmayıb
- THP/VFP modulu yazılmayıb

### Açıq suallar

| # | Sual |
|---|---|
| 1 | OneDrive-dakı `LAY-SIMULYATIR-MODELI-` qovluğu silinsinmi? |
| 2 | GitHub-da repo adı dəyişdirilsinmi? (`gh` yoxdur, sahibkar özü etməlidir) |

---

## 10 sentyabr 2026 — Seans 3: Sahibkarın planının kodla tutuşdurulması

### Tapşırıq

Sahibkar 3D 3-fazalı MPFA-O simulyatoru üçün tam plan verdi (texnologiya
stəki, fayl strukturu, 6 mərhələ) və soruşdu: **"bu plan bizim layihədə nə
yerindədir, nə yerində deyil — hamısını araşdır."**

Metod: sənədə deyil, **koda və işlədilən testlərə** əsaslanmaq.

### Ölçülmüş vəziyyət

| Göstərici | Dəyər |
|---|---|
| `imex2d/` mənbə kodu | 38 324 sətir (129 fayl) |
| `tests/` | 30 562 sətir (100 fayl) |
| Toplanan test | 2 210 |
| Python | 3.14.7 · VTK 9.7.0 işləyir |

### ƏSAS TAPINTI — öz sənədimiz səhv idi

`ROADMAP.md`, `AUDIT §8.1` və `ISH_HESABATI` (Seans 2) deyirdi:
**"MPFA-O yazılıb, amma mühərriyə qoşulmayıb."**

**Bu doğru deyil.** Kod sübutları:

- `implicit/jacobian.py` — "PHASE 5B-2" analitik ÇOXNÖQTƏLİ Jakobian
  (`_flux_multipoint`, `_build_pattern_diag_only`) MÖVCUDDUR
- `implicit/residual.py:89` — `_multipoint` yolu işlək
- `tests/test_phase_d_mpfa_integration.py` — 21 test, o cümlədən
  `test_end_to_end_...mpfa_pressure_solve`,
  `test_tpfa_path_is_never_reached_when_mpfa_selected`,
  `test_mpfa_analytic_jacobian_matches_finite_difference`

**İşlədildi (bu seansda):** `test_mpfa_o.py` + `test_mpfa_o_global_assembly.py`
+ `test_phase_d_mpfa_integration.py` → **145 test keçdi, 29.5 san.**

**Səhvin kökü tapıldı:** `simulation/discretization.py` sənəd sətri hələ də
*"Future MPFA-O (HƏLƏ YOXDUR)"* və *"jacobian.py BU FAZADA DƏYİŞMİR"*
yazırdı — Phase 5B-2-dən sonra köhnəlmiş şərh. Audit bu şərhi oxuyub
yanlış nəticəyə gəlmişdi.

**Düzəldildi:** `simulation/discretization.py` sənəd bloku yenidən yazıldı
(tarixi qeyd saxlanıldı) və `jacobian.py`-dəki mövcud olmayan
`tests/test_mpfa_jacobian.py` istinadı real fayl adı ilə əvəzləndi.

**Real qalan boşluq daha kiçikdir:** MPFA-O FIM-ə qoşulub, amma
`ModelAwareSimulationService.create_engine()` `flux_discretization`
ötürmür → istifadəçi UI-dən MPFA-nı SEÇƏ BİLMİR.

### A7 qaz fazası — cari vəziyyət ölçüldü

HEAD (`f2f5c74`) v69 addım 5b-ni geri qaytarıb: üç fazalı modullar +
qaz PVT/SCAL/provider/equilibrium kod bazasındadır.

**İşlədildi:** `test_three_phase*.py`, `test_stone_relperm.py`,
`test_gas_pvt.py`, `test_variable_switching.py` → **143 test keçdi.**

Qalan: v69 addım 1–4b (UI, application, rendering, well_state,
standard_well/coupled_newton) hələ geri qaytarılmayıb. Ölçüldü —
`application/` və `ui/` qatında üç fazalı mühərrikə HEÇ BİR istinad yoxdur
(`create_engine()`-də qaz budağı yoxdur), `tests/test_gas_ui_wiring.py`
hələ `berpa/`-dadır.

### Planın texnologiya bəndləri — hökm

| Planda | Bizdə | Hökm |
|---|---|---|
| Python 3.11 | 3.14.7 | ❌ geriyə addım |
| RBFInterpolator | öz kriginq + SGS + fasiya | ❌ RBF zəifdir (variance vermir) |
| PyKrige | öz kriginqimiz | ❌ artıq asılılıq |
| `spsolve`/SuperLU/AMG | CG+ILU + CPR | ❌ bizimki güclüdür |
| PyVista | birbaşa VTK 9.7 (895 sətir, 49 test) | ⚠️ eyni VTK-nın örtüyü — köçürmə xərci var, qazanc yalnız hazır slice/volume API |
| Plotly | matplotlib + PyQt5 masaüstü | ❌ Plotly veb üçündür |
| pytest, NumPy/SciPy | ✅ | ✅ uyğun |

**Fayl strukturu:** plandakı düz `src/{geostats,grid,numerical,...}` yığını
bizim heksaqonal (onion) arxitekturamızdan geridir. Modul ADLARI isə
bire-bir uyğun gəlir — yalnız yerləri fərqlidir.

### Planın GÖRMƏDİYİ — bizdə olan, itirilməməli işlər

Fault modeli və transmissivlik çarpanları · GRDECL/Eclipse `.DATA`
giriş-çıxışı və OPM idxalı · history matching + həssaslıq (tornado) ·
fasiya modelləşdirməsi, SIS, SGS ansamblı, cross-validation · tam vahid
sistemi və vahid-invariantlıq testləri · `.imx` layihə faylı · PDF
hesabat · 12 tablı PyQt5 UI.

### Həqiqətən görüləsi iş (prioritetlə)

| # | İş | Planda | Həcm |
|---|---|---|---|
| 1 | **THP / VFP modulu** | ✅ 3.8 | orta, sıfırdan |
| 2 | A7-nin servis + UI-yə qaytarılması (v69 addım 1–4b geri) | dolayı | orta |
| 3 | Qaz Nyuton rəqsi (trust-region / per-cell line search) | ❌ | orta-böyük |
| 4 | MPFA-nın UI-dən seçilə bilməsi (`create_engine`) | dolayı | **kiçik** |
| 5 | Slice plane + cəbhə animasiyası + GIF ixracı | ✅ 6.2–6.3 | orta |
| 6 | Pcog (qaz-neft kapilyar təzyiqi) | ✅ 2.7 | kiçik |
| 7 | Zaman sıralarının CSV/JSON ixracı | ✅ 5.6 | kiçik |

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Köhnəlmiş kod şərhləri düzəldildi** — `discretization.py` və
   `jacobian.py`. Səbəb: məhz bu şərh bir dəfə auditi yanlış nəticəyə
   apardı; düzəltməmək eyni səhvin təkrarına imkan verərdi. Kodun
   davranışı DƏYİŞMƏDİ, yalnız sənəd sətirləri.
2. **`ROADMAP.md` statusları real ölçmə ilə dolduruldu** — `Modul`
   sütununda plandakı fərzi adlar (`grid/`, `numerical/`) REAL fayl
   yolları ilə əvəzləndi.
3. **Köhnə hesabat bölmələrinə toxunulmadı** (layihə qaydası). Auditdəki
   §8.1 mətni saxlanıldı, üstünə düzəliş qeydi əlavə olundu.

### Yoxlama

```
tests/test_mpfa_o.py + test_mpfa_o_global_assembly.py
  + test_phase_d_mpfa_integration.py   →  145 keçdi   (29.5 san)

tests/test_three_phase*.py + test_stone_relperm.py
  + test_gas_pvt.py + test_variable_switching.py  →  143 keçdi  (2.9 san)
```

⏳ **Bütöv dəst bu seansda YENİDƏN İŞLƏDİLMƏDİ** — son tam icra Seans 2-dədir
(2 042 keçdi, 14 dəq 42 san), qaz testləri bərpa olunandan sonra toplanan
say 2 210-dur.

### Açıq suallar

| # | Sual |
|---|---|
| 1 | PyVista-ya keçmək, yoxsa mövcud VTK kodunu saxlamaq? (tövsiyəm: saxlamaq) |
| 2 | Növbəti iş hansıdır — THP/VFP, yoxsa A7-nin UI-yə qaytarılması? |
| 3 | MPFA-nın UI-dən seçilməsi (kiçik iş) indi edilsinmi? |

---

## 10 sentyabr 2026 — Seans 3 (davamı): İcra planı

### Tapşırıq

Sahibkar üç problemi qəbul etdi və **hədəf proqrama çatmaq üçün icra
planı** istədi.

### Görülən iş

[ICRA_PLANI.md](ICRA_PLANI.md) yazıldı: M1–M8 qəbul meyarları və
B1–B7 iş blokları (hər biri müstəqil commit-lənə bilən, testlə bitən).

### Planlaşdırma zamanı ölçülən YENİ fakt — B2 gözləniləndən KİÇİKDİR

v69 silmə commit-lərinin diff-ləri nəhəng görünürdü (məs. `main_window.py`
4 077 sətir). **Səbəb tapıldı:** silmə `tools/patch_*.py` skriptləri ilə
edilib və həmin skriptlər faylların **sətir sonluğunu (CRLF↔LF)
dəyişib**. `git diff --ignore-cr-at-eol` ilə real ölçü:

| v69 addımı | Görünən | Real semantik |
|---|---|---|
| 1 — UI | 3 306 | **97** |
| 2 — application | 785 | **20** (+133 test faylı) |
| 3 — rendering | 547 | **20** |
| 4a — `well_state` | 556 | **16** |
| 4b — `standard_well` + `coupled_newton` | 443 | **443** |
| `001cc12` | 34 | **34** |

Yəni A7-nin qaytarılmasında **yeganə real böyük addım 4b-dir**.

Əlavə ölçmə: `standard_well.py` və `well_state.py` v69-dan **sonra heç
dəyişməyib** → toqquşmasız götürülə bilər. `panels.py` isə 8 dəfə
dəyişib — orada əl ilə birləşdirmə lazımdır.

Yoxlanıldı: altı v69 commit-inin **heç biri** düz `git revert` ilə
tətbiq olunmur (`git apply -R --check` hamısında xəta verir) — bu,
Seans 2-dəki gözləntini təsdiqləyir.

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Blok sırası B1 → B2 → B3, paralel B4/B5/B6, sonda B7** — səbəbi
   plan sənədində yazılıb. Sahibkar sıranı dəyişə bilər.
2. **GIF ixracı üçün yeni asılılıq götürülmür** — Pillow 12.3.0 artıq
   quraşdırılıb və animasiyalı GIF yaza bilir (`imageio`/`ffmpeg` lazım
   deyil). Səbəb: SAC mühitində hər yeni imzasız paket risk deməkdir
   (bax Q-06/Q-07).
3. **VFP-nin I versiyasında slip modelləşdirilmir** (no-slip qarışıq
   sıxlığı) — Hagedorn-Brown/Beggs-Brill II versiyaya qalır və sənəddə
   ⏳ ilə açıq göstərilir. Səbəb: layihə qaydası — uydurma məzmun
   yazılmır, natamamlıq açıq işarələnir.
4. **Hər blokda dəyişməz qoruma:** 2 fazalı 5-spot etalonu
   (`test_regression.py`) pozulmamalıdır.

### Açıq suallar

| # | Sual |
|---|---|
| 1 | Blok sırası təsdiqlənirmi? (xüsusən B1-in birinci olması) |
| 2 | Q-08 texnologiya qərarları təsdiqlənirmi? |
| 3 | B7-dəki SPE1 benchmark-ı əhatəyə daxil edirikmi? |

---

## 10 sentyabr 2026 — Seans 4: B1 icra olundu (MPFA-O istifadəçiyə açıldı)

### Tapşırıq

Sahibkar: "icra planını başla." Plandakı sıraya uyğun olaraq **B1**
götürüldü — ən kiçik iş, amma layihənin əsas texniki fərqini
istifadəçiyə açan halqa.

### Problem

MPFA-O nüvəsi (1 699 sətir) və onun tam implicit mühərrikə qoşulması
(Phase 5B-2) ARTIQ mövcud idi və testlərlə doğrulanmışdı. Lakin
`ModelAwareSimulationService.create_engine()` mühərrikə
`flux_discretization` ÜMUMİYYƏTLƏ ötürmürdü — nəticədə istifadəçi
MPFA-O-nu heç bir yolla işə sala bilmirdi. Yəni layihənin reklam
etdiyi əsas üstünlük **istifadəçi üçün mövcud deyildi**.

### Görülən iş

| Fayl | Dəyişiklik |
|---|---|
| `application/config.py` | `TPFA`/`MPFA_O`/`FLUX_SCHEMES` sabitləri, `SimulationConfig.flux_scheme` (defolt TPFA), `uses_multipoint_flux`, `validate()`-də yoxlama |
| `application/simulation_service.py` | `_flux_discretization()` + `_reject_incompatible_engine()`; mühərrikə diskretizasiya ötürülür |
| `application/serialization.py` | `.imx`-ə `flux_scheme` açarı, oxunuşda `data.get(..., TPFA)` |
| `ui/panels.py` | `NumericalPanel`-ə "Axın diskretizasiyası" combobox-u + izah etiketi, `flux_scheme_choice()`, `set_flux_scheme()` |
| `ui/main_window.py` | Layihə açılanda sxem sonuncu işə salınmadan bərpa olunur |
| `tests/test_flux_scheme_selection.py` | **YENİ — 16 test** |

### İcra zamanı üzə çıxan, planda OLMAYAN məsələ

`MPFAODiscretization()`-un **defolt sərhəd bağlanışı `DIRICHLET`-dir**,
qalıq qatı isə onu hələ dəstəkləmir (Phase 5B-2 məhdudiyyəti — sərhəd
π dəyərləri ötürülmür). Defolt konstruktorla qurulsaydı, hər MPFA-O
seçimi dərhal xəta verərdi.

**Həll:** servis MPFA-O-nu `NEUMANN_ZERO` bağlanışı ilə qurur.

**Niyə bu düzgün seçimdir, "susdurma" deyil:** simulyator onsuz da
AXINSIZ (no-flow) xarici sərhəd tətbiq edir — TPFA yolunda da belədir.
Yəni seçim fizikanı DƏYİŞMİR, mövcud şərti təkrarlayır. Doğrulama
testləri (`test_phase_d_mpfa_integration.py::NEUMANN`) də məhz bunu
işlədir. Səbəb kodun içində sənədləşdirildi.

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Xəta emalı ikiyə bölündü.** İlk yazdığım variant hər
   `NotImplementedError`-u "MPFA-O yalnız implicit mühərriklə işləyir"
   mesajına çevirirdi — bu, YANILDICI idi: yuxarıdakı sərhəd problemi
   də həmin mesajı verirdi. İndi:
   - IMPES + MPFA-O → mühərrik qurulmazdan ƏVVƏL, dəqiq mesajla rədd;
   - digər hallar → `residual.py`-in ORİJİNAL mesajı saxlanılır (orada
     məhz hansı xüsusiyyətin maneə olduğu yazılıb; öz sözümüzlə əvəz
     etsək istifadəçi səbəbi itirərdi).
2. **`flux_scheme` `SimulationConfig`-də saxlanılır**, modeldə yox.
   Səbəb: diskretizasiya modelin xassəsi deyil, İŞƏ SALINMANIN
   xassəsidir — eyni model həm TPFA, həm MPFA-O ilə hesablana bilər.
3. **Defolt TPFA olaraq qalır** (`default_flux_discretization()`
   dəyişmədi) — köhnə modellərin və `.imx` fayllarının nəticəsi
   dəyişməməlidir.
4. **Layihə açılanda sxem sonuncu işə salınmadan bərpa olunur.**
   Mühərrik seçimi (IMPES/IMPLICIT) belə bərpa olunmur, amma orada
   nəticəyə yazılan bir dəyər yoxdur; burada isə var.

### Yoxlama — ölçülmüş

```
tests/test_flux_scheme_selection.py                       16 keçdi
test_regression.py + config + serialization + ui +
  phase_d_mpfa + initial_saturation_map + implicit_engine  131 keçdi (45.9 san)
```

**Ən vacib qoruma keçdi:** `test_regression.py` — 2 fazalı 5-spot
etalonu **dəyişmədi**.

Uc-uca test (`test_choosing_mpfa_changes_results_under_anisotropy_end_to_end`)
göstərir ki, seçim ƏDƏDLƏRƏ çatır: 30°-fırladılmış tenzorda TPFA və
MPFA-O fərqli kumulyativ neft verir.

**Bütöv dəst (B1-dən sonra):**

```
2 225 keçdi, 1 ötürüldü, 1 xfail, 0 UĞURSUZ  —  634.9 san (10:35)
```

Əvvəlki tam icra (Seans 2) 2 042 keçmişdi; fərq bərpa olunan qaz
testləri (143) və B1-in yeni testləridir (16).

Qeyd: dəst 15 dəqiqədən 10 dəqiqəyə düşüb — səbəb ölçülməyib, ehtimal
ki maşın yükü. Paralelləşdirmə (B7.4) hələ də faydalıdır.

### Proqramın REAL işlədilməsi — B1 uc-uca təsdiqləndi

Testlərdən sonra tətbiqin özü işə salındı və sürüldü (yalnız idxal
yox — pəncərə quruldu, seçicilər dəyişdirildi, «MODELİ İŞƏ SAL»
basıldı, ekran şəkilləri çəkildi).

**Nəticə — 10×10, 300 gün, IMPLICIT + MPFA-O:**

```
converged = True · 22 addım · 4.9 san
orta Δt = 13.6 gün · orta 3.2 Nyuton iterasiyası · 0 təkrar
OOIP 61.2 min m³ · RF 49.59 % · Water cut 83.7 % · Su gəlişi 155 gün
```

**Nəticə — 12×12×4 (ÜÇÖLÇÜLÜ), 200 gün, IMPLICIT + MPFA-O:**

```
converged = True · 17 addım · 5.2 san · RF 33.05 %
VTK səhnəsi: 1130×679 piksel, 14 aktyor, 10 416 fərqli rəng
```

Yəni B1 yalnız testdə deyil, **istifadəçinin gördüyü proqramda da**
işləyir: combobox görünür, seçim mühərriyə çatır, MPFA-O ilə həqiqi
3D simulyasiya yığılır.

### İki mühit tapıntısı (gələcək iş üçün vacib)

1. **GUI `QT_QPA_PLATFORM=offscreen` ilə ÇÖKÜR (segfault).**
   Səbəb: VTK-nın Win32 OpenGL pəncərəsi offscreen platformada
   "failed to get valid pixel format" verir və proses ölür.
   `MainWindow` VTK widget-ini dərhal qurduğu üçün proqram
   ümumiyyətlə açılmır.
   **Nəticə:** GUI-ni sürmək üçün `QT_QPA_PLATFORM=windows` işlədilir.
   Testlər bundan təsirlənmir — onlar `MainWindow` qurmur (AST təhlili
   və ayrıca panellər).
   ⏳ CI-də GUI sürmək lazım olsa, bu maneə həll edilməlidir.

2. **`QWidget.grab()` VTK sahəsini TUTMUR** — 3D tabının şəklində
   görüntü sahəsi QARA çıxır. Bu, proqramın səhvi DEYİL: native
   OpenGL alt-pəncərəni Qt özü çəkmir.
   **Düzgün üsul:** `vtkWindowToImageFilter` ilə birbaşa VTK
   pəncərəsindən almaq — belə edildikdə şəkil tam düzgündür (grid,
   Sw sahəsi, INJ-1/PROD-1 quyuları, oxlar, rəng şkalası).

### Öz səhvimin qeydi

İlk sürücü skriptində RF-i 4958.81 %, su kəsirini 8369.73 % kimi
çap etdim. Proqram düzgün idi — **`TimeSeries.recovery_factor` və
`water_cut` onsuz da FAİZDƏ saxlanılır**, mən üstünə 100 vurmuşdum.
Doğru dəyərlər: RF 49.59 %, WCT 83.7 %. Kod dəyişdirilmədi.

### Buraxılan iş

- B1 bitdi. Növbəti blok sahibkarın seçimindən asılıdır (plan sırası
  ilə **B2** — A7 qaz fazasının servis və UI-yə qaytarılması).
