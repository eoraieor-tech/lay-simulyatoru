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

---

## 10 sentyabr 2026 — Seans 5: Proqram açıldı, REAL SƏHV tapıldı (PVT + FIM)

### Necə tapıldı

Sahibkar "proqramı aç" dedi. Proqram açıldı (PID 23280, 12 tab, 402 MB)
və sahibkar özü iki simulyasiya işə saldı. Jurnal göstərdi ki, ikincisi
**yığılmayıb**. Bu, testlərin tutmadığı, yalnız real istifadədə üzə
çıxan səhvdir.

Jurnaldan (`logs/imex2d.log`):

```
18:09:14  RUN-001  pvt: statik (PVT yoxdur)
          Tamamlandı: 81 addım, t = 1500 gün, RF = 16.65 %

18:10:21  [XƏBƏRDARLIQ] PROD-1: BHP (150 bar) doyma təzyiqindən
          (240 bar) aşağıdır — quyudibində qaz ayrılacaq.

18:10:46  RUN-002  pvt: correlation(API=32, γg=0.75, T=70°C)
          t = 0.0 gün: zaman addımı minimal həddə də yığılmadı.
          RF = 0.00 %
```

Proqram Python istisnası ilə çökmədi (jurnalda `Traceback` yoxdur) —
pəncərə bağlandı.

### Təkrarlandı və TƏCRİD OLUNDU — ölçülmüş

Sahibkarın ayarları ilə eynilə təkrarlandı:

| Model | Nəticə |
|---|---|
| 41×41, **PVT YOX**, FIM | ✅ yığıldı, 81 addım, 2.1 san |
| 41×41, **PVT VAR**, FIM | ❌ **t=0-da yığılmadı** |
| 21×21, PVT VAR, FIM | ❌ yığılmadı |
| 11×11, PVT VAR, FIM | ❌ yığılmadı (0.4 san) |

**Grid ölçüsündən ASILI DEYİL.** Sxem TPFA idi — yəni B1-in MPFA
işi ilə ƏLAQƏSİ YOXDUR.

Doyma təzyiqi dəyişdirilərək səbəb təcrid olundu
(ilkin təzyiq 250 bar, istismarçı BHP 150 bar, 11×11):

| Doyma təzyiqi Pb | FIM nəticəsi |
|---|---|
| **240 bar** (panelin DEFOLTU, BHP-dən yuxarı) | ❌ **yığılmır** |
| 140 bar (BHP-yə yaxın) | ✅ 27 addım |
| 100 bar (BHP-dən aşağı) | ✅ 25 addım |
| 50 bar | ✅ 22 addım |
| **240 bar, amma IMPES mühərriki** | ✅ 2 499 addım |

### Diaqnoz

Tam implicit (Nyuton) mühərriki **doyma təzyiqindən AŞAĞI bölgə
yaranan kimi** yığılmır. IMPES eyni halda işləyir.

Bu, A7-dəki qaz Nyuton problemi ilə **eyni ailədəndir** (Pb-dən aşağı
qeyri-xəttilik), lakin burada **qaz fazası ümumiyyətlə yoxdur** — yəni
səbəb üç fazalı koddan asılı deyil.

**Praktiki nəticə:** "PVT modelini işlət" seçimi panelin ÖZ DEFOLT
dəyərləri ilə (Pb=240, ilkin təzyiq 250, BHP 150) tövsiyə olunan
mühərriklə İŞLƏMİR.

### Niyə 2 225 test bunu tutmadı

Testlər PVT-ni ya başqa parametrlərlə, ya da IMPES ilə işlədir. Panelin
defolt kombinasiyası (Pb yuxarı + FIM) heç bir testdə yoxdur. Bu boşluq
B3-də bağlanmalıdır.

### Plan üçün nə dəyişir — B3 XEYLİ ASANLAŞDI

Əvvəl B3-ün (Nyuton möhkəmliyi) təkrarlanması üçün A7 qaz mühərriki
lazım idi (B2-dən sonra). **İndi 2 fazalı, 0.4 saniyəlik təkrarlanma
halı var** — B3 artıq B2-ni GÖZLƏMİR və hər cəhd saniyələrlə ölçülür.
`ICRA_PLANI.md` → B3 buna uyğun yeniləndi.

### Yol boyu edilən kiçik düzəliş

`create_engine()` jurnalı artıq axın sxemini də yazır. Səbəb: RUN-002-ni
diaqnoz edərkən jurnaldan hansı sxemin işlədildiyini müəyyən etmək
MÜMKÜN OLMADI — B1-dən sonra jurnalda bu boşluq yaranmışdı.

### Buraxılan iş

- Səhv TAPILDI və təcrid olundu, **düzəldilmədi** — düzəlişi B3-ün
  işidir və sahibkarın prioritet qərarını gözləyir.

---

## 10 sentyabr 2026 — Seans 6: B2 icra olundu (A7 qaz fazası servisə və UI-yə qaytarıldı)

### Planın bir fərziyyəsi SƏHV çıxdı

`ICRA_PLANI.md` → B2 deyirdi: "4b ən çətin addımdır (443 sətir) — əvvəl
edilir ki, qalanı üstünə otursun."

**Ölçüldü: `standard_well.py` və `coupled_newton.py` YALNIZ
`tests/test_standard_well.py`-dan çağırılır.** Nə iki fazalı
`FullyImplicitEngine`, nə də üç fazalı `ThreePhaseSimulationEngine`
onları idxal etmir — bu, hələ mühərriyə qoşulmamış, PARALEL bir quyu
modeli yoludur (OPM tipli, BHP naməlum dəyişən kimi).

**Nəticə:** 4b heç nəyi bloklamırdı. Qazı proqrama qaytaran əsl iş
addım 2 (application) və 1 (UI) idi — cəmi ~120 semantik sətir.
Sıra dəyişdirildi: 2 → 1 → 3 → 4a/4b → 001cc12.

### Görülən iş

| Addım | Fayl | Nə edildi |
|---|---|---|
| 2 | `application/simulation_service.py` | qaz budağı: PVT-də qaz varsa `ThreePhaseSimulationEngine` |
| 2 | `application/model_builder.py` | `gas_scal` parametri |
| 2 | `application/serialization.py` | `.imx`-ə `gas_oil_contact` |
| 1 | `ui/panels.py` | PVT: "Qaz fazasını aktivləşdir" · SCAL: qaz-neft əyriləri · Ədədi: qaz papağı (GOC) |
| 1 | `ui/main_window.py` | `gas_scal` ötürülməsi, SCAL tabında 3-cü ox, müqayisədə `gas_saturation` |
| 3 | `rendering/renderers.py` | `draw_gas` (Stone II önizləməsi) |
| 4a/4b | `standard_well.py`, `well_state.py`, `coupled_newton.py` | üç fazalı formaya qaytarıldı |
| — | `tests/test_gas_ui_wiring.py` | `berpa/`-dan qaytarıldı — **8 test keçir** |
| 001cc12 | `domain/reservoir_model.py`, `version.py` | qaz-agah diaqnostika + FEATURES sətri |

### Bərpa zamanı TAPILAN VƏ DÜZƏLDİLƏN ÜÇ SƏHV

Bunlar bərpanın "mexaniki köçürmə" olmadığını göstərir:

1. **`time_stepping.py` üç fazalı kütlə balansını AÇA BİLMİRDİ.**
   `mb_water, mb_oil = getattr(result, "material_balance", ...)` dəqiq
   İKİ dəyər açırdı; üç fazalı Nyuton isə ÜÇ verir (su, neft, qaz).
   Bu kod v69-dan SONRA yazılıb, ona görə qazı heç vaxt görməmişdi.
   Nəticə: üç fazalı mühərrik minimal Δt-yə çatanda `ValueError` atırdı,
   mühərrikin qoruyucu bloku isə onu səssizcə "yığılmadı"ya çevirirdi.
   **Düzəliş:** faza sayından asılı olmayan `max(balances)`.

2. **ACTNUM daraltması üç fazalıda yanlış ölçü verirdi.**
   `ActiveDofReduction`-un defolt `variables_per_cell` dəyəri İKİ
   fazalıdır (=2); üç fazalıda 3 olmalıdır. Defoltla buraxılsaydı Nyuton
   "xətti həlledici uğursuz" verərdi. **Düzəliş:** açıq ötürülür.

3. **`self.reservoir` atributu üç fazalı sinifdə YOXDUR** (`self.model`
   var). ACTNUM kodu iki fazalı versiyadan gəldiyi üçün belə yazılmışdı.

Üçü də ancaq **testləri işlədəndə** üzə çıxdı — statik köçürmə ilə
tapılmazdı.

### B1 ilə toqquşma — həll olundu

B1 `create_engine()`-ə `flux_discretization` açar sözü əlavə etmişdi və
bütün mühərriklər onu eyni imza ilə alır. `ThreePhaseSimulationEngine`
belə parametr tanımırdı → `TypeError`. Üç fazalı qalıq/Jakobian yalnız
TPFA üçün yazılıb, ona görə:
- mühərrik parametri qəbul edir, çoxnöqtəli olanda **AÇIQ rədd edir**;
- servis daha əvvəl, aydın mesajla dayandırır ("MPFA-O hələ üç fazalı
  mühərriklə işləmir — TPFA seçin").

Səssizcə TPFA-ya keçmək istifadəçinin seçdiyindən BAŞQA bir hesab
aparmaq olardı — qadağandır.

### Orijinal koda nisbətən İKİ QƏSDƏN FƏRQ

1. **`engine_factory` daimi dəyişdirilmir.** Orijinal kod
   `self.engine_factory = ThreePhaseSimulationEngine` yazırdı. Servis
   təkrar işlədiləndə (history matching, həssaslıq) qaz söndürülsə belə
   üç fazalı mühərrik seçili qalır və `ValueError` verirdi. İndi seçim
   yalnız həmin çağırışa aiddir (`try/finally`). Test yazıldı.
2. **Üç fazalı mühərrik uğurla bitəndə YEKUN MESAJ yazır.** Bərpa
   olunan kodda bu blok yox idi: qaz aktiv olanda istifadəçi status
   sətrində BOŞ mesaj görürdü. İki fazalı ilə eyni format verildi.

### Davranış dəyişikliyi (qəsdən, testi ilə birlikdə)

`test_producer_below_bubble_point_warns_even_with_gas_in_the_pvt_table`
→ `test_producer_below_bubble_point_is_not_warned_when_gas_is_modelled`.

Xəbərdarlığın mətni "qaz fazası MODELLƏŞDİRİLMİR" deyir. v69 qaz
mühərrikini sildiyi üçün bu HƏMİŞƏ doğru idi. B2 mühərriki qaytardı —
PVT-də qaz varsa qaz REAL modelləşdirilir, yəni mətn FAKTİKİ OLARAQ
YANLIŞDIR. Şərt (`not has_gas_phase`) və test birlikdə dəyişdirildi.
Qaz sütunu olmayan cədvəl üçün xəbərdarlıq dəyişməz qalır.

### Ölçülmüş nəticə — qaz yolu işləyir

| Yoxlama | Nəticə |
|---|---|
| PVT-də qaz yox → mühərrik | `FullyImplicitEngine` ✅ |
| PVT-də qaz var → mühərrik | `ThreePhaseSimulationEngine` ✅ |
| qaz + `gas_scal` → relperm | `StoneRelativePermeabilityProvider` ✅ |
| qaz → sonra qazsız (eyni servis) | yenidən `FullyImplicitEngine` ✅ (yapışqan deyil) |
| qaz + MPFA-O | aydın `ModelValidationError` ✅ |
| uc-uca 8×8, 400 gün, Pb=100 | `converged=True`, **RF 62.72 %** |
| müqayisə: eyni model, 2 fazalı | RF 62.86 % |

**Sərbəst qaz olmayan halda üç fazalı model iki fazalı ilə demək olar
eyni cavabı verir (62.72 ↔ 62.86 %)** — bu, fizikanın sağlamlıq
yoxlamasıdır.

### AÇIQ QALAN İKİ MƏHDUDİYYƏT (uydurulmur, açıq yazılır)

1. **Qaz papağı (GOC) verilmədikdə neft "ölü" başlayır (Rs = 0).**
   Bu, A7-nin QƏSDƏN seçilmiş mühafizəkar defoltudur və kodda belə
   sənədləşib: "Domain modelində ayrıca 'ilkin Rs' sahəsi yoxdur; bu,
   ən mühafizəkar seçimdir — xəyali qaz yaratmır."
   **Praktiki nəticə:** GOC seçilməyibsə OGIP = 0 və qaz heç vaxt
   ayrılmır — çünki ayrılacaq həll olmuş qaz YOXDUR. Sahibkarın M4
   meyarı (P < Psat → qazın ayrılması) bunun üçün domain-də **ilkin Rs
   sahəsi** tələb edir. ⏳ Bu, B2-nin əhatəsində DEYİL (v69 onu
   silməmişdi — heç vaxt olmayıb); ayrıca iş kimi qeyd olunur.
2. **Doyma təzyiqi istismarçının BHP-sindən yuxarı olanda ÜÇ FAZALI
   mühərrik də yığılmır** — eynilə iki fazalı kimi (Seans 5 tapıntısı).
   Ölçüldü: Pb=240, BHP=150 → hər iki mühərrik t=0-da dayanır.
   Yəni **qazın görünəcəyi rejim məhz B3-ün bloklayıcı olduğu
   rejimdir.** B2 bağlantını qurdu, faydalı olması B3-dən asılıdır.

### Öz təşəbbüsümlə verilmiş qərarlar

1. **`001cc12`-dən versiya geri qaytarılması ALINMADI** (69 → 67).
   Səbəb: irəli gedirik; yalnız məntiqi düzəliş və FEATURES sətri
   götürüldü.
2. **ACTNUM daraltması (v69-dan sonrakı iş) SAXLANILDI** və üç fazalı
   formaya uyğunlaşdırıldı — bərpa sonrakı təkmilləşdirməni İTİRMƏDİ.
3. **`test_standard_well.py` və `test_well_state.py` üç fazalı
   versiyaya qaytarıldı**, çünki mənbə faylları da qaytarıldı — ikisi
   ayrılmazdır.

### Yoxlama — bütöv dəst

```
2 235 keçdi, 1 ötürüldü, 1 xfail, 0 UĞURSUZ  —  729.7 san (12:09)
```

Seans 4-də (B1-dən sonra) 2 225 idi; fərq bərpa olunan
`test_gas_ui_wiring.py` (8) və üç fazalı quyu testlərinin artımıdır.

**Ən vacib qoruma yenə keçdi:** `test_regression.py` — 2 fazalı 5-spot
etalonu **dəyişmədi**.

### Buraxılan iş

- Domain-də **ilkin Rs sahəsi** yoxdur → GOC-suz qaz ayrılması mümkün
  deyil (yuxarıda §1). `ICRA_PLANI.md`-yə **B4b** kimi əlavə olundu.
- Pb > BHP rejimində yığılmama qalır → **B3**.

---

## 10 sentyabr 2026 — Seans 7: B3-A — kök səbəb tapıldı və düzəldildi

### Tapşırıq

Sahibkar: "B3 başla." Planın öz sırası ilə əvvəlcə ÖLÇMƏ aparıldı —
sərhədin harada olduğunu bilmirdik.

### Ölçmə 1 — sərhəd kəskindir

11×11, ilkin təzyiq 250 bar, istismarçı BHP 150 bar, tam implicit:

| Pb | Nəticə |
|---|---|
| ≤ 200 bar | ✅ yığılır (25–31 addım, 0.1 san) |
| 220 bar | ❌ 28 addımdan sonra ilişir |
| 240 bar | ❌ t = 0-da |

### Ölçmə 2 — üç fərziyyə RƏDD edildi

| Fərziyyə | Yoxlama | Nəticə |
|---|---|---|
| "Quyu BHP həddində rəqs" (A7_PLAN) | ən pis qalıq quyuda deyil, VURUCUDA və Pb-yə yaxın hüceyrələrdə | ❌ |
| "Zaman addımı böyükdür" | Δt 1000 dəfə azaldıldı (1.0 → 0.001) | ❌ eyni yerdə donur |
| "Cədvəl kobuddur" | nöqtə sayı 20 → 200 | ❌ heç bir fərq yoxdur |

### Ölçmə 3 — Nyuton RƏQS ETMİR, DONUR

Line search izlənildi. Pb=200-də klassik kvadratik yığılma:
`2.16 → 1.12 → 9.7e-2 → 6.7e-3 → 2.9e-5 → 4.6e-10`.

Pb=240-da isə 4-cü iterasiyadan sonra:

```
it   norm         |dP| qəbul edilən   TAM Nyuton addımı
 4   6.847e-03    0.81 bar            413 bar
 5   6.837e-03    0.62 bar            639 bar
 6   6.835e-03    0.64 bar            656 bar
```

**Tam Nyuton addımı 656 bar dəyişim istəyir — halbuki bütün lay
176–299 bar aralığındadır.** Qoruyucular (Appleyard chopping + line
search) onu düzgün olaraq 0.6 bara kəsir, ona görə hər iterasiyada
qalıq cəmi 0.05 % azalır. Yəni qoruyucular İŞLƏYİR; problem
Jakobiandadır.

### KÖK SƏBƏB

Neft tənliyinin təzyiq üzrə diaqonalı `−So·B'o/Bo²`-yə mütənasibdir.
Korrelyasiya cədvəlində Bo doyma təzyiqinə qədər ARTIR (qaz həll olur,
neft şişir), ondan yuxarı AZALIR (sıxılma) — yəni `dBo/dp` işarə
dəyişir və **yolda SIFIRDAN keçir**:

| Pb | Sıfırdan keçid | İlkin lay təzyiqi |
|---|---|---|
| 200 bar | 205–210 bar | 250 bar → uzaqdır ✅ |
| 240 bar | **245–250 bar** | 250 bar → **düz üstündə** ❌ |

Pb=240 halında bütün lay t=0-da məhz həmin sıfır nöqtəsində oturur →
diaqonal ≈ 0 → Jakobian təkləşir → kiçik qalıq nəhəng addıma çevrilir.

**Bu, hər müşahidəni izah edir:** niyə t=0-da sınır, niyə Δt kömək
etmir (təkləşmə zaman həddində deyil), niyə YALNIZ neft tənliyi pisdir
(su qalığı 7e-03, neft 1e+01), niyə Pb=220 sonradan sınır (lay boşalıb
həmin təzyiqə çatanda), niyə cədvəl həlledicilik dərəcəsi təsirsizdir.

### Bu, ədədi qüsur DEYİL — model natamamlığıdır

Doymuş budaq məhz **ayrılan qazı** təsvir edir. İki fazalı model o qazı
modelləşdirmir, ona görə tənliklərdə "mənfi neft sıxılması" kimi
görünür. Üç fazalı mühərrikdə qaz tənliyi bunu kompensasiya edir və
doymuş budaq DÜZGÜNDÜR.

Bu, üç fazalı mühərrikin niyə eyni yerdə sındığını da izah edir:
**ilkin Rs = 0 olduğu üçün o da faktiki olaraq iki fazalıdır.**

### Düzəliş

Qaz sütunları OLMAYAN cədvəldə doyma təzyiqindən aşağı Bo **ölü-neft
budağı** ilə əvəz olunur: undersaturated meyl aşağı uzadılır, yəni
"neftdən qaz ayrılmır, tərkib sabit qalır". Belədə `dBo/dp < 0` hər
yerdə qalır, sıxılma müsbətdir, Jakobian təkləşmir.

**Yerləşmə — QƏSDƏN application qatında.** `BlackOilPVTProvider`
defolt olaraq TƏMİZ interpolyatordur (mövcud müqavilə və testlər
dəyişmir); `dead_oil_below_bubble_point=True` bayrağını
`ModelAwareSimulationService._build_pvt_provider()` verir. Səbəb: "qaz
modelləşdirilirmi" sualı iş axınının qərarıdır, cədvəl
interpolyasiyasının deyil. Bayraq qaz sütunlu cədvəldə TƏSİRSİZDİR.

### Ölçülmüş nəticə

| Pb | xam cədvəl | düzəlişlə |
|---|---|---|
| 150 | ✅ 28 addım, RF 59.75 % | ✅ **eyni** (28, 59.75) |
| 200 | ✅ 31, 61.94 % | ✅ 31, 62.04 % |
| 220 | ❌ yığılmır | ✅ 32, **62.66 %** |
| 240 | ❌ yığılmır | ✅ 34, **63.21 %** |
| 300 | ❌ yığılmır | ✅ 36, **64.36 %** |

**Sahibkarın RUN-002-si (41×41, Pb=240, tam implicit):**

```
əvvəl : t = 0.0 gündə yığılmadı, RF = 0.00 %
indi  : converged=True, 89 addım, 3.8 san
```

Heç bir hüceyrə Pb-dən aşağı düşmürsə nəticə BİTƏ-BİT eynidir
(Pb=150 sətri) — bu, ayrıca testlə kilidləndi.

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Yalnız Bo düzəldilir, özlülük TOXUNULMUR.** Ölçüldü: özlülüyü də
   ölü-neft budağına keçirmək yığılmaya TƏSİR ETMİR, yalnız nəticəni
   dəyişir (Pb=300: RF 62.92 % → 64.36 %). Minimal müdaxilə prinsipi.
   ⏳ Sahibkar istəsə özlülük də keçirilə bilər — açıq sual kimi
   koda yazıldı.
2. **Düzəliş provider-in defoltu DEYİL, açıq seçimdir.** İlk yazdığım
   variant provider-i defolt olaraq dəyişirdi və `test_pvt.py`-ın iki
   testini sındırdı. Həmin testlər HAQLI idi — provider təmiz
   interpolyator kimi sənədləşib. Seçim application qatına köçürüldü.
3. **Log mesajı əlavə olundu** — düzəliş tətbiq olunanda istifadəçi
   jurnalda görür; səssiz məlumat dəyişikliyi olmur.

### Yoxlama

```
tests/test_dead_oil_below_bubble_point.py   11 keçdi (YENİ)
Bütöv dəst: 2 246 keçdi, 1 ötürüldü, 1 xfail, 0 UĞURSUZ (12:25)
```

Seans 6-da 2 235 idi; fərq 11 yeni testdir. `test_regression.py`
(2 fazalı 5-spot etalonu) yenə **dəyişmədi**.

### B3-dən NƏ QALIR

B3-A bitdi. Qalan iki sual yalnız B4b-dən sonra ölçülə bilər:

1. Üç fazalı mühərrik ilkin Rs verilən kimi Pb-dən aşağı düzgün
   işləyəcəkmi? (nəzəri olaraq bəli — qaz tənliyi kompensasiya edir)
2. A7_PLAN-dakı "quyu BHP həddində Nyuton rəqsi" ayrıca problem idi,
   yoxsa elə BU degenerasiyanın özü? Ölçülməyib.

---

## 10 sentyabr 2026 — Seans 8: B4b — ilkin həll olmuş qaz (Rs)

### Problem

Üç fazalı mühərrik qaz papağı (GOC) verilmədikdə nefti "ölü"
başladırdı: `third_variable = np.zeros(n)`, yəni Rs = 0. Koddakı şərh
bunu belə əsaslandırırdı: *"domain modelində ayrıca ilkin Rs sahəsi
yoxdur; bu, ən mühafizəkar seçimdir — xəyali qaz yaratmır."*

Nəticə zəncirvari idi: OGIP = 0 → təzyiq doyma təzyiqindən aşağı düşsə
BELƏ qaz ayrılmırdı (ayrılacaq həll olmuş qaz yox idi) → GOR = 0.
Yəni üç fazalı mühərrik **faktiki olaraq iki fazalı işləyirdi** və
sahibkarın M4 meyarı ("P < Psat olduqda qazın ayrılması") ödənmirdi.

### Həll

`InitialConditions.solution_gor` sahəsi əlavə olundu:

- `None` (defolt) — PVT cədvəlindən çıxarılır: `Rs = Rs_sat(min(P, Pb))`.
  Cədvəldə Rs onsuz da Pb-dən yuxarı sabitdir, lakin `min()` AÇIQ
  yazıldı ki, qeyri-standart (Eclipse idxalı) cədvəldə də düzgün işləsin.
  Sənaye standartı — Eclipse `EQUIL`/`RSVD` ilə eyni məntiq.
- Ədəd verilsə, bütün hüceyrələrdə həmin sabit Rs (laboratoriya ölçməsi).

Toxunulan fayllar: `domain/initial.py`, `implicit/three_phase_engine.py`
(`_initial_solution_gor`), `application/serialization.py`,
`ui/panels.py` (PVT tabında "İlkin Rs-i əl ilə ver"), `ui/main_window.py`.

**Qaz papağı halı da düzəldildi:** əvvəl GOC verilsə də, papağın
ALTINDAKI hüceyrələrdə Rs = 0 qalırdı. İndi orada da həll olmuş qaz var.

### Ölçülmüş nəticə — M4 ARTIQ İŞLƏYİR

8×8, 400 gün, istismarçı BHP 150 bar:

| Pb | yığıldı | OGIP | ilkin Rs | maks Sg | GOR (son) |
|---|---|---|---|---|---|
| 100 | ✅ | 2 070 955 | 52.8 | 0.0000 | 52.8 |
| 150 | ✅ | 3 098 350 | 86.1 | 0.0000 | 86.1 |
| **200** | ✅ | 4 038 447 | 121.9 | **0.0807** | **124.9** |
| 240 | ❌ | 4 721 237 | 152.1 | — | — |
| 300 | ❌ | 4 936 519 | 161.6 | — | — |

- **OGIP artıq müsbətdir** (əvvəl 0 idi).
- **Pb = 200-də qaz HƏQİQƏTƏN ayrılır:** sərbəst qaz doyumluluğu 0.081-ə
  çatır, GOR 121.9-dan 124.9-a qalxır. Bu, M4-ün tam olaraq tələb etdiyi
  davranışdır.
- **Pb = 100/150-də qaz ayrılmır** — və bu DÜZGÜNDÜR: ən aşağı təzyiq
  ~150 bar, heç bir hüceyrə doyma təzyiqinə çatmır. Fiziki nəzarət kimi
  ayrıca test yazıldı.

### B3-B ölçüldü — problem QALIR

Pb = 240 və 300-də üç fazalı mühərrik hələ də yığılmır.

Yəni B4b **M4-ü bağladı, lakin B3-B-ni bağlamadı.** Qaz tənliyi
indi işləyir, amma yüksək doyma təzyiqi rejimində yığılma hələ də
alınmır. Bu, B3-A-dakı Bo degenerasiyasının üç fazalı analoqu ola
bilər (orada cədvəl QƏSDƏN düzəldilmir, çünki qaz tənliyi kompensasiya
etməlidir) — **lakin bu, hələ ölçülməyib**, iddia etmirəm.

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Defolt "PVT-dən çıxar" seçildi**, "sıfır" yox. Səbəb: istifadəçi
   qaz xassələri olan PVT cədvəli veribsə, Rs = 0 heç vaxt onun
   istədiyi deyil. Köhnə `.imx` faylları da bu yolla açılır — orada
   nəticə onsuz da OGIP = 0 idi, yəni itirilən mənalı bir şey yoxdur.
2. **UI-da defolt AVTOMATİKDİR**, əl ilə vermək isteğe bağlıdır. Səbəb:
   düzgün dəyər cədvəldən çıxır; istifadəçini məcbur etmək səhv riski
   yaradar.
3. **`min(P, Pb)` açıq yazıldı**, sadəcə `pvt.solution_gor(P)` yox.
   Korrelyasiya cədvəlində fərq yoxdur, amma idxal olunan cədvəldə ola
   bilər.

### Yoxlama

```
tests/test_initial_solution_gor.py   13 keçdi (YENİ)
Bütöv dəst: 2 259 keçdi, 1 ötürüldü, 1 xfail, 0 UĞURSUZ (10:42)
```

Seans 7-də 2 246 idi; fərq 13 yeni testdir. `test_regression.py`
(2 fazalı 5-spot etalonu) yenə **dəyişmədi**.

---

## 11 sentyabr 2026 — Seans 9: Dəyişikliklərin ikinci maşına endirilməsi

### Nə edildi

Sahibkarın ikinci maşınında (`C:\Users\Acer\IMEX2D`) GitHub-dan son
dəyişikliklər endirildi. Bu seansda **yeni kod yazılmadı** — yalnız
mövcud işin bu maşına gətirilməsi və qeydə alınması.

### Tapıntı: son iş `main`-də DEYİL

`git fetch` göstərdi ki, `origin/main` dəyişməyib (`2c2efd8` — Seans 2-nin
sonu). Seans 3–8-in bütün işi **`a7-berpa` budağındadır** və `main`-ə
birləşdirilməyib. `main..origin/a7-berpa` fərqi:

| Commit | Nə |
|---|---|
| `479bba3` | docs: sahibkarın planı kod bazası ilə tutuşduruldu (Seans 3) |
| `118c752` | docs: hədəf proqrama çatmaq üçün icra planı (B1–B7) |
| `fed364f` | faza B1: MPFA-O istifadəçiyə açıldı (konfiqurasiya + servis + UI + `.imx`) |
| `7ef8cd2` | docs: B1 real tətbiqdə uc-uca yoxlanıldı + iki mühit tapıntısı |
| `1168a58` | tapıntı: PVT + tam implicit mühərrik doyma təzyiqindən aşağı bölgədə yığılmır |
| `b509d7f` | faza B2: A7 qaz fazası servisə və UI-yə qaytarıldı |
| `b4d9fbb` | faza B3-A: doyma təzyiqindən aşağı Bo degenerasiyası düzəldildi |
| `97c576b` | faza B4b: ilkin həll olmuş qaz (Rs) — M4 bağlandı |

Ölçü: **49 fayl, +11 948 / −4 644 sətir.** Onlardan 13 fayl yeni test
(`test_three_phase_residual.py`, `test_gas_pvt.py`, `test_stone_relperm.py`,
`test_initial_solution_gor.py`, `test_dead_oil_below_bubble_point.py`,
`test_flux_scheme_selection.py` və s.).

### Bu maşında edilən əməliyyat

```
git fetch origin
git checkout -b a7-berpa origin/a7-berpa      # yerli budaq yaradıldı
```

`main` **toxunulmadı** — hələ də `2c2efd8`-dədir. Birləşdirmə (`merge`)
edilmədi, çünki bu, sahibkarın qərarıdır (aşağıda açıq sual).

### Yoxlama

```
python -c "import app"      → OK (idxal zənciri təmiz, PyQt5 + matplotlib yüklənir)
imex2d.version.VERSION      → 69  (buraxılış: 2026-08-28)
```

### Açıq suallar ⏳

1. **`a7-berpa` → `main` birləşdirilsinmi?** B1–B4b işi 8 commit-dir və
   `main`-də yoxdur. Sahibkar qərar verməlidir: birləşdirilsin, yoxsa
   B7-yə qədər ayrı budaqda qalsın.
2. **Bütöv test dəsti BU maşında qaçırılmadı.** Seans 8-in qeydinə görə
   2 259 test keçir, lakin həmin ölçmə digər maşında (`C:\Dev\LSM`)
   aparılıb. Bu maşının `venv`-i ilə təsdiq edilməyib.
3. **Proqram bu maşında GUI ilə açılıb sınaqdan keçirilmədi** — yalnız
   idxal yoxlanıldı.

---

## 11 sentyabr 2026 — Seans 10: B3-B həll olundu (doymamış Bo qolu)

### Nə istənildi

Sahibkar B3-B-ni `xfail` kimi kilidləməyi RƏDD ETDİ: "qaz fazasının
Pb ≥ 240 bar şəraitində tam problemsiz işləməsini istəyirəm".

### Diaqnozun DƏQİQLƏŞDİRİLMƏSİ — əvvəlki hipotez natamam idi

Seans 9-da ölçmüşdüm ki, ölü-neft düzəlişini qaz cədvəlinə MƏCBURİ
tətbiq etsək Pb=240/300 yığılır — yəni kök səbəb B3-A ilə eynidir.
Bu doğru idi, LAKİN düzəlişin ÖZÜ yanlış olardı: 3 fazalı modeldə
ölü-neft Bo-su neftin şişməsini silir və qaz kütlə balansını pozar.

Bu seansda əsl səbəb tapıldı — `three_phase_newton.build_fluid()`:

```python
bo=self.pvt.oil_fvf(pressure)      # HƏR hüceyrədə DOYMUŞ qol
```

Doymamış hüceyrədə bu termodinamik olaraq YANLIŞDIR. Orada neftin
tərkibi sabitdir (`Rs` sərbəst primary dəyişəndir və `Rs < Rs_sat(p)`),
ona görə Bo həmin `Rs`-in doyma təzyiqindən başlayan SIXILMA qoluna
aiddir — `Rs_sat(p)`-in doymuş qoluna YOX.

İki nəticəsi vardı:

1. `dBo/dp` doymuş qolda Pb-də sıfırdan keçir → neft tənliyinin təzyiq
   diaqonalı (`−So·B'o/Bo²`) itir → Jakobian kilidlənir.
2. `∂N_o/∂Rs = 0` (`three_phase_residual.py`, `blocks[:, 1, 2]`) —
   neft tənliyi 3-cü dəyişəndən TAM qopmuşdu. Kodda yazılmışdı ki
   "qaz tənliyi kompensasiya edir"; **ölçmə bunu təkzib etdi** —
   kompensasiya edə BİLMƏZDİ, çünki əlaqə heç yox idi.

### Həll — sənaye standartı, saxtakarlıq yox

Doymamış qol Eclipse `PVTO` məntiqi ilə bərpa olundu:

    Bo(p, Rs) = Bo_sat(Pb(Rs)) · exp(c_o · (Pb(Rs) − p))

`c_o` cədvəlin ÖZ doymamış qolundan çıxarılır (`ln Bo` orada p-yə görə
xəttidir). `Rs` TOXUNULMUR — qaz kütlə balansı pozulmur.

Xassələri:

* `dBo/dp = −c_o·Bo < 0` HƏR YERDƏ → kilidlənmə yoxdur;
* `∂Bo/∂Rs ≠ 0` → neft tənliyi 3-cü dəyişənə bağlanır;
* doymuş hüceyrədə `Pb(Rs) = p` → `Bo = Bo_sat(p)`, yəni keçid
  KƏSİLMƏZ (ölçüldü: fərq 8·10⁻⁷);
* Pb-dən yuxarı düstur cədvəlin qurulduğu düsturun eynisidir → köhnə
  nəticələr DƏYİŞMİR.

### Toxunulan fayllar

| Fayl | Nə |
|---|---|
| `simulation/pvt/black_oil.py` | `saturation_pressure`, `oil_fvf_undersaturated` və törəmələri (yeni, additiv) |
| `implicit/three_phase_newton.py` | `_oil_fvf()` — vəziyyətdən asılı Bo |
| `implicit/three_phase_residual.py` | `bo_p`/`bo_rs` flüid sahələri; akkumulyasiya, axın və quyu Jakobianlarında `∂Bo/∂Rs` hədləri |

Köhnə çağırışlar üçün geri-dönüş saxlanıldı: flüid `bo_p`/`bo_rs`
verməzsə Jakobian əvvəlki davranışa qayıdır (`coupled_newton` və
flüidi əl ilə quran testlər toxunulmadı).

### Ölçülmüş nəticə — 8×8, ilkin 250 bar, BHP 150 bar, 400 gün

| Pb | əvvəl | indi | maksSg | GOR |
|---|---|---|---|---|
| 100 | ✅ RF 62.72 % | ✅ **bitə-bit eyni** | 0.0000 | 52.8 |
| 150 | ✅ RF 64.72 % | ✅ **bitə-bit eyni** | 0.0000 | 86.1 |
| 200 | ✅ RF 65.61 % | ✅ RF 65.71 % | 0.0807 | 124.9 |
| 220 | ✅ RF 64.98 % | ✅ RF 65.56 % | 0.0890 | 155.6 |
| **240** | ❌ **yığılmırdı** | ✅ **23 addım, RF 65.24 %** | 0.1033 | 199.4 |
| **300** | ❌ **yığılmırdı** | ✅ **22 addım, RF 64.41 %** | 0.1119 | 283.6 |

Pb = 100/150 dəyişmir, çünki orada heç bir hüceyrə Pb-yə çatmır və
`Rs` platoda qalır — yəni Bo elə cədvəlin öz qoludur.

Fizika monotondur: Pb artdıqca həll olmuş qaz, ayrılan sərbəst qaz
(maksSg) və GOR hamısı artır.

### Jakobianın doğrulanması

* `∂R_neft/∂Rs` diaqonalı: **0.0 → 1.31** (doymamış hüceyrələrdə)
* `∂Bo/∂p`, `∂Bo/∂Rs` mərkəzi fərqlə uyğun (5·10⁻⁹ və 1.6·10⁻²;
  ikincisi parçalı-xətti tərsin sınıq nöqtələrindən gəlir)
* akkumulyasiya bloku tam sonlu fərqlə yoxlanıldı — yeni testdə

### Yeni testlər

`tests/test_three_phase_high_bubble_point.py` — **23 test**: qolun
kəsilməzliyi, Pb-dən yuxarı dəyişməzlik, `dBo/dp` işarəsi, törəmələrin
sonlu fərqlə uyğunluğu, `Rs_sat` tərsinin dəqiqliyi, neft tənliyinin
3-cü dəyişənə bağlanması, Pb=240/300-də yığılma, M4 (qaz ayrılması),
monotonluq və REQRESSİYA qorumaları (Pb=100 dəyişməzliyi, 2 fazalı
yolun toxunulmazlığı).

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Ölü-neft düzəlişi 3 fazalıya TƏTBİQ EDİLMƏDİ.** Sahibkarın sualı
   "onu qaz balansını pozmadan necə tətbiq etmək olar?" idi. Cavab:
   tətbiq etmək LAZIM DEYİL — o, iki fazalı modelin natamamlığını
   örtən yamaqdır. Üç fazalı modeldə düzgün fizika onsuz da mövcuddur,
   sadəcə YANLIŞ QOLDAN oxunurdu.
2. **`c_o` cədvəldən çıxarılır**, korrelyasiya parametrlərindən yox —
   provider yalnız cədvəli tanıyır; idxal olunan (Eclipse) cədvəldə də
   işləsin deyə.
3. **Geri-dönüş yolu saxlanıldı** (`bo_p=None` → köhnə davranış) ki,
   `IPVTProvider` müqaviləsi genişlənməsin və öz provider-ini yazan
   kod sınmasın.

### Açıq qalan — AYRI qüsur ⏳

Tam Jakobianın sonlu fərqlə yoxlanışında bir hüceyrədə (istismarçı
quyusu, **qaz tənliyi ↔ Sw** elementi) **87 % xəta** var.

⚠️ Bu, BU seansın düzəlişindən GƏLMİR — ölçüldü: düzəlişdən əvvəl də
eyni idi (−9.22·10⁴ vs −9.25·10⁴). Yəni B3-B-dən ayrı, əvvəldən
mövcud qüsurdur. Yığılmanı bloklamır (22–23 addım), lakin Nyutonun
yaxınsama sürətini aşağı salır.

Ehtimal olunan mənbə: `ThreePhaseWellJacobian`-da RATE rejiminin açıq
sənədləşmiş sadələşdirməsi və ya sərbəst qaz debitinin `∂/∂Sw`
həddinin olmaması. **Ölçülməyib** — iddia etmirəm.

---

## 11 sentyabr 2026 — Seans 11: Proqram işə salındı, `run.bat` düzəldildi

### Kontekst

Bu seans **ikinci maşında** (yerli nüsxə) aparıldı. Başlayanda yerli
`main` commit `2c2efd8`-də ilişib qalmışdı — B1…B4b fazalarının
12 commit-i (Seans 3–10) yalnız remote-da idi. Sahibkarın tələbi
sadə idi: "proqramı run elə".

### Görülən iş

#### 1. Proqram işə salındı və doğrulandı

Kod dəyişikliyi aparılmadı — mövcud vəziyyət işlək halda yoxlanıldı.
Yoxlama **iki dəfə** edildi: əvvəl köhnə yerli kodla, sonra remote-dan
endirilən B1–B4b kodu ilə. **Hər ikisi təmiz açıldı.**

| Yoxlanılan | Nəticə |
|---|---|
| `.venv` asılılıqları (PyQt5, matplotlib, numpy, scipy) | ✅ mövcud |
| Proqramın başlaması | ✅ `IMEX-2D v69 başladıldı` |
| Əsas pəncərə | ✅ açıldı, 12 tab |
| Jurnalda xəta / istisna | ✅ yoxdur |

İşə salma əmri:

```
.venv\Scripts\pythonw.exe app.py
```

Pəncərənin başlığı: `IMEX-2D v69 · Geoloji modelləşdirmə və rezervuar
simulyasiyası`. Yaddaş: ~427 MB.

⚠️ **Yalnız açılış yoxlanıldı** — simulyasiya işlədilmədi, B1–B4b-nin
funksional nəticələri bu seansda **ölçülmədi**. Onlar Seans 4–10-da
öz maşınında artıq yoxlanılıb.

#### 2. `run.bat` — `.venv` yolu əlavə edildi

Skript venv-i `..\venv` və `venv`-də axtarırdı, faktiki mühit isə
kökdə `.venv`-dədir — yəni skript **heç vaxt işləmirdi**, həmişə
"XETA: venv tapilmadi" verirdi. Bax: `QARARLAR.md` → Q-10.

### Tapıntı 1 — GUI-ni fon prosesi kimi başlatmaq olmur

İlk cəhddə `python.exe app.py` **fon (background) prosesi** kimi
başladıldı. Nəticə:

- proses qalxdı, jurnala həm `başladıldı`, həm də `12 tab` yazıldı —
  **yəni Qt işləyirdi və `window.show()` çağırılmışdı**
- lakin `EnumWindows` ilə heç bir görünən pəncərə tapılmadı
- yaddaş 10 MB-da ilişib qalmışdı (normal işləyən instans ~427 MB)

Pəncərə interaktiv masaüstünə düşmür. Düzgün üsul — `Start-Process`
ilə birbaşa istifadəçinin sessiyasında başlatmaq.

### Tapıntı 2 — venv `pythonw.exe` yalnız ötürücüdür

`.venv\Scripts\pythonw.exe` əsl interpretator deyil (launcher stub);
o, əsas Python-u **ayrı proses** kimi çağırır:

```
pythonw.exe (venv, PID N)
  └── pythonw3.12.exe (Microsoft Store Python 3.12, PID M)  ← pəncərə BUNUNDUR
```

Prosesi axtaranda və ya dayandıranda **uşaq prosesə** baxmaq lazımdır —
`python.exe` / `pythonw.exe` adına görə filtr uşağı tapmır.

**Əlavə qeyd:** proqram `.venv` (Python 3.12) altında problemsiz işləyir.
Bu, Q-05-in "venv-də SciPy bloklanır" ölçməsinə **zidd deyil** — Q-07-də
blokun aradan qalxdığı artıq qeydə alınıb. VTK ayrı məsələdir və
açılışda idxal olunmur.

### Buraxılan iş

Bu seansda kod dəyişikliyi olmadı (yalnız `run.bat`) — `ROADMAP.md` və
`ARCHITECTURE.md` toxunulmadı, mərhələ statusu dəyişməyib.

### Açıq suallar

Seans 2-nin iki açıq sualı (OneDrive qovluğu, GitHub repo adı) **hələ
də açıqdır** — cavab verilməyib.

---

## 11 sentyabr 2026 — Seans 12: B4 (A variantı) — THP hesabatı

### Sahibkarın seçimləri (müzakirədən sonra)

| Sual | Qərar |
|---|---|
| Əhatə | **A variantı**: BHP məlum → yuxarı traverse → THP. B ayrı commit-də |
| Sürüşmə | v1-də **no-slip**; interfeys Beggs-Brill üçün açıq qalsın |
| RATE quyuları | v1-də **yalnız BHP rejimi**; Peaceman-ın tərsi sonraya |
| Lülə həndəsəsi | **tam şaquli** (TVD = MD) |

Hagedorn-Brown RƏDD EDİLDİ — onun `CNL`/`ψ` qrafiklərinin rəqəmsal
datası bizdə yoxdur və uydurulmayacaq.

### Kodda tapılan iki əsassız iddia

1. **`ROADMAP.md`-dəki "M5: BHP ✅" həqiqətdən zəif idi.** BHP çıxış
   deyil, GİRİŞ idi: `ControlMode.BHP`-də istifadəçi onu yazır,
   `SimulationResult`-da isə `well_bhp` sahəsi ÜMUMİYYƏTLƏ YOX idi.
   RATE rejimində mühərrik BHP-ni heç hesablamır. İndi `well_bhp`
   əlavə olundu və BHP həqiqətən qrafikə düşür.
2. **"M6: GOR ✅" — GOR heç bir renderer-də çəkilmirdi.** `grep`
   boş qayıtdı. İndi dashboard-a ayrıca panel əlavə olundu.

### Riyazi struktur

    dp/dz = [ ρ_qrav·g + f·ρ_ns·v_m²/(2d) ] / 1e5      [bar/m]
    THP   = BHP − Σ (Δp_qravitasiya + Δp_sürtünmə)

* **Sürtünmə:** Chen (1979) açıq düsturu (Colebrook-un iterasiyasız
  yaxınlaşması). Laminar budaq `64/Re`; 2000–4000 aralığında XƏTTİ
  keçid — sıçrayış B variantında Nyutonu pozardı.
* **Qarışıq sıxlığı:** kütlə axını lülə boyunca sabitdir, ona görə
  `ρ = kütlə/həcm` DƏQİQDİR. Qravitasiya `H_L` ilə, sürtünmə `λ_L`
  ilə çəkilir (sənaye standartı); no-slip halda ikisi üst-üstə düşür.
* **Sürətlənmə həddi DAXİL DEYİL** — ⏳ v2.

### Ölçülmüş: TƏK SEQMENT YETƏRSİZDİR

Plandakı ilkin eskiz tək seqment (orta təzyiqdə) nəzərdə tuturdu.
Ölçüldü (qazlı quyu, GOR 150, 2000 m):

| Seqment | THP, bar |
|---|---|
| 1 | 115.58 |
| 5 | 118.20 |
| 20 | **118.34** |
| 200 | 118.35 |

Tək seqmentin xətası **2.8 bar**-dır. Defolt 20 seçildi (200-dən fərq
0.006 bar). Bax `QARARLAR.md` → Q-10.

### Testin tutduğu SƏHV — B3-B-nin lülədəki analoqu

İlk versiyada `Bo` doymuş cədvəldən oxunurdu. Nəticə fiziki olaraq
QEYRİ-MONOTON çıxdı:

    GOR=  0 → THP 113.7
    GOR= 50 → THP 108.0   ← AŞAĞI düşdü, halbuki qaz sütunu
    GOR=150 → THP 118.4      yüngülləşdirməlidir
    GOR=400 → THP 155.0

Səbəb B3-B-dəki ilə EYNİ idi: doymamış neftə doymuş qolun Bo-su
verilirdi. Düzəliş: axının hasilat GOR-u `Rs`-in HƏQİQİ dəyərini
verir; `Rs < Rs_sat(p)` olduqda `oil_fvf_undersaturated()` işlədilir.
Düzəlişdən sonra: 78.8 → 90.6 → 118.3 → 155.0 (monoton).

### Sərhəd halları — uydurma dəyər YOXDUR

| Hal | Davranış |
|---|---|
| Quyu dayanıb | `thp = nan` — axan traverse təyin olunmayıb |
| Sütun BHP-ni üstələyir | `thp = nan` + "süni qaldırma lazımdır" |

Sıfıra "qısaldılmış" dəyər QƏSDƏN qaytarılmır: qrafikdə 0 bar real
ölçmə kimi görünərdi.

Ölçüldü — dərinliyə görə axan addımların sayı (24 addımdan):
1200 m → 24, 2000 m → 20, 3000 m → 11, 3500 m → 1, 5000 m → 0.
Keçid kəskin deyil, tədricidir: qaz ayrıldıqca sütun yüngülləşir.

### Toxunulan fayllar

| Fayl | Nə |
|---|---|
| `domain/tubing.py` | **YENİ** — `TubingGeometry` |
| `domain/wells.py` | `Well.tubing` (defolt `None`) |
| `interfaces/providers.py` | `IWellHydraulicsProvider` |
| `simulation/wellbore/friction.py` | **YENİ** — Chen (1979) |
| `simulation/wellbore/holdup.py` | **YENİ** — `IHoldupCorrelation`, `NoSlipHoldup` |
| `simulation/wellbore/traverse.py` | **YENİ** — çoxseqmentli marş |
| `simulation/wellbore/hydraulics.py` | **YENİ** — nəticəyə THP/BHP yazır |
| `simulation/results.py` | `well_bhp`, `well_thp` |
| `application/simulation_service.py` | post-proses qoşuldu |
| `application/serialization.py` | `.imx` açarı + geriyə uyğunluq |
| `ui/panels.py` | "Lülə — quyu başı təzyiqi (THP)" qrupu |
| `rendering/renderers.py` | dashboard 2×2 → 3×2: quyu təzyiqləri + GOR |
| `ui/main_window.py` | `_figure(3, 2)` |

### Öz təşəbbüsümlə verilmiş qərarlar

1. **THP POST-PROSESDİR, mühərrik toxunulmadı.** BHP rejimli quyuda
   quyu dibi təzyiqi onsuz da sabitdir, debit sıraları isə nəticədə
   var. Nyutona SIFIR risk. Test bunu kilidləyir: lüləli və lüləsiz
   qaçışın RF-i BİTƏ-BİT eynidir.
2. **Lülə parametrləri quyu cədvəlinə DEYİL, ayrıca qrupa qoyuldu.**
   Cədvəl artıq 10 sütundur; v1 onsuz da bütün istismarçılara eyni
   həndəsə tətbiq edir. Per-quyu lülə VFP cədvəlləri gələndə mənalı
   olacaq (⏳).
3. **Renderer düzümü ADAPTİVDİR** (2×2 və 3×2) — mövcud testlər və
   çağırışlar dəyişmədən keçir.

### Yoxlama

```
tests/test_wellbore_thp.py   26 keçdi (YENİ)
```

Sərhəd testləri analitik etalonla tutuşdurulur: sürtünməsiz hədd
hidrostatik sütuna 1e-9 dəqiqliklə bərabərdir; Chen ↔ Colebrook
(iterasiya ilə həll olunan müstəqil etalon) fərqi < 0.5 %.

### Qalan ⏳

* B variantı — `ControlMode.THP` (açıq birləşmə)
* RATE quyularında BHP (Peaceman-ın tərsi)
* Beggs-Brill sürüşməsi
* Sürətlənmə həddi
* Eclipse `VFPPROD` cədvəl idxalı (`domain/vfp.py`, `io/vfp_io.py`)
* Əyri (deviated) quyu — MD→TVD profili

---

## 11 sentyabr 2026 — Seans 13: Son dəyişikliklər əsas maşına endirildi

### Nə istənildi

Sahibkar: "Ən son dəyişiklikləri məndə də et." Yəni əsas maşında
(`C:\Dev\LSM`) uzaq repodakı bütün yeni iş tətbiq olunsun. Bu seansda
**yeni funksional kod yazılmadı** — endirmə, mühit yoxlanışı və iki
köhnəlmiş ROADMAP sətrinin düzəldilməsi.

### Vəziyyət: yerli repo ÜÇ commit geri idi

Seans başlayanda yerli `a7-berpa` `97c576b`-də (Seans 8-in sonu) idi,
iş ağacı təmiz. `git fetch --all` göstərdi:

| Budaq | Əvvəl | Sonra |
|---|---|---|
| `a7-berpa` | `97c576b` | `ad6104e` (+3 commit) |
| `main` | `2c2efd8` | `3b3e66b` (+12 commit) |
| `b4-thp` | *yox idi* | `963358f` (**YENİ budaq**) |

Diqqət: Seans 9-un açıq sualı (`a7-berpa` → `main` birləşdirilsinmi)
**həll olunub** — `main` artıq B1–B3-B işini daşıyır (`3b3e66b`).

### Endirilən iş (`97c576b` → `963358f`)

**26 fayl, +2 542 / −29 sətir.**

| Blok | Nə gətirdi |
|---|---|
| **B3-B** (`8a5ecce`) | üç fazalı mühərrik yüksək doyma təzyiqində yığılır — `three_phase_newton.py`, `three_phase_residual.py`, `black_oil.py` + `test_three_phase_high_bubble_point.py` (23 test) |
| **B4-A** (`97979be`) | quyu başı təzyiqi (THP) — YENİ `simulation/wellbore/` paketi: `friction.py` (Chen 1979), `holdup.py`, `traverse.py`, `hydraulics.py`; `domain/tubing.py`; UI + dashboard 3×2 + `.imx` açarı + `test_wellbore_thp.py` (26 test) |
| **run.bat** (`3b3e66b`) | `.venv` yolunu tanıyır (Seans 11 düzəlişi) |
| `ad6104e` | yerli qrafik faylı `.gitignore`-a |

### Bu maşında edilən əməliyyat

```
git fetch --all
git merge --ff-only origin/a7-berpa      # a7-berpa: 97c576b -> ad6104e
git fetch . origin/main:main             # main:      2c2efd8 -> 3b3e66b
git checkout -B b4-thp origin/b4-thp     # ish burada davam edir
```

Üç budağın hamısı `origin` ilə eynidir. Aktiv budaq: **`b4-thp`**.

### Mühit yoxlanışı — Seans 9-un 2-ci və 3-cü açıq sualı bağlandı

Bu maşında `.venv` **yoxdur**; `run.bat`-ın gözlədiyi üç yerin heç
birində venv tapılmır. Lakin sistem Python-u (3.14.7) bütün asılılıqları
daşıyır:

```
PyQt5 OK · vtk OK · numpy OK · scipy OK · matplotlib OK · pytest OK
import app  → OK,  imex2d.version.VERSION = 69
```

**Bütöv test dəsti BU maşında qaçırıldı** (Seans 9-da bu edilməmişdi):

```
python -m pytest -q
2308 keçdi · 1 skip · 1 xfail · 0 XƏTA   (648 s)
```

Seans 8-də 2 259 idi; fərq 49 testdir — məhz B3-B (23) və B4-A (26)
ilə gələnlər. Yəni iki maşının nəticəsi uc-uca uyğundur.

### Düzəldilən: ROADMAP-da iki köhnəlmiş sətir

B4-A THP-ni gətirsə də, `ROADMAP.md`-də status `❌` qalmışdı:

| Sətir | Əvvəl | İndi |
|---|---|---|
| **3.8** VFP: hidrostatik sütun + sürtünmə → BHP ↔ THP | `❌ YOXDUR` | `🟡` BHP → THP ✅ (B4-A), tərs istiqamət ⏳ B4-B |
| **6.5** Dashboard: debet, BHP, RF | `🟡 ... THP ❌` | `✅` THP/BHP paneli (3×2 düzüm) |

3.8 QƏSDƏN `✅` edilmədi: `ControlMode.THP` (tərs istiqamət — THP-dən
BHP tapmaq) hələ yoxdur, o B4-B-dir.

### Açıq suallar ⏳

1. **`b4-thp` → `main` birləşdirilsinmi?** B4-A `main`-də yoxdur.
   Sahibkarın qərarıdır.
2. **`.venv` bu maşında qurulsunmu?** Hazırda sistem Python-u işləyir,
   amma `run.bat` venv axtarır və tapmayanda "XETA: venv tapilmadi"
   verib dayanır. Yəni `run.bat` ilə proqram BU maşında açılmır —
   `python app.py` ilə açılır.
3. **Proqram bu maşında GUI ilə açılıb sınanmadı** — yalnız idxal və
   testlər yoxlanıldı.
## 11 sentyabr 2026 — Seans 14: "Parametrlər təsir etmir" — görünürlük düzəlişi

### Bildiriş

Sahibkar: interfeysdən neft lözlüyünü və φ/K dəyişdikdə RF tərpənmir.

### Audit nəticəsi — MÜHƏRRİK SAĞLAMDIR

Ölçüldü, fərziyyə qurulmadı. **Sahibkarın fərziyyəsi qismən yanlış idi:**
defolt (sintetik) yolda φ və K **əla işləyir**:

| Parametr | Dəyər | RF % | Kum. neft |
|---|---|---|---|
| məsaməlilik | 0.10 → 0.35 | 63.23 → 53.57 | 11 260 → 33 392 |
| keçiricilik | 10 → 1000 mD | 10.72 → 66.82 | 4 201 → 26 181 |
| neft lözlüyü (PVT söndürülü) | 1 → 30 cP | 64.92 → 22.49 | — |

**Hardcode rəqəm TAPILMADI.** Yoxlanıldı:

* `implicit/residual.py:112-128` — təmiz `pvt is None` budaqlanması
* `discretization.py:102` — `model.units.darcy_constant` (vahid sistemi
  sabiti, UNITS.md sənədləşdirib, audit əl hesabı ilə təsdiqləyib)
* `three_phase_newton.py:131,177` — `fluids.water_viscosity` yalnız
  `pvt is None` halında

### Tapılan ÜÇ səssiz üstələmə mexanizmi

**A — PVT açıq olanda panel lözlüyü ÖLÜDÜR.**
`μo = 1 cP → RF 64.8206`, `μo = 30 cP → RF 64.8206` (bitə-bit eyni).
Səbəb dizayn üzrədir: PVT cədvəli paneli üstələyir. PVT açıq olanda
əsl düymə **API və temperaturdur**:

| Düymə | Dəyər | μo(250 bar) | RF % |
|---|---|---|---|
| panel μo | 1 → 30 | 1.03 (dəyişmir) | 64.82 (dəyişmir) |
| PVT: API | 15 → 50 | 6.36 → 0.34 | 49.01 → 67.61 |
| PVT: temperatur | 40 → 120 °C | 1.83 → 0.66 | 62.15 → 66.15 |

**B — interpolyasiya olunmuş geologiya keşlənir.**
`_geology_model_from_wells` YALNIZ "İnterpolyasiya et" düyməsi və
layihə açılışı ilə təyin olunur; φ/K dəyişəndə etibarsızlaşdırılmırdı.
Keşin özü DÜZGÜNDÜR (böyük gridi hər klikdə interpolyasiya etmək
olmaz) — problem istifadəçiyə heç bir işarə verilməməsi idi.

`mark_stale()` mexanizmi ARTIQ MÖVCUD idi, lakin yalnız geologiya
CƏDVƏLİ redaktə olunanda çağırılırdı.

**C — GRDECL idxalında φ/K faylın xəritələrindən gəlir.**
Geologiya panelinin hesabatında yazılırdı, lakin süxur panelinin
sahələri açıq və redaktə edilə bilən qalırdı.

Sahibkarın halı: **A + B** (GRDECL idxal etməyib).

### Edilən — MÜHƏRRİYƏ TOXUNULMADI

| # | Düzəliş |
|---|---|
| 1 | PVT açıq olanda μw/μo/lözlük vahidi/Bo **bozarır** + izah: dəyərlər PVT tabından gəlir, API/temperaturu redaktə edin |
| 2 | φ/K və ya grid dəyişəndə interpolyasiya **köhnəlmiş** işarələnir → mövcud banner ("Nəticə köhnəlib — 'İnterpolyasiya et' basın") çıxır |
| 3 | GRDECL idxalında φ/K sahələri bozarır + izah |
| 4 | **Həssaslıq test dəsti** — əlaqə gələcəkdə səssizcə qırılsa tutulsun |

Toxunulan fayllar: `ui/panels.py` (`RockFluidPanel.set_context()`,
yeni `geology_changed` siqnalı, `context_note` etiketi),
`ui/main_window.py` (kontekst yenilənməsi + `_mark_geology_stale`).

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Ayrıca `geology_changed` siqnalı quruldu**, mövcud `changed`
   işlədilmədi. Səbəb: `changed` lözlük dəyişəndə də atəş açır, o isə
   geologiyaya TƏSİR ETMİR — interpolyasiyanı köhnəlmiş saymaq yanlış
   siqnal olardı.
2. **Keş avtomatik SİLİNMİR, yalnız köhnəlmiş işarələnir.** Orijinal
   niyyət (böyük gridi hər klikdə yenidən hesablamamaq) qorunur.
3. **Grid paneli də köhnəlmə siqnalına qoşuldu** — eyni qüsur ona da
   aiddir (nx/dx dəyişəndə keşlənmiş model qalırdı).
4. **Bozarma dəyəri SİLMİR** — söndürülmüş widget-dən `fluids()` hələ
   də oxuyur, yəni heç bir hesablama davranışı dəyişmir. Ayrıca testlə
   kilidləndi.

### Yoxlama

```
tests/test_parameter_sensitivity.py   14 keçdi (YENİ)
```

9 fizika testi (parametr → RF əlaqəsi) + 5 UI testi (bozarma, banner,
siqnal seçiciliyi, dəyərin itməməsi).

---

## 12 sentyabr 2026 — Seans 15: Seans 14 işi əsas maşına endirildi

### Nə istənildi

Sahibkar: "Bütün dəyişiklikləri məndə et." Yəni bu maşında
(`C:\Dev\LSM`) uzaq repodakı bütün yeni iş tətbiq olunsun. Seans 13-ün
eyni tapşırığıdır, bir mərhələ sonra. Bu seansda **yeni funksional kod
yazılmadı** — endirmə, mühit yoxlanışı və bütöv test dəsti.

### Vəziyyət: yerli `main` altı commit geri idi

Seans başlayanda aktiv budaq `b4-thp` (`18a1b88`) idi, iş ağacı təmiz.
`git fetch --all --prune` göstərdi:

| Budaq | Əvvəl | Sonra |
|---|---|---|
| `main` | `3b3e66b` | `9d3ac51` (+6 commit) |
| `b4-thp` | `18a1b88` | `18a1b88` (dəyişməyib) |
| `a7-berpa` | `ad6104e` | `ad6104e` (dəyişməyib) |
| `fix-ui-parameter-visibility` | *yox idi* | `9d3ac51` (**yeni uzaq budaq**, `main` ilə eyni commit) |

### Seans 13-ün 1-ci açıq sualı BAĞLANDI

Sual idi: "`b4-thp` → `main` birləşdirilsinmi?" Cavab artıq verilib —
digər maşında birləşdirilib:

```
git merge-base --is-ancestor b4-thp origin/main   ->  BƏLİ
git log origin/main..b4-thp                       ->  BOŞ
```

`b4-thp` tam şəkildə `origin/main`-in içindədir. Yəni B4-A (THP) artıq
`main`-dədir. Bu səbəbdən **aktiv budaq `b4-thp`-dən `main`-ə keçirildi**
— ayrıca budaqda qalmağın daha bir səbəbi yoxdur.

### Endirilən iş (`3b3e66b` → `9d3ac51`)

Altı commit-in beşi `b4-thp` tərəfindən gələn və bu maşında onsuz da
mövcud olan işdir (Seans 12–13). **Həqiqətən YENİ olan tək commit:**

| Commit | Nə gətirdi |
|---|---|
| `8895993` | **Seans 14** — `fix(ui)`: işləməyən panel sahələri bozarılır + geologiya köhnəlmə banneri |
| `9d3ac51` | birləşdirmə commit-i (`b4-thp` → `fix-ui-parameter-visibility`) |

`8895993` — **5 fayl, +459 sətir, 0 silinmə**:
`ui/panels.py`, `ui/main_window.py`, `tests/test_parameter_sensitivity.py`
(YENİ, 14 test), `ISH_HESABATI.md`, `QARARLAR.md` (Q-12).

Mühərrikə toxunulmayıb — dəyişən yalnız `ui/` altındadır.

### Bu maşında edilən əməliyyat

```
git fetch --all --prune
git checkout main
git merge --ff-only origin/main     # main: 3b3e66b -> 9d3ac51
```

Yalnız fast-forward. Heç bir birləşdirmə münaqişəsi, heç bir yerli
dəyişiklik itkisi olmadı (iş ağacı əvvəlcədən təmiz idi). Hər üç yerli
budaq indi `origin` ilə eynidir.

### Mühit yoxlanışı

```
Python 3.14.7 (sistem interpretatoru, .venv YOXDUR)
PyQt5 OK · vtk 9.7.0 · numpy 2.5.2 · scipy 1.18.1 · matplotlib OK · pytest OK
import app  -> OK,  imex2d.version.VERSION = 69
```

### Bütöv test dəsti BU maşında qaçırıldı

```
python -m pytest -q
2322 keçdi · 1 skip · 1 xfail · 0 XƏTA   (742 s)
```

Seans 13-də bu maşında 2 308 idi; fərq **tam olaraq 14 testdir** —
məhz Seans 14-ün gətirdiyi `test_parameter_sensitivity.py`. Rəqəm
digər maşının nəticəsi (2 322) ilə **bitə-bit uyğundur**.

### Sənədlərə edilən dəyişiklik

Yalnız bu bölmə. `ROADMAP.md`, `QARARLAR.md`, `ARCHITECTURE.md`
**qəsdən toxunulmadı**: bu seansda nə mərhələ statusu dəyişdi, nə
texniki qərar verildi, nə də struktur. Seans 14-ün öz qərarı (Q-12)
onsuz da endirilən commit-lə birlikdə gəldi.

### Açıq suallar ⏳

1. **`.venv` bu maşında qurulsunmu?** (Seans 13-dən qalır, hələ açıq.)
   Hazırda sistem Python-u hər şeyi daşıyır və `python app.py` işləyir,
   lakin `run.bat` üç yerdə venv axtarır (Q-10), tapmayanda
   `XETA: venv tapilmadi` verib dayanır — yəni proqram bu maşında
   `run.bat` ilə AÇILMIR. İki yol var: (a) burada `.venv` qurmaq,
   (b) `run.bat`-a "venv yoxdursa sistem Python-u ilə davam et" qolu
   əlavə etmək. Seçim sahibkarındır — öz təşəbbüsümlə edilmədi.
2. **Proqram bu maşında GUI ilə açılıb sınanmadı** (Seans 13-dən qalır)
   — yalnız idxal və testlər yoxlanıldı. Seans 14 məhz UI dəyişikliyidir,
   ona görə gözlə yoxlanması xüsusilə mənalıdır.
3. **`fix-ui-parameter-visibility` uzaq budağı silinsinmi?** `main` ilə
   eyni commit-dədir, yəni işi bitib.

---

## 12 sentyabr 2026 — Seans 15 (davamı): `.venv` quruldu, proqram gözlə yoxlandı

Yuxarıdakı Seans 15 bölməsi iki açıq sual qoymuşdu. Sahibkar hər
ikisinə cavab verdi: **(1) `.venv` burada qurulsun**, **(2) proqram
açılıb gözlə yoxlansın**. Bu bölmə onların icrasıdır.

### 1 · `.venv` — adi yol UĞURSUZ oldu

`python -m venv .venv` + `pip install -r requirements-dev.txt`
işlədi (22 paket, `pip check` təmiz), amma proqram AÇILMADI:

```
ImportError: DLL load failed while importing vtkCommonCore:
An Application Control policy has blocked this file.
```

Əvvəl eyni xəta `scipy._ufuncs_cxx`-də çıxdı, bir neçə dəqiqə sonra
scipy keçdi, **vtk isə bloklanmış qaldı**. Yəni Windows Application
Control (ROADMAP 0.7-dəki Smart App Control maneəsinin eyni ailəsi)
**təzə endirilmiş imzasız DLL-ləri** bloklayır — paketin özünü yox.

Təsdiq: eyni paketlərin **sistem Python-undakı** nüsxələri işləyir
(bu maşında 2 322 test keçir, VTK-nın 49 testi daxil).

### 1b · Həll — nazik venv

Sistem Python-u tələb olunan **18 paketin hamısını dəqiq pinlənmiş
versiyalarda** daşıyır (`requirements.txt` + `requirements-dev.txt`
ilə bir-bir tutuşduruldu: **uyğunsuzluq 0**). Ona görə:

```
python -m venv --system-site-packages .venv
```

`run.bat` artıq venv tapır, idxallar isə etibarlı sayılan sistem
DLL-lərinə düşür. Tam əsaslandırma: [QARARLAR.md](QARARLAR.md) → **Q-13**.

⚠️ **Bundan sonra yeni asılılıq SİSTEM Python-una qurulmalıdır** —
`.venv`-ə pip ilə qurulan hər şey yenidən bloklana bilər.

### 1c · Nəticə — `run.bat` işləyir

Proqram açıldı: `IMEX-2D v69 · 12 tab: Layihə, Model, Nəticələr,
Nisbi keçiricilik, 3D görüntü, PVT, Validasiya (B-L), Müqayisə,
Tarixçə, Uyğunlaşdırma, Həssaslıq, Jurnal`.
Seans 13-ün 2-ci açıq sualı **BAĞLANDI**.

### 2 · Seans 14-ün UI düzəlişi GÖZLƏ təsdiqləndi

**Mexanizm A (bozarma) — canlı proqramda görüldü.** PVT korrelyasiya
rejimində (`pvt: correlation(API=32, γg=0.75, T=70°C)`) süxur/flüid
panelində **Su lözlüyü μw, Neft lözlüyü μo, Lözlük vahidi, Bo**
bozarmış vəziyyətdədir, altında narıncı izah:

> Lözlük və Bo PVT tabından gəlir — buradakı dəyərlər işlədilmir.
> Onları dəyişmək üçün PVT tabındakı API, temperatur və doyma
> təzyiqini redaktə edin.

**Mexanizm B (köhnəlmə banneri)** — kor-kora klikləmə akkordeonu
sürüşdürdüyü üçün etibarsız oldu; əvəzinə real `MainWindow` qurub
widget-ləri birbaşa idarə edən sürücü skript yazıldı. Altı quyu
(φ 0.16–0.28, k 60–310 mD) qurulub **"İnterpolyasiya et" BASILDI** —
`GeologicalModel` həqiqətən quruldu (layihə ağacında "Quyu
cədvəlindən geoloji model", K-1…K-6):

| Addım | Banner | Gözlənilən |
|---|---|---|
| interpolyasiyadan sonra | *(boş)* | təzə ✅ |
| neft lözlüyü μo dəyişdi | *(boş)* | geologiyaya təsir etmir ✅ |
| **məsaməlilik φ dəyişdi** | "Nəticə köhnəlib — 'İnterpolyasiya et' basın." | köhnəlir ✅ |
| yenidən interpolyasiya | *(boş)* | təzələnir ✅ |
| **grid NX dəyişdi** | "Nəticə köhnəlib — ..." | köhnəlir ✅ (Q-12 qərar 3) |

Yəni Q-12-nin hər dörd qərarı canlı proqramda işləyir.

### Yol boyu tapılan, düzəliş TƏLƏB ETMƏYƏN davranışlar

* **Boş quyu cədvəli ilə interpolyasiya haqlı olaraq RƏDD EDİLİR**:
  "Heç bir xassə üçün seçilmiş üsula kifayət qədər quyu yoxdur
  (tələb olunan: 3)". İlk sürücü cəhdim məhz buna düşdü — proqram
  düzgün davrandı, **test məlumatım yanlış idi**.
* Konsola çıxışda Azərbaycan hərfləri `cp1252`-də sınır
  (`UnicodeEncodeError`). Bu YALNIZ skriptlə sürüklənəndə görünür,
  proqramın öz jurnalına (`logs/imex2d.log`) aid deyil. Sürücü
  `PYTHONUTF8=1` ilə işlədildi; **koda toxunulmadı**.

### Bu bölmədə kod dəyişmədi

`imex2d/` altında bir sətir də dəyişməyib. Dəyişən: `.venv` (repoya
düşmür, `.gitignore`-da) və bu sənədlər.

### Açıq suallar ⏳

1. **`fix-ui-parameter-visibility` uzaq budağı silinsinmi?** `main`
   ilə eyni commit-dədir (Seans 15-dən qalır).
2. **Yeni asılılıq lazım olanda hara qurulsun?** Q-13-ə görə sistem
   Python-una — amma bu, maşına qlobal təsir edir. Sahibkar başqa
   yol istəyirsə (məsələn Application Control siyasətinin
   dəyişdirilməsi), qərar onundur.

---

## 12 sentyabr 2026 — Seans 16: B5-a — nəticələrin CSV / JSON ixracı

### Nə üçün

Layihədə yalnız **PDF hesabat** vardı (`reporting/report.py`) — insana
baxmaq üçün. B4-A `well_bhp` və `well_thp` sıralarını yaratdıqdan sonra
istifadəçinin əlində qrafikdən başqa heç nə yox idi: rəqəmləri Excel-ə
və ya başqa alətə çıxarmaq mümkün deyildi.

Sahibkarın qərarı: **B5-a əvvəl, B5-b (Pcog) ayrıca commit-də** — Pcog
mühərriyə toxunur, B5-a isə yalnız oxuyur.

### Nə edildi

| Fayl | Nə |
|---|---|
| `reporting/results_export.py` | **YENİ** — `write_csv`, `read_csv`, `write_json`, `write` |
| `ui/main_window.py` | menyuya "Nəticələri ixrac et (CSV/JSON)…" + `export_results()` |
| `tests/test_results_export.py` | **YENİ** — 14 test |

İki format, iki məqsəd:

* **CSV** — TƏMİZ cədvəl (başlıq + sətirlər), şərh sətri yoxdur. Excel və
  `pandas.read_csv` heç bir parametr olmadan açır.
* **JSON** — eyni məlumat + metadata (model adı, OOIP/OGIP, addım sayı,
  yığılma vəziyyəti).

Metadata CSV-yə QƏSDƏN salınmadı: `#` şərh sətirləri ciddi CSV
oxuyucularını sındırır və gedər-gələr müqaviləsini pozardı.

### Ölçülmüş nəticə

8×8 5-spot, 400 gün, qaz aktiv, lülə verilmiş:

```
CSV : 17 sütun, 30 sətir, 9.4 KB
JSON: 13.9 KB, metadata + 12 seriya sütunu + 1 quyu
```

### Testin tutmadığı, GÖZLƏ görünən qüsur

İlk versiyada vahid **vergüllə** yazılırdı (`"t, gün"`). Gedər-gələr
testi keçirdi, çünki `csv` modulu dırnaqlanmış sahəni düzgün oxuyur.
LAKİN faylı açıb baxanda göründü: **hər başlıq dırnağa düşür**, çünki
vergül CSV-nin öz ayırıcısıdır. `awk`/`cut` kimi sadə alətlər belə
faylı sındırır.

Düzəliş: vahid **kvadrat mötərizədə** — `t [gün]`, `RF [%]`,
`GOR [sm³/sm³]`. Başlıqlar artıq dırnaqsızdır. Test də yeniləndi və
indi açıq şəkildə `"," not in header` yoxlayır.

> Qeyd: bu qüsuru test tapmadı — **faylı açıb baxmaq** tapdı.

### Dörd tələ (hamısı testlə kilidləndi)

1. **`nan` BOŞ xana kimi yazılır, `0` kimi YOX.** `well_thp`-də quyunun
   səthə axa bilmədiyi addımlar `nan`-dır (B4-A, qəsdən belədir).
   `0 bar` yazsaydıq, qrafikdəki kimi real ölçmə kimi oxunardı.
   JSON-da isə `null` — `NaN` JSON spesifikasiyasına ziddir.
2. **İki fazalı nəticədə qaz sıraları boşdur** → sütun ümumiyyətlə
   yaradılmır. Beləcə faylın özü "burada qaz fazası yoxdur" deyir.
3. **`well_bhp`/`well_thp` yalnız bəzi quyularda var** (BHP rejimli,
   lüləsi verilmiş istismarçılar) — `.get(name, [])` naxışı.
4. **Vahid hər başlıqdadır** (`UNITS.md` prinsipi).

### Öz təşəbbüsümlə verilmiş qərarlar

1. **CSV `utf-8-sig` (BOM ilə) yazılır.** Windows Excel BOM-suz UTF-8-i
   tanımır və `ə/ş/ğ` hərfləri korlanır. `csv`/`pandas` BOM-u özləri
   atır, yəni gedər-gələr pozulmur. Ayrıca testlə kilidləndi.
2. **Rəqəmlər `repr(float)` ilə, tam dəqiqliklə yazılır** (məs.
   `8.952477630417802e-05`). Yuvarlaqlaşdırma gedər-gələr müqaviləsini
   pozardı; ixrac faylı hesabat deyil, MƏLUMATDIR.
3. **`read_csv()` modulun özündə saxlanıldı**, testə yazılmadı —
   formatın müqaviləsi ("geri oxunanda eyni ədədlər") koddan
   görünməlidir, yalnız testdən yox.
4. **Qısa sütun uydurma dəyərlə doldurulmur**, boş qalır.

### Qalan ⏳

* **B5-b — Pcog** (qaz-neft kapilyar təzyiqi + `dpcog_dsg`). Mühərriyə
  toxunur, ayrıca commit-də ediləcək.

---

## 12 sentyabr 2026 — Seans 17: `max_dt` üç fazalı mühərrikdə hörmət olunmurdu

### Bildiriş

Sahibkar: interfeysdə "Maks. Δt = 20 gün" qoyulur, status sətri isə
orta addımı 26.3 gün göstərir.

### Diaqnoz — səbəb CFL DEYİL

Sahibkarın fərziyyəsi "CFL tənzimləməsi limiti qulaqardına vurur" idi.
**Ölçmə başqa şey göstərdi:** limitləmə məntiqi tamamilə düzgündür —
`time_stepping.py:126` `dt = min(self.dt, config.max_dt, remaining)`,
`:221` isə `np.clip(..., min_dt, max_dt)` edir.

Problem odur ki, **`max_dt`-nin ÖZÜ şişirdilir**:

```python
# three_phase_engine.py:110
max_dt=max(stepping.max_dt, 30.0)      # ← süni döşəmə
```

Yəni istifadəçi 20 desə də, həlledici 30 alır və limitləmə düzgün
işləyərək 30-a qədər addım atır.

### ⚠️ Bu səhv BİR DƏFƏ ARTIQ DÜZƏLDİLİB

`engine.py::_time_config` sənədi hərfən yazır:

> "`max_dt` istifadəçinin sorduğu kimi hörmət edilir. ƏVVƏLLƏR
> (TAPILAN SƏHV) bura süni minimum (30 gün) tətbiq olunurdu —
> istifadəçi 0.5 və ya 2 gün desə də, mühərrik səssizcə 30 günə
> keçirdi."

İki fazalı mühərrikdə düzəldilmiş, üç fazalıda qalmışdı. Səbəb aydındır:
üç fazalı yol **A7 bərpasında (B2) git tarixçəsindən qaytarılıb** və
köhnə kodu özü ilə geri gətirib.

### Ölçülmüş — əvvəl / sonra

8×8 5-spot, 400 gün:

| Mühərrik | İstənilən | ƏVVƏL işlədilən | SONRA | Maks Δt (sonra) |
|---|---|---|---|---|
| iki fazalı | 5 / 20 / 50 | 5 / 20 / 50 ✅ | eyni | 5.00 / 20.00 / 50.00 |
| üç fazalı | 5 | **30** ❌ | **5** ✅ | 5.00 |
| üç fazalı | **20** | **30** ❌ | **20** ✅ | 20.00 |
| üç fazalı | 50 | 50 ✅ | 50 | 50.00 |

Simptom kəskin idi: **5 və 20 TAM EYNİ nəticə verirdi** (31 addım,
maks Δt 30.00), çünki hər ikisi eyni həddə qaldırılırdı.

### Düzəlişin doğruluğunun ƏN GÜCLÜ sübutu

Pb = 100 bar-da qaz ayrılmır (maks S_g = 0.000000), yəni üç fazalı
mühərrik faktiki olaraq iki fazalı məsələni həll edir — nəticələr
ÜST-ÜSTƏ DÜŞMƏLİDİR:

| | RF | Addım |
|---|---|---|
| iki fazalı | 62.861505 | 31 |
| **üç fazalı (düzəlişdən sonra)** | **62.861429** | 31 |
| fərq | **7.6·10⁻⁵** | — |
| üç fazalı (düzəlişdən əvvəl) | 62.719661 | — |

Düzəlişdən əvvəl bu uyğunluq YOX İDİ — 30 günlük addım nəticəni
kobudlaşdırırdı.

### Bir ETALON DƏYİŞDİ

`test_three_phase_high_bubble_point.py::test_low_bubble_point_results_are_unchanged`
Pb = 100 üçün **RF 62.72** gözləyirdi. Bu dəyər Seans 8-də ölçülmüşdü —
yəni **səhvin öz dəyərini kilidləyirdi**. Etalon 62.86-ya yeniləndi və
səbəbi testin sənədinə yazıldı.

> Bu, layihənin "etalon pozulursa dəyişiklik səhvdir" qaydasının
> İSTİSNASIDIR: burada pozulan şey səhvin özü idi. Ona görə etalon
> dəyişdirildi və əvəzinə DAHA GÜCLÜ invariant qoyuldu — iki və üç
> fazalı mühərriklərin qazsız rejimdə uyğunluğu.

### Toxunulan fayllar

| Fayl | Nə |
|---|---|
| `simulation/implicit/three_phase_engine.py` | `max(stepping.max_dt, 30.0)` → `stepping.max_dt` + sənəd |
| `tests/test_max_timestep_respected.py` | **YENİ** — 16 test, hər iki mühərrik üçün kilid |
| `tests/test_three_phase_high_bubble_point.py` | etalon 62.72 → 62.86, səbəbi sənədləşdi |

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Testlər hər İKİ mühərrik üçün yazıldı**, təkcə üç fazalı üçün yox.
   Səbəb: səhv əvvəl iki fazalıda düzəldilib, sonra üç fazalıda geri
   qayıdıb. Ortaq kilid üçüncü dəfəni dayandırır.
2. **`soft_failure_*` fərqinə TOXUNULMADI.** İki fazalı `_time_config`
   onları verir, üç fazalı vermir. Bu, AYRI fərqdir və yığılma
   davranışına toxunur — ⏳ kimi qeyd olundu.
3. **Etalon dəyişdirildi, testə `xfail` qoyulmadı.** Yeni dəyər fiziki
   olaraq daha düzgündür və iki fazalı mühərriklə uyğunluqla
   təsdiqlənir.

### Açıq qalan ⏳

* Üç fazalı `_time_config`-də `soft_failure_cnv_tolerance` və
  `soft_failure_mb_tolerance` yoxdur (iki fazalıda var). Onları əlavə
  etmək yığılma davranışını dəyişər — ayrıca ölçülməlidir.

---

## 12 sentyabr 2026 — Seans 18: B6-b — oynatma idarəsi

### Planın B6 bölməsi vəziyyəti OLDUĞUNDAN ZƏİF göstərirdi

İşə başlamazdan əvvəl dörd bəndin hamısı yoxlanıldı:

| Plan bəndi | Real vəziyyət |
|---|---|
| 2. Cəbhə animasiyası (`QTimer`) | 🟢 **əsasən HAZIR İDİ** |
| 4. Dashboard-a THP | ✅ B4-A-da edilmişdi |
| 1. İnteraktiv kəsik | ❌ yoxdur |
| 3. GIF / PNG ixracı | ❌ yoxdur |

`main_window`-da `▶ Oynat` düyməsi, zaman slider-i və 140 ms-lik
`QTimer` ARTIQ vardı. Slider həm 2D xəritəni, həm də **3D VTK
görüntüsünü** yeniləyir (`_update_volume_vtk` sənədi: *"səhnə
keşlənir, zaman slider-i sürüşdürəndə yalnız DƏYƏRLƏR yenilənir"*).

Yəni cəbhənin hərəkəti onsuz da izlənilirdi. Çatışmayan: sürət
tənzimi, kadr-kadr addımlama, sonda dayanma seçimi.

### Nə edildi

| Fayl | Nə |
|---|---|
| `ui/playback.py` | **YENİ** — saf məntiq: `interval_ms`, `step_value`, `advance` |
| `ui/main_window.py` | ◀ / ▶ düymələri, sürət seçicisi (0.25×…4×), "Dövrə" qutusu |
| `tests/test_playback_controls.py` | **YENİ** — 20 test |

Yeni idarə elementləri:

* **◀ / ▶** — kadr-kadr geri/irəli. Oynatma gedirsə **dayandırılır**
  (yoxsa taymer istifadəçinin seçdiyi kadrı dərhal üstələyərdi).
* **Sürət** — 0.25× / 0.5× / **1×** / 2× / 4×. Əmsal intervalı bölür.
* **Dövrə** — işarəli (defolt): sonda əvvələ qayıdır; söndürülmüş:
  son kadrda dayanır.

### TESTİN QAÇIRICINI ÇÖKDÜRMƏSİ — dizaynı dəyişdirdi

İlk versiyada testlər `MainWindow()` qururdu. Nəticə: **`pytest` fatal
xəta ilə çökdü** (C stack trace, 14 test birdən). Səbəb: pəncərənin
qurulması VTK səhnəsi və bir neçə matplotlib kanvası yaradır.

Diqqətəlayiq: **mövcud 2350+ testin heç biri `MainWindow` qurmur** —
hamısı panelləri ayrıca sınayır. Bu, təsadüf deyilmiş.

Ona görə oynatma məntiqi `ui/playback.py`-yə çıxarıldı: Qt-dən tam
asılısız, saf funksiyalar. `main_window` yalnız taymeri idarə edir və
onları çağırır. Bu, layihənin öz qaydası ilə də uyğundur
(`ARCHITECTURE.md` → 1.3 "UI biznes məntiqini daşımamalıdır").

### Öz təşəbbüsümlə verilmiş qərarlar

1. **Kadr-kadr addımlama DÖVRƏ ETMİR**, avtomatik oynatma edir.
   Səbəb: əl ilə addımlayan istifadəçi son kadrdan birdən başlanğıca
   atılmağı gözləmir; avtomatik oynatmada isə dövrə faydalıdır.
2. **Defolt dəyərlər köhnə davranışı BİTƏ-BİT saxlayır** — 1× = tam
   140 ms, dövrə açıq. Yəni B6-b-yə toxunmayan istifadəçi heç bir
   fərq görmür. Ayrıca testlə kilidləndi.
3. **İntervalın aşağı həddi 10 ms** qoyuldu. Çox yüksək sürətdə Qt
   hadisə növbəsi boğulur və interfeys cavab verməz olur.
4. **Yararsız sürət dəyəri 1×-ə düşür** (`None`, mənfi, mətn) —
   taymer heç bir halda partlamamalıdır.

### Qalan ⏳

* **B6-a** — interaktiv kəsik müstəvisi (`vtkPlaneWidget` + `vtkCutter`)
* **B6-c** — GIF / PNG ixracı (Pillow 12.3.0 hazırdır)

⚠️ B6-a üçün bilinməli: `vtkPlaneWidget` **interaktor tələb edir**,
ekransız testdə mümkün olmaya bilər. Ölçülməlidir; mümkün olmasa,
kəsik məntiqi (`vtkCutter`) widget-dən ayrı test olunmalıdır.

---

## 13 sentyabr 2026 — Seans 19: GitHub-dakı son dəyişikliklər (Seans 16–18) bu nüsxəyə endirildi

### 1 · Başlanğıc vəziyyət

`git fetch --all --prune` nəticəsi:

| Budaq | Vəziyyət |
|---|---|
| `origin/main` | `8a85e64` — lokal `main` ilə **eyni** (0 / 0) |
| `origin/a7-berpa` | **GitHub-da silinib** (`[gone]`) |
| `origin/b4-thp` | **GitHub-da silinib** — `main`-in içindədir ✅ |

İş bu nüsxədə `a7-berpa` budağında gedirdi — `main`-dən **5 commit geridə**.

### 2 · GitHub-dan gələn 5 commit

| Commit | Tarix | Nə |
|---|---|---|
| `2448ab9` | 12.09 15:29 | chore: `B5_PROMPT.md` `.gitignore`-a |
| `e3cd8a3` | 12.09 15:48 | **B5-a** — nəticələrin CSV / JSON ixracı (Seans 16) |
| `6f4827d` | 12.09 15:48 | chore: `NUMUNE_NETICELER.*` `.gitignore`-a |
| `f63d361` | 12.09 16:12 | **fix** — `max_dt` üç fazalı mühərrikdə hörmət olunmurdu (Seans 17) |
| `8a85e64` | 12.09 16:37 | **B6-b** — oynatma idarəsi: sürət, kadr-kadr, dövrə (Seans 18) |

Cəmi: **13 fayl, +1 345 / −121**.

| Sahə | Fayllar |
|---|---|
| **B5-a ixrac** | `reporting/results_export.py` (**yeni**, 216) · `ui/main_window.py` (menyu) · `tests/test_results_export.py` (**yeni**) |
| **`max_dt` düzəlişi** | `simulation/implicit/three_phase_engine.py` (`max(max_dt, 30.0)` → `max_dt`) · `tests/test_max_timestep_respected.py` (**yeni**) · `tests/test_three_phase_high_bubble_point.py` (etalon 62.72 → 62.86) |
| **B6-b oynatma** | `ui/playback.py` (**yeni**, saf məntiq) · `ui/main_window.py` (◀ / ▶, 0.25×…4×, "Dövrə") · `tests/test_playback_controls.py` (**yeni**) |
| **Sənədlər** | `ISH_HESABATI.md` (Seans 16–18) · `QARARLAR.md` (**Q-14** — ixrac formatı) · `ROADMAP.md` (5.6 → ✅) · `ICRA_PLANI.md` (M8 ✅, B5-a ✅, B6-b ✅) |

Təfərrüat həmin seansların öz bölmələrindədir (yuxarıda) — burada təkrarlanmır.

### 3 · Yoxlama — gələn testlər bu nüsxədə qaçırıldı

`main` (`8a85e64`) üzərində, `.venv` (Python 3.12.10, pytest 9.1.1):

```
tests/test_results_export.py                  14 keçdi
tests/test_max_timestep_respected.py          16 keçdi
tests/test_playback_controls.py               20 keçdi
tests/test_three_phase_high_bubble_point.py   23 keçdi
============== 73 passed in 22.67s ==============
```

**Bütöv test dəsti QAÇIRILMADI** — yalnız dəyişikliklə gələn 4 fayl.

### 4 · TAPINTI — iki fərqli "Seans 16" var, biri GitHub-da yoxdur

Bu nüsxədə `a7-berpa` budağında GitHub-a çatmamış bir commit qalıb:

```
b2a795e  12.09 21:07  docs: Seans 16 - GitHub deyishiklikleri bu nusxeye endirildi
```

O, `ISH_HESABATI.md`-yə **"Seans 16: GitHub-dakı dəyişikliklər bu nüsxəyə
endirildi"** bölməsini (+147 sətir) əlavə edir. Eyni vaxtda başqa maşında
`main`-ə **"Seans 16: B5-a"** yazılıb. Nəticə:

* seans nömrələri toqquşur (iki Seans 16);
* `origin/a7-berpa` silindiyi üçün `b2a795e` **yalnız bu diskdə** yaşayır —
  nüsxə itsə, qeyd də itir (CLAUDE.md: "tarixçə itmir").

Bu seansda `a7-berpa`-nı `main`-ə birləşdirmək cəhdi **icazə sistemi
tərəfindən bloklandı**. Ona görə commit-ə **toxunulmadı**: lokal
`a7-berpa` budağı olduğu kimi qalır, bu bölmə isə ondan asılı deyil.

### 5 · Yol boyu müşahidə — proqram köhnə kodla işlədildi

Bu seansdan əvvəl (12.09 21:12 – 22:05) proqram `a7-berpa` üzərində
açılmışdı — yəni **B5-a, B6-b və `max_dt` düzəlişi OLMADAN**. Jurnalda
4 hesablama var (RUN-001…004; TPFA, üç fazalı, MPFA-O), xəta yoxdur.
Proses 127 kodu ilə bağlandı, traceback yoxdur — səbəb **müəyyən edilmədi**.

Qeyd: RUN-002/003 üç fazalı mühərriklə getdi və `max_dt` düzəlişindən
**əvvəlki** koddur — həmin nəticələr (RF 1.53 %) Seans 17-dəki səhvdən
təsirlənmiş ola bilər. ⏳ ölçülməyib.

### Bu bölmədə nə dəyişmədi və niyə

* `imex2d/` altında bir sətir də dəyişməyib — yalnız çəkmə, test, qeyd.
* `QARARLAR.md` — yeni texniki seçim yoxdur.
* `ROADMAP.md`, `ARCHITECTURE.md`, `ICRA_PLANI.md` — `main`-dən gələn
  halda artıq aktualdır.

### Açıq suallar ⏳

1. **`b2a795e` (lokal Seans 16) nə olsun?** Variantlar: `main`-ə
   birləşdirmək (başlığı dəyişmədən, toqquşma bu bölmədə izah olunub)
   və ya yalnız lokal saxlamaq.
2. **Bütöv test dəsti** hələ də bu nüsxədə qaçırılmayıb.
3. **Seans nömrələməsi** iki maşında paralel işdə toqquşur — hər
   seansdan əvvəl `git fetch` qaydası CLAUDE.md-yə əlavə olunsunmu?
4. Seans 16-dan qalan: **bu nüsxə hansı maşındır?**

---

## 13 sentyabr 2026 — Seans 20: git işləri + Sg xassəsi görüntüyə əlavə olundu

### 1 · Git

* `git pull --ff-only` → lokal `main` `8a85e64` → **`3e357bb`** (ikinci
  nüsxənin Seans 19 qeydi gəldi, kod dəyişikliyi yoxdur).
* **`b2a795e` XİLAS EDİLƏ BİLMƏDİ — bu kompüterdə yoxdur.** Yoxlanıldı:
  `git cat-file -t b2a795e` → "Not a valid object name"; reflog-da da yoxdur;
  `C:\Dev\LSM` bu kompüterdə mövcud deyil. Seans 19-a görə commit ikinci
  maşının diskindədir (orada `.venv` / Python 3.12, burada `venv` / 3.14).
  Xilas yalnız HƏMİN maşından mümkündür — ⏳ sahibkara əmrlər verildi.
* Seans 19 qaçışlarının (RF 1.53 %) yenidən ölçülməsi — sahibkarın qərarı
  ilə **təxirə salındı**.

### 2 · Bildiriş

Sahibkar üç fazalı mühərriki işlətdi: status sətrində GOR görünür, amma
"Xassə" siyahısında **Qaz doyumluluğu (Sg) yox idi** — qazın cəbhəsini
ekranda izləmək mümkün deyildi.

### 3 · Səbəb

Məlumat ARTIQ VAR idi: `three_phase_engine._record_snapshot` hər snapshot-a
`gas_saturation` yazır. Çatışmayan yalnız görüntü qatı idi — iki açılan
siyahıda (Model tabı, 3D tabı) və renderer seçicilərində Sg açarı yox idi.

### 4 · GİZLİ TƏLƏ — sadəcə siyahıya əlavə etmək YANLIŞ olardı

Hər iki seçicinin (`_select_volume`, `_select`) sonunda tanınmayan açar
üçün **səssiz defolt** var:

```python
return model.rock.porosity.values.reshape(shape3d), POROSITY_CMAP, None, None
```

Yəni "Sg" yalnız siyahıya salınsaydı, ekranda **Sg adı altında məsaməlilik**
görünərdi — heç bir xəta vermədən. Ona görə hər iki seçiciyə AÇIQ budaq
yazıldı və bu, ayrıca testlə kilidləndi.

### 5 · Edilən

| Fayl | Nə |
|---|---|
| `rendering/theme.py` | `GAS_SATURATION_CMAP` — su/neft xəritəsindən fərqli rəng ailəsi |
| `rendering/renderers.py` | `GAS_SATURATION = "SG"`, etiket, `_gas_saturation()`, `_gas_saturation_limits()`, iki seçici budağı |
| `ui/main_window.py` | Sg **həm Model, həm 3D** tabının siyahısına |
| `tests/test_gas_saturation_view.py` | **YENİ** — 10 test |

### 6 · Öz təşəbbüsümlə verilmiş qərarlar

1. **Sg hər İKİ taba əlavə olundu**, təkcə 3D-yə yox. Səbəb: B6-b-nin
   oynatma düymələri **Model** tabındadır — qazın cəbhəsini kadr-kadr
   izləmək orada mümkündür.
2. **Snapshot yoxdursa `NaN`**, sıfır yox. Simulyasiyadan əvvəl ilkin qaz
   papağı ola bilər; renderer onu bilmir — uydurma 0 verilmir.
3. **İki fazalı snapshot-da sıfırlar.** Modeldə qaz fazası yoxdur, Sg ≡ 0
   fiziki olaraq DƏQİQDİR.
4. **Rəng şkalası sabitdir: `(0, 1 − Swc)`.** Avtomatik şkala hər kadrda
   dəyişərdi: erkən kadrda cüzi Sg tünd qırmızı görünər, cəbhənin böyüməsi
   gözə çarpmazdı. `1 − Swc` Sg-nin fiziki yuxarı həddidir.

### 7 · Yoxlama

```
tests/test_gas_saturation_view.py + görüntü testləri   127 keçdi
```

Real üç fazalı qaçışda (Pb = 240 bar) son kadrda Sg > 0 hüceyrələr görünür
və qaz zonası zamanla kiçilmir — testlə kilidləndi.

### Açıq qalan ⏳

* **`b2a795e`** — ikinci maşından xilas edilməlidir.
* **B6-a** — interaktiv kəsik (növbəti).

---

## 13 sentyabr 2026 — Seans 21: B6-a — interaktiv kəsik müstəvisi

### 1 · Tətbiqdən ƏVVƏL ölçülən risk

B6-a planında əsas risk qeyd olunmuşdu: `vtkPlaneWidget` interaktor tələb
edir, ekransız test mühitində işləməyə bilər. Kod yazmadan ölçüldü:

| Yoxlama | Nəticə |
|---|---|
| VTK versiyası | 9.7.0 |
| ekransız `Render()` | işləyir |
| `vtkImplicitPlaneWidget2` ekransız `On()` | **işləyir** |
| `vtkCutter` hüceyrə skalyarlarını daşıyır | bəli — eyni rəng cədvəli işlədilə bilər |
| `vtkCutter` gizlədilmiş (blank) hüceyrələrə hörmət edir | bəli — K-filtri ilə 18 → 6 hüceyrə |

Nəticə: widget-i məntiqdən ayırmağa **ehtiyac olmadı**, birbaşa sınanır.
Kəsim həddi, K aralığı və status filtri kəsiyə **avtomatik** tətbiq olunur.

### 2 · Edilən

| Fayl | Nə |
|---|---|
| `rendering/vtk_volume.py` | `VtkViewSettings.slice_axis` / `slice_position`; saf `slice_plane()` və `slice_fraction()`; `update_slice()`, `slice_output()`, `attach_slice_widget()` |
| `ui/main_window.py` | 3D tabında "Kəsik: Yox / X / Y / Z — dərinlik" + mövqe sürgüsü; sürüklənən müstəvi sürgü ilə sinxron |
| `tests/test_vtk_slice.py` | **YENİ** — 28 test |

Davranış:

* Kəsik aktiv olanda əsas həcm gizlədilir, kəsik səthi görünür (quyular və
  faylar qalır). Kəsik söndürüləndə həcm geri qayıdır.
* Müstəvi iki yolla hərəkət edir: **sürgü** ilə və ya 3D görüntüdə
  **birbaşa sürükləyərək**. Sürükləmə sürgünü yeniləyir; proqramla qoyulan
  mövqe widget hadisəsini təkrar yaratmır — dövrə yoxdur.
* Normal seçilmiş oxa kilidlənir — müstəvi yalnız ox boyunca sürüşür.
* Kəsik **yalnız VTK motorunda** — matplotlib seçiləndə idarə söndürülür.
* Z oxunda koordinat `−dərinlik` olduğu üçün sürgüdə 0 = **ən dərin** təbəqə.

### 3 · İKİ YANLIŞ FƏRZİYYƏ — ölçmə ilə təkzib olundu

**Fərziyyə 1 — "müstəvi hüceyrə sərhədinə düşəndə üst-üstə poliqon yaranır".**
İlk probda 8×6×3 modelin X-ortasında 18 əvəzinə 36 poliqon görüldü və bu,
iki qonşu hüceyrənin ortaq üzü (z-fighting) kimi yozuldu. Buna görə müstəvini
10⁻⁶ qədər sürüşdürən `_SLICE_NUDGE` yazıldı.

Testlər düşdü: sürüşdürmədən sonra da 36 idi. Ayrıca ölçüldü:

| Müstəvinin yeri | poliqon | tip | **unikal hüceyrə** |
|---|---|---|---|
| sərhəddə dəqiq (x = 100) | 36 | üçbucaq | **18** |
| sərhəddən 10⁻⁴ kənar | 36 | üçbucaq | **18** |
| hüceyrənin ortası | 36 | üçbucaq | **18** |

Həqiqət: `vtkCutter` **hər hüceyrəni iki üçbucaq** kimi verir. Üst-üstə
düşmə **yoxdur**. `_SLICE_NUDGE` və onun şərhi **çıxarıldı**; testlər
poliqon sayını deyil, unikal hüceyrə sayını yoxlayır.

**Fərziyyə 2 — "`SetNormalToZAxis(1)` widget-i Z-yə çevirir".** Ölçüldü:
bayraq normal **vektorunu dəyişmir** (1, 0, 0 qalır). Vektor indi açıq
`SetNormal(...)` ilə verilir, əvvəlki oxun bayrağı söndürülür.

> Dərs (Seans 16-dakı kimi): **testin keçməsi yetərli deyil, ölçmə lazımdır**
> — burada isə əksinə, testin DÜŞMƏSİ yanlış izahı üzə çıxardı.

### 4 · Öz təşəbbüsümlə verilmiş qərarlar

1. **Kəsik zamanı həcm gizlədilir.** Şəffaf həcm + kəsik birlikdə oxunmaz
   olur; quyular və faylar kontekst üçün qalır.
2. **Normal oxa kilidlidir** — maili müstəvi lay kəsiyində oxunmaz olur.
3. **Həndəsə saf funksiyalardadır** (`slice_plane`, `slice_fraction`) —
   VTK obyektsiz sınanır.

### 5 · Yoxlama

```
tests/test_vtk_slice.py + VTK + Sg görüntü testləri   154 keçdi
```

### Açıq qalan ⏳

* **B6-c** — GIF / PNG ixracı (`vtkWindowToImageFilter` + Pillow).
* **`b2a795e`** — ikinci maşından xilas edilməlidir.

---

## 13 sentyabr 2026 — Seans 22: B6-c — PNG kadr və animasiyalı GIF ixracı

### 1 · Tətbiqdən ƏVVƏL ölçülən

| Yoxlama | Nəticə |
|---|---|
| ekransız VTK pəncərəsindən kadr (`vtkWindowToImageFilter`) | işləyir — 10 kadr 0.11 san, RGB |
| Pillow 12.3 ilə animasiyalı GIF | yazılır; geri oxunanda 10 kadr, 140 ms, sonsuz dövrə |
| yeni asılılıq | **yoxdur** (`imageio` / `ffmpeg` lazım deyil) |

### 2 · Yol boyu TAPILAN SƏHV — "Şəkli saxla…" yanlış görüntünü yazırdı

3D tabındakı "Şəkli saxla…" düyməsi HƏMİŞƏ `volume_fig.savefig(...)`
çağırırdı. VTK motoru aktiv olanda həmin matplotlib fiquru ekranda deyil,
**gizli kanvasdadır** — yəni saxlanılan şəkil istifadəçinin gördüyü VTK
görüntüsü DEYİLDİ. İndi kadr aktiv motordan tutulur; matplotlib yolu olduğu
kimi qalıb (vektor PDF, dpi 200).

### 3 · Edilən

| Fayl | Nə |
|---|---|
| `rendering/animation_export.py` | **YENİ** — `capture_render_window`, `capture_figure`, `write_png`, `write_gif` (Qt-siz) |
| `rendering/vtk_volume.py` | `set_caption()` / `caption()` — səhnənin içində "t = … gün · xassə" yazısı |
| `ui/main_window.py` | "Şəkli saxla…" düzəlişi; yeni **"Animasiyanı GIF saxla…"** düyməsi |
| `tests/test_animation_export.py` | **YENİ** — 16 test |

GIF ixracı:

* hər snapshot bir kadr; 3D zaman sürgüsü hər kadra qoyulur, görüntü çəkilir
  və **aktiv motordan** tutulur — yəni xassə, kəsik (B6-a), baxış bucağı
  ekrandakı kimi saxlanılır;
* kadr müddəti Model tabındakı **oynatma sürətindən** (B6-b) gəlir — GIF
  ekrandakı oynatma ilə eyni tempdə gedir;
* irəliləmə pəncərəsində "Dayandır" var — yarımçıq ixrac da yazılır və
  status sətrində neçə kadr olduğu göstərilir;
* sonda sürgü istifadəçinin qoyduğu kadra **qaytarılır**.

### 4 · Öz təşəbbüsümlə verilmiş qərarlar

1. **Zaman yazısı VTK səhnəsinin İÇİNDƏDİR.** Pəncərədən tutulan kadr Qt
   etiketlərini daxil etmir — yazı olmasa GIF-də hansı kadrın hansı günə aid
   olduğu görünməzdi. Testlə təsdiqləndi: yazı tutulan piksellərə düşür.
2. **GIF kadr müddətinin aşağı həddi 20 ms.** Brauzerlərin çoxu daha qısa
   müddəti 100 ms kimi oynadır — yəni "daha sürətli" GIF əslində yavaş olardı.
3. **Kadr ölçüsü dəyişsə aydın xəta.** İxrac zamanı pəncərə ölçüsü dəyişərsə,
   sıçrayan GIF səssizcə yazılmır.
4. **İxrac məntiqi Qt-siz modulda.** `MainWindow()` testdə `pytest`-i
   çökdürdüyü üçün (Seans 18) kadr tutma və yazma ayrıca sınanır.

### 5 · Yoxlama

```
tests/test_animation_export.py + kəsik + VTK + Sg + oynatma   145 keçdi
```

Real üç fazalı qaçışda (6 snapshot) GIF-in kadr sayı snapshot sayına bərabərdir
və ilk/son kadr fərqlidir — testlə kilidləndi.

### Açıq qalan ⏳

* **B6 bloku BİTDİ** (B6-a, B6-b, B6-c). Növbəti: **B4-B** (`ControlMode.THP`),
  **B5-b** (Pcog — mühərriyə toxunur), **B7** (yekun doğrulama + SPE1).
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.


## 13 sentyabr 2026 — Seans 23: B4-B — quyunu THP ilə idarə etmək (`ControlMode.THP`)

Sahibkarın tapşırığı: UI-dakı lülə parametrlərindən istifadə edərək quyunu
yalnız BHP ilə deyil, **birbaşa THP ilə** idarə etmək. B5-b (Pcog) üç fazalı
mühərrikə toxunduğu üçün sonraya saxlanıldı; `b2a795e` gözləmədədir.

### 1 · Tətbiqdən ƏVVƏL tapılan iki TƏHLÜKƏLİ boşluq

| Yer | Nə olardı |
|---|---|
| `connection.mode is ControlMode.BHP` — `residual.py`, `jacobian.py`, `three_phase_residual.py`, `impes_engine.py`, `standard_well.py`, `well_state.py` | Yeni rejim bağlantıya ötürülsəydi **hamısı onu SƏSSİZCƏ RATE kimi** işlədərdi — THP ədədi (məs. 20 bar) debit hədəfinə (20 m³/gün) çevrilərdi, heç bir xəta olmadan |
| `io/eclipse_export.py` WCONPROD | BHP olmayan istismarçını `LRAT` kimi yazır — THP ədədi deck-ə debit kimi düşərdi |

Hər ikisi dizaynla bağlandı və testlə kilidləndi (bax aşağıda, Q-15).

### 2 · Dizayn — açıq (explicit) birləşmə

```
THP (istifadəçi) + son addımın debitləri ──tərs traverse──► BHP (növbəti addım)
```

* THP quyusu bağlantı səviyyəsində **adi BHP bağlantısıdır**
  (`mode = BHP`, `thp_target = THP`). Qalıq və Jakobian **toxunulmadı** —
  onlar yalnız BHP görür.
* `ThpController` hər qəbul olunmuş addımdan sonra traversi **tərsinə**
  həll edir (ikiqat bölmə: `bhp_from_thp`) və bağlantının `target`-ini
  yerində yeniləyir.
* İlk addımdan əvvəl debit məlum deyil — statik (cüzi axınlı) neft sütunu
  ilə təxmin edilir.
* Rəqsə qarşı iki qoruyucu: relaksasiya ω = 0.5 və addım başına maksimal
  BHP dəyişməsi 25 bar.
* Hər addımda **işlədilən** BHP `result.well_bhp`-yə yazılır; THP hesabatı
  (B4-A) artıq sabit hədəfdən deyil, həmin addım-addım BHP-dən hesablanır —
  yəni qrafikdə nəzarətçinin həqiqi dəqiqliyi görünür.

### 3 · Edilən

| Fayl | Nə |
|---|---|
| `simulation/wellbore/thp_control.py` | **YENİ** — `bhp_from_thp` (tərs traverse), `ThpController` |
| `domain/wells.py` | `ControlMode.THP`; THP hədəfi təzyiq kimi yoxlanılır |
| `simulation/well_model.py` | `WellConnection.thp_target`; THP quyusu `mode = BHP` ilə qurulur |
| `simulation/implicit/engine.py`, `three_phase_engine.py` | nəzarətçi quruluşda işə düşür, hər addımdan sonra `record` + `update` |
| `simulation/impes_engine.py`, `application/simulation_service.py` | IMPES + THP **açıq rədd** (istifadəçi dilində mesaj) |
| `domain/reservoir_model.py` | "ən azı bir BHP quyusu" şərti THP-ni də sayır; THP üçün lülə tələbi, yalnız istismarçı, THP ≥ lay təzyiqi xəbərdarlığı; PVT diapazonuna THP daxil |
| `simulation/wellbore/hydraulics.py` | THP quyularında addım-addım BHP işlədilir |
| `io/eclipse_export.py` | THP quyusu ixracda **aydın xəta** verir (VFPPROD yoxdur) |
| `ui/panels.py` | Rejim siyahısına **THP**; izah sətri |
| `tests/test_thp_control.py` | **YENİ** — 21 test |

### 4 · Real 8×8 qaçışda ölçülən (5-nöqtə, lülə Ø 62 mm, 20 seqment)

| Qaçış | Yığılma | Addım | BHP ilk → son, bar | ΔBHP maks, bar | İşarə dəyişməsi | \|THP − hədəf\| median / maks, bar |
|---|---|---|---|---|---|---|
| 2 faza, THP = 20, 1000 gün | ✅ | 66 | 86.4 → 142.7 | 16.5 | 2 | 0.18 / 18.4 |
| 2 faza, THP = 60, 1000 gün | ✅ | 65 | 127.0 → 180.5 | 10.9 | 2 | 0.19 / 21.7 |
| 3 faza, THP = 20, 600 gün | ✅ | 46 | 117.3 → 142.1 | 14.5 | 2 | 1.85 / 16.3 |

Oxunuşu:

* **Rəqs yoxdur** — BHP bütün qaçışda cəmi 2 dəfə istiqamət dəyişir.
* THP ilk bir neçə addımdan sonra hədəfə oturur (son addımlarda 19.98 / 59.98 bar).
* Maksimal sapma **ilk addımlardadır**: statik təxmin həqiqi axını bilmir və
  açıq birləşmə bir addım gecikir. Bu, açıq birləşmənin gözlənilən qiymətidir,
  gizlədilmir.
* BHP zamanla ARTIR — su payı böyüdükcə lülədəki sütun ağırlaşır, eyni THP
  üçün daha böyük BHP lazımdır. Fiziki cəhətdən düzgün istiqamətdir.
* 2 fazalı THP = 20 qaçışında ilk addımın THP-si `nan`-dır: statik təxmin
  (neft sütunu) həqiqi axın üçün kifayət etmədi. Növbəti addımda düzəlir.
* Monotonluq testlə yoxlandı: THP = 60 bar → THP = 10 bar-dan az neft.

### 5 · Öz təşəbbüsümlə verilmiş qərarlar

1. **THP bağlantıya ötürülmür** — yuxarıdakı səssiz-RATE təhlükəsinə görə.
2. **IMPES rədd edilir**, səssizcə sabit BHP = THP işlədilmir.
3. **Eclipse ixracı rədd edir** — WCONPROD THP rejimi VFPPROD cədvəli tələb edir,
   o isə modeldə yoxdur; yanlış deck yazmaq əvəzinə aydın xəta.
4. **THP ≥ lay təzyiqi — xəbərdarlıq**, xəta deyil (BHP rejimindəki qayda ilə eyni).

### 6 · Yoxlama

```
tests/test_thp_control.py        21 keçdi
tam dəst                         2447 keçdi, 1 buraxıldı, 1 xfailed (əvvəl 2426 + 21 yeni)
```

### Açıq qalan ⏳

* **Tam implicit THP birləşməsi** (THP qalıq tənliyinə daxil) — açıq birləşmənin
  ilk addım sapması onu tələb edərsə.
* **RATE quyusunda THP**, **Beggs-Brill sürüşməsi**, **VFPPROD idxalı** — dəyişmədi.
* **Vurucu quyuda THP** — V1-də yoxdur (diaqnostika xəta verir).
* Növbəti bloklar: **B5-b** (Pcog — sahibkar sonraya saxladı), **B7** (yekun
  doğrulama + SPE1).
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.


## 14 sentyabr 2026 — Seans 24: B5-b — üç fazalı kapilyar təzyiq (Pcow qoşuldu, Pcog əlavə olundu)

Sahibkarın tapşırığı: B4-B yoxlanıldı və təsdiqləndi; növbəti mərhələ B5-b
(kapilyar təzyiq funksiyaları). Üç fazalı yol sabitləşdiyi üçün ehtiyatla
davam edildi.

### 1 · Gözlənilməyən tapıntı — Pcow üç fazalı rejimdə SƏSSİZCƏ ATILIRDI

Plan yalnız Pcog əlavə etməyi nəzərdə tuturdu ("hazırda yalnız `pcow` var").
Kodu oxuyanda məlum oldu ki, üç fazalı yolda **Pcow da işləmir**:

* `ThreePhaseNewtonSolver.build_fluid` içində `pc = None` **SABİT** yazılmışdı;
* servis Pc provider-ini mühərrikə ötürürdü, mühərrik isə Nyutona vermirdi;
* `ThreePhaseFluidState.pc` sahəsi və qalıqdakı istifadəsi VAR idi — yəni
  boşluq yalnız bu bir sətirdə idi.

**Ölçmə ilə sübut olundu** (kod dəyişməzdən əvvəl): Pc = 1 bar ilə və Pc-siz
üç fazalı qaçış **bit-bit eyni** nəticə verirdi (RF 65.871437 %, 46 addım).
Yəni istifadəçi SCAL tabında Pc verirdi və qazlı modeldə heç nə dəyişmirdi.

Ona görə B5-b iki addıma bölündü: (1) Pcow-un qoşulması, (2) Pcog.

### 2 · Edilən

| Fayl | Nə |
|---|---|
| `simulation/implicit/three_phase_newton.py` | `capillary` / `gas_capillary` qəbul edir; `_capillary_terms()` — Pcow, ∂Pcow/∂Sw, Pcog, ∂Pcog/∂Sg |
| `simulation/implicit/three_phase_residual.py` | flüid vəziyyətinə `dpc_dsw`, `pcog`, `dpcog_dsg`; potensialda `Pg = Po + Pcog`; Jakobianda hər iki kapilyar hədd |
| `simulation/capillary.py` | **YENİ** `BrooksCoreyGasCapillaryProvider` — `pcog`, `dpcog_dsg` |
| `simulation/implicit/three_phase_engine.py` | Pcow Nyutona ötürülür; Pcog provider-i modeldən qurulur |
| `domain/reservoir_model.py` | `gas_capillary_parameters` + validasiya + "qaz yoxdursa Pcog işləməyəcək" xəbərdarlığı |
| `application/model_builder.py`, `serialization.py` | `gas_capillary` ötürücüsü və layihə faylında saxlanması |
| `ui/panels.py`, `ui/main_window.py` | SCAL tabının qaz bölməsində üç yeni sahə: Pcog Pe, λ, yuxarı hədd |
| `tests/test_three_phase_capillary.py` | **YENİ** — 21 test |

### 3 · Riyazi tərəf

```
Pw = Po − Pcow(Sw)          (əvvəldən iki fazalı yolda var idi)
Pg = Po + Pcog(Sg)          (YENİ)

Pcog(Sg) = Pe · S_L,n^(−1/λ),   S_L,n = (1 − Sg − Swc − Sorg) / (1 − Swc − Sorg)
```

Su-neft modeli ilə eyni Brooks-Corey forması, yalnız islanan faza mayedir.
Sg = 0-da Pcog = Pe — sabit əlavə, yəni qazsız hüceyrələr arasında SÜNİ AXIN
YARATMIR. Jakobianda kapilyar hədlər hər iki hüceyrəyə düşür (upstream-dən
asılı deyil):

```
∂F_su /∂Sw_a = −T·m_w,up·Pcow'_a      ∂F_su /∂Sw_b = +T·m_w,up·Pcow'_b
∂F_qaz/∂Sg_a = +T·m_g,up·Pcog'_a      ∂F_qaz/∂Sg_b = −T·m_g,up·Pcog'_b
```

Doymamış hüceyrədə 3-cü dəyişən Rs-dir və Pcog ondan asılı deyil — ona görə
orada törəmə **sıfırlanır**.

### 4 · Ölçmələr

**Jakobian (sonlu fərqə qarşı, maksimal nisbi xəta):**

| Vəziyyət | Kapilyarsız | Pcow | Pcog | Hər ikisi |
|---|---|---|---|---|
| tam doymuş, heterogen | — | < 10⁻⁵ | < 10⁻⁵ | < 10⁻⁵ |
| qarışıq (doymuş + doymamış) | 2.711×10⁻⁵ | 2.719×10⁻⁵ | 2.711×10⁻⁵ | 2.719×10⁻⁵ |

Qarışıq vəziyyətdəki 2.7×10⁻⁵ **B5-b-dən əvvəl də var idi** (dəyişən keçidi
olan hüceyrələrdə Rs sütununun mövcud dəqiqliyi). Kapilyarın əlavə etdiyi
xəta ~8×10⁻⁸-dir; test məhz bu FƏRQİ yoxlayır, mütləq həddi yox.

**Real 8×8 üç fazalı qaçış (600 gün):**

| Qaçış | Yığılma | Addım | RF, % | p_orta, bar |
|---|---|---|---|---|
| kapilyarsız | ✅ | 46 | 65.871437 | 239.91583 |
| Pcow = 1 bar | ✅ | 45 | 65.882462 | 239.92538 |
| Pcog = 1 bar | ✅ | 46 | 65.872458 | 239.91583 |
| Pcow = 3, Pcog = 6 bar | ✅ | 45 | 65.893475 | 239.94202 |

Oxunuşu: nəticə artıq kapilyar parametrlərdən ASILIDIR (əvvəl deyildi).
Bu bircins 5-nöqtəli modeldə təsir kiçikdir (RF-də 0.01–0.02 faiz bəndi) —
bu, gözlənilən nəticədir: kapilyar təzyiq özünü təbəqəli/heterogen modeldə
və keçid zonasında göstərir. Güclü kapilyarla da Nyuton yığılır, addım sayı
praktik olaraq dəyişmir.

### 5 · Öz təşəbbüsümlə verilmiş qərarlar

1. **Pcow düzəlişi Pcog-dan AYRI aparıldı** — biri mövcud funksiyanın
   qoşulması, digəri yeni fizikadır.
2. **Pcog Pcow ilə EYNİ dataklassdan** (`CapillaryParameters`) qidalanır —
   yeni parametr tipi yaradılmadı, UI və layihə faylı eyni üç sahəni işlədir.
3. **Son nöqtələr relperm provider-indən götürülür** (Sgc, Sorg, Swc) — kapilyar
   əyri ilə nisbi keçiricilik əyrisi eyni son nöqtələrə istinad etməlidir.
4. **Qaz olmayan modeldə Pcog verilməsi XƏBƏRDARLIQDIR** — səssizcə atılsaydı,
   istifadəçi yenə "parametr verdim, heç nə dəyişmədi" vəziyyətinə düşərdi.
   Bu, məhz indi düzəldilən səhvin təkrarı olardı.

### 6 · Yoxlama

```
tests/test_three_phase_capillary.py   21 keçdi
tam dəst                              2468 keçdi, 1 buraxıldı, 1 xfailed (əvvəl 2447 + 21 yeni)
```

### Açıq qalan ⏳

* **İlkin tarazlıqda Pcog yoxdur:** `equilibrium.py` GOC-da kəskin sərhəd
  qurur. Qaz papağı olan modeldə ilk addımlarda kiçik süni kapilyar axın olur.
  Keçid zonasının Pcog-dan qurulması ayrı işdir.
* **Eclipse ixracında SGOF/Pcog yoxdur** (SWOF-da Pcow var).
* **Cədvəldən Pcog** (`SaturationTableSet`) — hazırda yalnız Brooks-Corey
  analitik forması var.
* **Quyu həddində kapilyar** nəzərə alınmır (iki fazalı yolda da belədir).
* Növbəti blok: **B7** (yekun doğrulama + SPE1).
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.


## 15 sentyabr 2026 — Seans 25: THP idarəsinin SABİTLİYİ — nodal analiz

Sahibkar 41×41×3 heterogen, üç fazalı modeldə qeyri-stabillik bildirdi: qaz və
su sıçrayışından sonra təzyiq/debit şiddətli rəqs edir, Δt 0.03 günə düşür.
İlkin diaqnoz "açıq THP birləşməsi çökür" idi. **Ölçmə bunu qismən təkzib etdi**
və əsl səbəbi üzə çıxardı.

### 1 · Ölçmə ilə səbəbin ayrılması

Əvvəlcə sintetik 8×8×3 modeldə dörd variant müqayisə olundu (900 gün):

| Variant | Kəsilmə | Addım |
|---|---|---|
| THP + kapilyar | 27 | 131 |
| sabit BHP + kapilyar | 47 | 174 |
| THP, kapilyarsız | 42 | 135 |
| sabit BHP, kapilyarsız | **61** | 196 |

Yəni THP çıxarılanda vəziyyət PİSLƏŞİRDİ — sintetik modeldə THP günahkar deyildi.

Sonra sahibkarın öz `layihe.imx` faylı yükləndi. İlk tapıntı: **faylda qaz fazası
söndürülüb** (PVT-də Bg/qaz lözlüyü yoxdur), yəni saxlanmış model iki fazalıdır və
Pcog işləmir; sahibkar faylı fərqli testlər apararkən saxlamışdı. O fayl olduğu kimi
1000 gün qaçırıldı: **65 addım, kəsilmə yoxdur, 46 saniyə** — problem görünmür.

Eyni grid üzərində qaz fazası aktivləşdirildikdə simptom TAM təkrarlandı:
**300 addım, 125 uğursuz cəhd, Δt medianı 0.104 gün, minimum 0.0048** (sahibkarın
gördüyü 0.03 ilə eyni diapazon), 3575 saniyə.

### 2 · Kök səbəb — VLP əyrisi debitdən asılıdır, nəzarətçi isə bunu bilmirdi

Addım-addım jurnal göstərdi ki, PROD-1-in BHP-si iki dəyər arasında gedib-gəlir:
**183.8 bar** (durğun NEFT sütunu) və **~130 bar** (axan, qazlı yüngül sütun).

Ölçülmüş VLP əyrisi (PROD-1, 2015 m, THP = 20 bar):

| neft debiti, m³/gün | tələb olunan BHP, bar |
|---|---|
| 0 (durğun, qazsız) | 183.8 |
| 5 (qazla) | 71.4 |
| 125 | 79.6 |
| 500 | 144.6 |
| 1000 | 244.8 |

Nəzarətçi BHP-ni "keçən addımın debiti" ilə hesablayırdı, halbuki debitin özü
BHP-dən asılıdır. Nəticədə VLP əyrisi boyunca sıçrayış: quyu bağlanır → debit
sıfır → sütun ağırlaşır → BHP qalxır → quyu daha da bağlanır.

### 3 · Həll — üç addım (heç biri qalıq/Jakobiana toxunmur)

1. **NODAL ANALİZ.** İş nöqtəsi IPR ∩ VLP kəsişməsindən tapılır:
   `q(BHP) = J·(p_lay − BHP)` düz xətti ilə traverse birləşdirilir və
   `THP(BHP, q(BHP)) = THP_hədəf` tənliyi BHP üzrə ikiqat bölmə ilə həll olunur.
   BHP artdıqca debit azalır və sürtünmə düşür — alınan THP monotondur, ona görə
   kəsişmə yeganədir. Xərc əvvəlki ilə eynidir (~40 traverse).
2. **YARI-İMPLİCİT TƏKRAR.** Addım həll olunandan sonra BHP həmin addımın öz
   debitləri ilə yenilənir; fərq 1 bardan böyükdürsə addım eyni Δt ilə yenidən
   həll olunur (ən çox 2 təkrar). Təkrar yığılmasa əvvəlki həll saxlanılır —
   dövrə nəticəni heç vaxt pisləşdirə bilməz.
3. **AVTOMATİK BAĞLANMA.** Quyu HEÇ BİR debitdə THP hədəfinə çatmırsa, bağlantıların
   quyu indeksi sıfırlanır — hasilat dəqiq sıfır olur və quyu həddi Jakobiandan
   çıxır. Şərait düzələndə 2 barlıq ehtiyatla yenidən açılır (aç/bağla rəqsinə qarşı).

### 4 · Yol boyu tapılan iki SƏHV (hər ikisi ölçmə ilə)

* **Başlanğıcda səhv bağlanma.** İlk variantda PROD-1 t = 0-da bağlandı, çünki
  qərar DURĞUN (ən ağır) sütuna əsaslanırdı. Real quyu işə düşüb qaz verəndə
  sütun yüngülləşir. İndi bağlanma qərarı yalnız HƏQİQİ axan tərkib məlum
  olandan sonra verilir; başlanğıcda BHP lay təzyiqinin 5 bar altına qoyulur ki,
  quyu işə düşsün (`STARTUP_DRAWDOWN_BAR`).
* **Yalnız sıxma (clamp) kifayət etmir.** BHP-ni lay təzyiqində saxlamaq rəqsi
  dayandırsa da, sıfıra yaxın drawdown ədədi cəhətdən ən pis vəziyyətdir —
  ölçüldü: Δt 0.003 günə düşdü. Bağlanma məhz buna görə lazım oldu.

### 5 · Nəticə — sahibkarın modeli, 400 gün, qaz aktiv, THP = 20 bar

| | Köhnə | Yeni (nodal) |
|---|---|---|
| addım | 300 | **93** |
| uğursuz cəhd | 125 | **23** |
| Δt medianı | 0.104 gün | **3.958 gün** |
| Δt minimumu | 0.0048 gün | **0.125 gün** |
| vaxt | 3575 san | **421 san** |
| bağlanan quyu | — | yoxdur |

Hasilat: PROD-1 296.8 m³/gün (BHP 122.6), W-1 427.9 m³/gün (BHP 117.0),
kumulyativ neft 386 937 m³, RF 19.24 %, GOR 370.5, su payı 0 %, p_orta 180.3 bar.
Yarı-implicit dövrə 33 əlavə həll aparıb.

**Ən vacib cəhət:** köhnə qaçış "sürətli" görünən yerlərdə quyunu faktiki olaraq
bağlı saxlayırdı (fizika dayanmışdı). Yeni qaçışda hər iki quyu real hasilat verir
VƏ hesablama 8.5 dəfə sürətlidir.

### 6 · Yoxlama

```
tests/test_thp_control.py (7 yeni test)   28 keçdi
tam dəst                                  2475 keçdi, 1 buraxıldı, 1 xfailed
```

### Açıq qalan ⏳

* **Tam implicit THP birləşməsi LAZIM OLMADI** — ölçmə göstərdi ki, problem
  gecikmədə deyil, iş nöqtəsinin səhv tapılmasında idi. Qalıq/Jakobian toxunulmadı.
* IPR düz xətt kimi götürülür (J son addımdan). Vogel tipli əyri IPR ⏳ gələcəyə.
* Sürüşmə (Beggs-Brill), VFPPROD idxalı, RATE quyusunda THP — dəyişmədi.
* Növbəti blok: **B7** (yekun doğrulama + SPE1).
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.


## 15 sentyabr 2026 — Seans 26: B7 addım 1 — QAZ VURULMASI

B7 (SPE1 etalonu) üçün mühərrikin imkanları tələblərlə tutuşduruldu və
**üç boşluq** tapıldı:

| SPE1 tələbi | Vəziyyət |
|---|---|
| künc quyusuna 100 MMscf/gün **qaz vurulması** | ❌ vurucular yalnız su vururdu → **BU SEANSDA BAĞLANDI** |
| istismarçıda **debit + BHP limiti** (20 000 STB/gün, min 1000 psia) | ❌ rejim keçidi yoxdur → növbəti addım |
| **səth** neft debiti hədəfi | ⚠️ bizdə RATE = maye, lay həcmi |

Sahibkarın seçimi: ardıcıl getmək — əvvəl qaz vurulması, ölçmək, commit etmək;
sonra BHP limiti; sonra SPE1. Səbəb: iki yeni xüsusiyyət eyni anda üç fazalı
qalıq/Jakobiana girsəydi, səhv çıxanda mənbəyi ayırmaq çətinləşərdi.

### 1 · Tapıntı — `injected_phase` VAR İDİ, lakin OXUNMURDU

`WellControl.injected_phase` domendə mövcud idi və layihə faylında saxlanılırdı,
lakin heç bir simulyasiya kodu ona baxmırdı: `ThreePhaseWellModel.well_rates`
vurucunu həmişə su kimi işlədirdi (`endpoint_water_mobility / mu_w` → `water[cell]`).
Yəni istifadəçi qaz seçsə belə, quyu SU vururdu və heç bir xəbərdarlıq yox idi.

### 2 · Edilən

| Fayl | Nə |
|---|---|
| `simulation/well_model.py` | `WellConnection.injected_phase` — faza bağlantıya çıxarıldı |
| `simulation/implicit/three_phase_residual.py` | qaz vurucusu **qaz tənliyinə** yazır; mobillik = `krg_end / μ_g`, həcm `Bg`; Jakobianda vurulan fazaya görə SƏTİR seçilir (su → 0, qaz → 2) |
| `simulation/implicit/three_phase_newton.py` | qazın son nöqtəsi relperm provider-indən (`gas.krg_end`) |
| `domain/reservoir_model.py` | qaz fazası olmayan modeldə qaz vurulması — **bloklayıcı XƏTA** |
| `ui/panels.py` | quyular cədvəlində «Vurulan faza» sütunu (SU / QAZ), layihə faylında saxlanılır |
| `tests/test_gas_injection.py` | **YENİ** — 15 test |

### 3 · Ölçmələr

**Jakobian (sonlu fərqə qarşı, maksimal nisbi xəta):**

| Yol | Xəta |
|---|---|
| qaz vurucusu, BHP rejimi | < 10⁻⁵ |
| qaz vurucusu, RATE rejimi | < 10⁻⁵ |
| su vurucusu (TOXUNULMAYIB) | 1.343141869382806×10⁻⁵ |

Su yolunun xətası **B7-dən əvvəlki kodda ÖLÇÜLDÜ** (HEAD-dən ayrıca iş nüsxəsi
qurularaq) və rəqəm **rəqəminə eynidir** — yəni su vurulması ədədi cəhətdən
zərrə qədər dəyişməyib. Həmin 1.34×10⁻⁵ vurucu mobilliyinin doyumluluq
törəməsinin nəzərə alınmamasından gəlir (köhnə sadələşdirmə, açıq sənədlənib).

**Uc-uca:** qaz vuran 5×5 model yığılır, vurucunun ətrafında qaz doymuşluğu
0.1-dən yuxarı qalxır və lay təzyiqi vurulmayan haldan YÜKSƏK qalır.

### 4 · Öz təşəbbüsümlə verilmiş qərar

**Qaz fazası olmayan modeldə qaz vurulması XƏTADIR** (xəbərdarlıq deyil).
Səbəb: iki fazalı mühərrik vurulan fazanı ümumiyyətlə oxumur, yəni quyu
səssizcə su vurardı. Bu, Pcog-dakı (Q-16) eyni prinsipin davamıdır, lakin
orada nəticə "parametr işləmir", burada isə "TAMAMİLƏ BAŞQA flüid vurulur" —
ona görə xəbərdarlıq kifayət etmir.

### 5 · Yoxlama

```
tests/test_gas_injection.py   15 keçdi
tam dəst                      2490 keçdi, 1 buraxıldı, 1 xfailed
```

### Açıq qalan ⏳

* **B7 addım 2:** RATE rejimində **BHP limiti** və rejim keçidi (SPE1-in istismarçısı).
* **B7 addım 3:** səth debiti hədəfi + SPE1 modelinin qurulması və nəşr olunmuş
  nəticələrlə tutuşdurulması.
* Qaz vurucusunda mobillik SON NÖQTƏ yaxınlaşmasıdır (su vurulmasında olduğu kimi).
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.

## 15 sentyabr 2026 — Seans 27: RATE hədəfinin perforasiyalara bölünməsi (B7 addım 2-yə hazırlıq)

B7 addım 2-yə (RATE rejimində BHP limiti) başlayarkən RATE budaqları oxundu
və **səssiz səhv** tapıldı. BHP limiti quyunun debitindən hesablanacağı üçün
əvvəlcə bu səhv bağlandı.

### 1 · Tapıntı — RATE hədəfi HƏR perforasiyaya TAM yazılırdı

Qalıqda, Jakobianda və IMPES-də RATE budağı bağlantı üzrə dövrdə
`total = -abs(connection.target)` yazırdı — yəni hər perforasiya quyunun
BÜTÜN hədəfini alırdı. Ölçüldü (5×5 model, istismarçı RATE = 50 m³/gün,
lay həcmi):

| Perforasiya sayı | Faktiki maye debiti | Nisbət |
|---|---|---|
| 1 | 50.0000 | 1.0000 |
| 3 | 150.0000 | **3.0000** |

UI-də «Perf üst/alt» boş qalanda quyu bütün təbəqələrdə perforasiya olunur,
yəni çox təbəqəli modeldə UI-dən yaradılan hər RATE quyusu bu səhvdən
təsirlənirdi. SPE1-ə təsiri yox idi (orada quyular tək perforasiyalıdır).

Sahibkara üç variant təqdim olundu (düzəlt / yalnız sənədləşdir / açıq xəta
ver). **Sahibkarın seçimi: əvvəl düzəlt**, sonra BHP limiti.

### 2 · Edilən

| Fayl | Nə |
|---|---|
| `simulation/well_constraints.py` | **YENİ** — `assign_rate_shares`, `needs_rate_allocation` |
| `simulation/well_model.py` | `WellConnection.rate_share` (defolt 1.0); ilkin pay WI ilə |
| `simulation/implicit/residual.py`, `jacobian.py` | RATE hədəfi paya vurulur; `connection_mobilities` |
| `simulation/implicit/three_phase_residual.py` | eyni — `well_rates` və quyu Jakobianı; `connection_mobilities` |
| `simulation/implicit/engine.py`, `three_phase_engine.py` | pay hər addımdan ƏVVƏL yığılmış vəziyyətin λ-sı ilə yenilənir |
| `simulation/impes_engine.py` | pay təzyiq addımında, BHP budağındakı eyni mobilliklə |
| `tests/test_rate_allocation.py` | **YENİ** — 13 test |

Qayda: `pay_c = WI_c·λ_c / Σ WI·λ`, λ addımın əvvəlindən, Nyuton daxilində
sabit (bax Q-19). Tək perforasiyada pay dəqiq 1.0-dır.

### 3 · Ölçmələr

**Düzəlişdən sonra** eyni ölçmə: 1 perforasiya → 50.0000, 3 perforasiya →
50.0000 (nisbət 1.0000).

**Jakobian (sonlu fərqə qarşı, maksimal nisbi xəta):**

| Yol | Xəta |
|---|---|
| iki fazalı, iki perforasiyalı RATE istismarçısı + vurucusu, pay 0.25/0.75 və 2/7 / 5/7 | < 10⁻⁸ (dəqiq) |
| üç fazalı RATE istismarçısı, **tək** perforasiya (düzəlişə dəxli yoxdur) | 0.4893141064217388 |
| üç fazalı RATE istismarçısı, iki perforasiya | 0.3239019497740145 |

⚠️ Üç fazalı RATE istismarçısının Jakobianı **əvvəldən** kobud təqribidir:
tək perforasiyada da xəta 0.49-dur (RATE budağında sərbəst qazın təzyiq/Sg
törəmələri sıfır qoyulub). Pay bunu pisləşdirmir; düzəlişi ayrıca işdir ⏳.

**Uc-uca:** 3 və 1 perforasiyalı quyunun 10 gündə cəmi maye hasilatı üç
mühərrikdə (tam implicit, IMPES, üç fazalı) ±5 % daxilində üst-üstə düşür —
testlə kilidlənib.

### 4 · Yoxlama

```
tests/test_rate_allocation.py   13 keçdi
tam dəst                        2503 keçdi, 1 buraxıldı, 1 xfailed
```

### 5 · Nəticəsi dəyişən modellər

Tək perforasiyalı quyular bit-bit eynidir (tam dəst dəyişməz keçdi). **Çox
perforasiyalı RATE quyusu olan modellərin nəticəsi DƏYİŞİR** — düzgün
tərəfə: quyu artıq hədəfin özünü hasil edir/vurur.

### Açıq qalan ⏳

* **B7 addım 2:** RATE rejimində BHP limiti — bu seansın davamı.
* Üç fazalı RATE istismarçısının Jakobian xətası (0.49) — ayrıca iş.
* Təbəqələr arası paylanma bir addım gecikir (λ əvvəlki addımdandır); quyunun
  CƏMİ debiti dəqiqdir.
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.

## 15 sentyabr 2026 — Seans 27 (davamı): Eclipse ixracında quyu sətirləri

Təhvil sənədi B7 addım 2 üçün `WCONPROD` BHP limit sahəsinin yoxlanılmasını
tələb edirdi. Yoxlanıldı və **iki səssiz səhv** tapıldı. BHP limitindən ayrı
commit kimi düzəldildi.

### 1 · Tapıntılar

| # | Nə yazılırdı | Niyə səhvdir |
|---|---|---|
| 1 | istismarçı: `'LRAT' 2* debit 2* 1.0` | WCONPROD sütunları 4 ORAT · 5 WRAT · 6 GRAT · 7 LRAT · 8 RESV · 9 BHP — `2*` iki sütunu buraxır, debit **6-cı (QAZ debiti)** sütununa düşürdü |
| 1b | rejim `LRAT` / `RATE` | bizim RATE hədəfi **lay həcmidir** (UI: «RATE → m³/gün (rezervuar həcmi)»), LRAT və WCONINJE `RATE` isə səth debitidir |
| 2 | qaz vurucusu: `'WATER' ...` | deck iki fazalıdır (`OIL`/`WATER`, qaz PVT-si yoxdur) — qaz vuran quyu səssizcə SU vurucusu kimi yazılırdı (B7 addım 1-dən sonra yaranan boşluq) |

Repoda `WCONPROD`/`WCONINJE`-ni oxuyan kod yoxdur (idxal bu açar sözləri
tanımır), ona görə düzəliş heç bir dövrə testinə təsir etmir.

### 2 · Edilən

| Fayl | Nə |
|---|---|
| `io/eclipse_export.py` | istismarçı `'RESV' 4* debit 1.0` (debit 8-ci, BHP 9-cu sütunda); vurucu `'RESV' 1* debit 1000.0` (debit 6-cı, BHP 7-ci sütunda); qaz vurucusu — `ValueError` |
| `tests/test_eclipse_well_controls.py` | **YENİ** — 4 test, sətir `n*` defoltları açılaraq sütun-sütun yoxlanır |
| `ECLIPSE_IO.md` | quyu rejimlərinin xəritəsi |

### 3 · Yoxlama

```
tests/test_eclipse_well_controls.py                        4 keçdi
test_eclipse_io.py + test_scal_tables.py + test_thp_control.py   81 keçdi
```

Tam dəst bu commit üçün ayrıca işlədilmədi — dəyişiklik yalnız ixrac
mətnidir, simulyasiya koduna toxunmur. Tam dəst BHP limiti commit-indən
sonra işlədiləcək.

### Açıq qalan ⏳

* Deck-in üç fazalı ixracı (`GAS`, `PVDG`/`PVTO`, `SGOF`) — mövcud backlog.
* BHP limiti ixrac sahəsi — B7 addım 2 commit-ində.

## 15 sentyabr 2026 — Seans 27 (davamı): B7 addım 2 — RATE rejimində BHP LİMİTİ

Təhvil sənədindəki tövsiyə olunan dizayn qəbul edildi: rejim keçidi
**addımlar arasında**, qalıq və Jakobiana toxunulmadan (Q-15-dəki THP
yanaşmasının təkrarı). Bu işə başlamazdan əvvəl iki səssiz səhv bağlandı
(RATE payı — Q-19, Eclipse sütunları — Q-20), çünki limit məhz RATE debitindən
hesablanır və ixrac məhz limit sahəsini yazır.

### 1 · Edilən

| Fayl | Nə |
|---|---|
| `domain/wells.py` | `WellControl.bhp_limit: Optional[float] = None`; mütləq təzyiq kimi doğrulanır |
| `application/serialization.py` | `.imx`-də `bhp_limit`; açarı olmayan köhnə fayl → `None` |
| `simulation/well_model.py` | `WellConnection.bhp_limit` — yalnız RATE rejimli quyuda |
| `simulation/well_constraints.py` | **`BhpLimitController`**, `RATE_RESTORE_MARGIN_BAR = 2.0` |
| `simulation/implicit/engine.py`, `three_phase_engine.py` | `_bhp_limit_loop` (THP dövrəsindən sonra), `bhp_limit_resolves`, `_connection_mobilities` |
| `simulation/results.py` | `well_control_mode` (hər addımda "RATE"/"BHP"); limitli quyuda `well_bhp` doldurulur |
| `simulation/impes_engine.py`, `application/simulation_service.py` | IMPES + BHP limiti — **açıq imtina** (texniki və istifadəçi dilində) |
| `domain/reservoir_model.py` | limitli RATE quyusu "təzyiq idarəsi" sayılır; limit BHP/THP rejimində və ya lay təzyiqinin səhv tərəfində — xəbərdarlıq |
| `io/eclipse_export.py` | limit `WCONPROD` 9-cu / `WCONINJE` 7-ci sütuna yazılır (verilməyibsə əvvəlki 1.0 / 1000.0) |
| `ui/panels.py` | quyular cədvəlində «BHP limiti» sütunu (boş = limit yoxdur) |
| `tests/test_bhp_limit.py` | **YENİ** — 21 test |

### 2 · Qaydalar (bax Q-21)

Yığılmış addımdan sonra hədəf debiti verən BHP Peaceman-ın tərsindən tapılır:

    istismarçı:  BHP = (Σ WI·λ·p − q) / Σ WI·λ
    vurucu:      BHP = (Σ WI·λ·p + q) / Σ WI·λ

* RATE quyusunda bu BHP limiti pozursa (və ya quyu axa bilmirsə) → `mode = BHP`,
  `target = limit`, addım eyni Δt ilə yenidən həll olunur;
* limitdə olan quyu RATE-ə yalnız hədəf **2 bar ehtiyatla** əldə olunanda qayıdır;
* hər quyu bir addımda ən çox bir dəfə rejim dəyişir.

### 3 · Ölçmələr

**Kəsilməzlik:** tələb olunan BHP-də BHP rejimi eyni cəmi debiti verir —
nisbi fərq < 10⁻⁹ (testlə kilidlənib).

**Tükənən lay** (8×8×2, vurucu yoxdur, istismarçı RATE = 150 m³/gün lay həcmi,
ilkin təzyiq 250 bar, 600 gün, tam implicit):

| Qaçış | Addım | Δt medianı | Δt min | Δt kəsilməsi | Əlavə həll | Keçid |
|---|---|---|---|---|---|---|
| iki fazalı, limit 1 bar (praktik limitsiz) | 36 | 20.000 | 1.0000 | 0 | 1 | 1 (t ≈ 69 gün) |
| iki fazalı, limit 180 bar | 36 | 20.000 | 1.0000 | 0 | 1 | 1 (addım #5, t = 20.8 gün) |
| üç fazalı, limit 1 bar | 228 | 0.634 | 0.1691 | 17 | 1 | 1 (t = 155.7 gün) |
| üç fazalı, limit 180 bar | 36 | 20.000 | 1.0000 | 0 | 1 | 1 (addım #5, t = 20.8 gün) |

* Limit 180 bar olan qaçışlarda keçiddən sonra BHP min = maks = **180.000000**;
  keçiddən əvvəl tələb olunan BHP-nin minimumu 183.903 bar — yəni qeydə alınan
  BHP heç vaxt limitdən aşağı deyil.
* Keçid addımında səth maye debiti 120.25 → 43.64 m³/gün (addım BHP = 180 ilə
  yenidən həll olundu), sonra lay təzyiqi limitə yaxınlaşdıqca sıfıra enir.
* Rejim seriyada **bir dəfə** dəyişir — rəqs yoxdur; keçid Δt-ni kəsmədi.
* Üç fazalı "limit 1 bar" qaçışının 228 addımı limitdən DEYİL: təzyiq doyma
  təzyiqindən (140 bar) aşağı düşür və qaz ayrılır. Limit 180 bar təzyiqi Pb-dən
  yuxarı saxladığı üçün həmin qaçış 36 addımdır — bu iki sətir mexanizmin
  xərcini deyil, fizikanın fərqini göstərir.

**Toxunulmayan limit:** vurucu (BHP 255 bar) təzyiqi 239 barda saxlayanda limit
180 bar heç vaxt işə düşmür və nəticə limitsiz qaçışla **eynidir** (iki fazalı:
cəmi maye 75936, orta təzyiq 238.84; üç fazalı: 76163, 238.85 — hər iki
halda ölçüldü, testdə massivlər bərabərliklə yoxlanılır).

### 4 · Yoxlama

```
tests/test_bhp_limit.py   21 keçdi
tam dəst                  2528 keçdi, 1 buraxıldı, 1 xfailed
```

### 5 · Öz təşəbbüsümlə verilmiş qərarlar

* **Limitli RATE quyusu "ən azı bir təzyiq idarəli quyu" qaydasını ödəyir.**
  Əks halda SPE1 (iki RATE quyusu, hər ikisində BHP həddi) bloklanardı.
* **BHP/THP rejimində verilmiş limit — xəbərdarlıq, xəta deyil.** Hesab yanlış
  deyil, sadəcə parametr işləmir (Q-16, Qərar 4 ilə eyni prinsip). IMPES-də isə
  limit səssizcə atılardı və quyu BHP-ni pozardı — orada **xəta**.

### Açıq qalan ⏳

* **B7 addım 3:** səth neft debiti hədəfi + SPE1 modeli (parametrlər mənbədən
  yoxlanılmalıdır — təhvil sənədi §6.2).
* Vurucu limiti yalnız vahid testlərlə yoxlanılıb; uc-uca vurucu ssenarisi yoxdur.
* Tərs düstur BHP rejimindəki `min(q, 0)` kəsməsini nəzərə almır (çox təbəqəli
  quyuda təbəqələr arası təzyiq fərqi böyük olanda təqribidir).
* THP quyusu ilə BHP limitli quyu eyni modeldə olanda: limit dövrəsi addımı
  yenidən həll edirsə, THP quyusunun BHP-si həmin yeni həllə görə yenilənmir
  (THP dövrəsi limitdən ƏVVƏL işləyir) — bir addımlıq gecikmə, ölçülməyib.
* `well_control_mode` CSV/JSON ixracına və dashboard-a çıxarılmayıb.
* RATE quyusunda THP hələ hesablanmır (hidravlika yalnız BHP/THP rejimini tanıyır).
* Üç fazalı RATE istismarçısının Jakobian xətası (0.49) — Seans 27-dən.
* **`b2a795e`** — sahibkar digər maşında özü xilas edəcək.
