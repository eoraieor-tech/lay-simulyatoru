# 3D görüntü

## Niyə matplotlib

`pyvista`/`vtk` daha güclüdür, lakin əlavə quraşdırma tələb edir və
bu mühitdə mövcud deyil. `matplotlib` onsuz da layihənin asılılığıdır
və Qt ilə birlikdə interaktiv fırlatma verir.

Əvəzində performans üçün xüsusi iş görülüb (aşağıda).

## Əsas texniki qərar — yalnız görünən üzlər

Naiv yanaşmada hər hüceyrə üçün 6 üz çəkilir. 54 000 hüceyrəli
modeldə bu, **324 000 çoxbucaqlıdır** və matplotlib onu dəqiqələrlə
çəkir.

Həll: üz yalnız o halda çəkilir ki,
- modelin sərhədindədirsə, **yaxud**
- qonşusu filtr ilə gizlədilibsə.

### Ölçülmüş nəticə

| Model | Hüceyrə | Bütün üz | Çəkilən | Azalma | Vaxt |
|---|---|---|---|---|---|
| 20×20×5 | 2 000 | 12 000 | 1 200 | 10× | 0.02 san |
| 40×40×10 | 16 000 | 96 000 | 4 800 | 20× | 0.05 san |
| 60×60×15 | 54 000 | 324 000 | 10 800 | **30×** | 0.20 san |

Azalma ölçü ilə **yaxşılaşır**, çünki səth həcmə nisbətdə azalır
(kub qanunu). Bu, üsulun böyük modellərdə daha da faydalı olduğunu
göstərir.

Doğruluq testlə yoxlanılıb: filtrsiz halda çəkilən üz sayı analitik
düsturla üst-üstə düşür — `2(nx·ny + ny·nz + nx·nz)`.

## İşıqlandırma

3D görüntünün ən çox itirdiyi məlumat **formadır**. Bütün üzlər eyni
parlaqlıqda olanda model yastı bir ləkə kimi görünür və hansı hissənin
öndə, hansının arxada olduğu bilinmir.

Həll — Lambert işıqlandırması: üzün parlaqlığı onun işıq mənbəyinə
baxma bucağından asılıdır.

    parlaqlıq = 1 − güc + güc · max(0, n⃗ · l⃗)

`n⃗` üzün normalı, `l⃗` işıq istiqaməti (yuxarı-sol-öndən). Altı üz
istiqaməti → altı fərqli parlaqlıq (test ilə yoxlanılıb).

| İşıq gücü | Nəticə |
|---|---|
| 0 | düz rəng — **rəqəm oxumaq üçün ən dəqiq** |
| 0.45 (defolt) | balanslı |
| 1.0 | güclü kölgə — **forma ən yaxşı görünür** |

Vacib: işıqlandırma rəng xəritəsinin ÜZƏRİNƏ tətbiq olunur, onu əvəz
etmir. Yəni dəyər-rəng uyğunluğu qorunur, sadəcə parlaqlıq dəyişir.

## İdarəetmə elementləri

| Element | Nə edir |
|---|---|
| **Xassə** | Sw, təzyiq, keçiricilik, məsaməlilik, dərinlik |
| **Zaman sürgüsü** | Simulyasiya anları arasında hərəkət |
| **Kəsim həddi** | Faiz — aşağı dəyərli hüceyrələri gizlədir |
| **K … K** | Təbəqə aralığı (şaquli kəsim) |
| **Z×** | Şaquli mübaliğə (1–30) |
| **Quyular** | Lülə və perforasiyaların göstərilməsi |
| **Kənarlar** | Hüceyrə sərhədlərinin xətti |
| **Baxış** | Hazır bucaqlar: İzometrik, Üstdən, Yandan (X/Y), Künc |
| **İşıq** | Kölgənin gücü (0–100 %) |
| **Şəffaflıq** | 20–100 %; daxili qatları göstərir |
| **Görünüşü sıfırla** | İzometrik + defolt parametrlər |

"Baxış" siyahısında **Sərbəst** seçilibsə, siçanla qurulan bucaq
qorunur — xassə və ya zaman dəyişəndə kamera sıfırlanmır.

Siçanla sürüşdürərək model fırladılır (matplotlib-in öz mexanizmi).

### Kəsim həddi niyə vacibdir

Filtr sadəcə gizlətmir — o, **daxili strukturu açır**. Hüceyrə
gizlədiləndə qonşusunun üzü "görünən" olur və modelin içi görünür.

Praktik istifadə: `Sw` xassəsində həddi 60-70 %-ə qaldırsan, yalnız
su ilə süpürülmüş zona qalır və cəbhənin 3D forması aydın görünür.

### Şaquli mübaliğə

Real rezervuarlar yastıdır: 800 × 800 m sahədə 30 m qalınlıq.
Həqiqi nisbətdə model kağız vərəqi kimi görünür. `Z×` onu şaquli
istiqamətdə uzadır — sənaye vizualizatorlarında standart vasitədir.

Diqqət: mübaliğə yalnız **görüntüyə** təsir edir, hesablamaya yox.
Test bunu qoruyur: X və Y koordinatları toxunulmaz qalır.

## Şaquli mübaliğə nəyə tətbiq olunur

Real rezervuarlar yastıdır: 2000 × 2000 m sahədə 100 m qalınlıq.
Həqiqi nisbətdə model kağız vərəqi kimi görünür.

**Mübaliğə yalnız görüntü nisbətinə (`box_aspect`) tətbiq olunur,
koordinatlara YOX.** İlk versiyada əksinə idi: hüceyrələr
`dz × mübaliğə` hündürlüyündə çəkilirdi, ox etiketləri isə həqiqi
dərinliyi göstərirdi — ikisi bir-birinə zidd idi.

Nümunə: tavan 1599 m, bir təbəqə, DZ = 100 m. Lay 1599–1699 m
arasındadır. `Z×3` ilə ox 1349–1949 göstərirdi — üç dəfə yanlış
interval.

İndi ox həmişə **1594–1704 m** göstərir (lay + 5 % marja) və `Z×`
ona təsir etmir; yalnız qutunun nisbətini dəyişir (0.067 → 0.635).

Test bunu qoruyur: mübaliğə 1-dən 8-ə dəyişəndə çəkilən
çoxbucaqlıların koordinatları **eyni qalır**.

## Dərinlik oxu

İşarələr **təbəqə sərhədlərində** qoyulur, ixtiyari addımlarda yox.
5 təbəqə × 80 m, tavan 2000 m üçün: **2000, 2080, 2160, 2240, 2320,
2400**. Beləliklə hansı hüceyrənin harada bitdiyi birbaşa görünür.

12-dən çox sərhəd olanda işarələr seyrəldilir.

## Siçanla idarəetmə

CMG Builder ilə eyni təyinat:

| Düymə | Nə edir |
|---|---|
| **Sol** basıb sürüşdür | Modeli fırladır |
| **Orta** basıb sürüşdür | Səhnəni sürüşdürür (Pan) |
| **Sağ** basıb sürüşdür | Yaxınlaşdırır |
| **Çarx** | Yaxınlaşdırır / uzaqlaşdırır |

### Alət paneli niyə çıxarıldı

3D tabında matplotlib naviqasiya paneli **zərərlidir**. `Axes3D._on_move`
başlanğıcda yoxlayır:

    get_navigate_mode() is not None  ->  return

Yəni panelin "pan" və ya "zoom" rejimi aktiv olan kimi siçanla
fırlatmaq **tamamilə dayanır**. İstifadəçi düyməni basır və modelin
niyə fırlanmadığını anlamır.

Ona görə 3D tabında yalnız "Şəkli saxla…" düyməsi qalıb. 2D xəritədə
panel faydalıdır və orada saxlanılıb. Test bunu qoruyur.

### Çarxın sürəti

Çarx hər addımda səhnəni yenidən çəkmir — yalnız kamera nisbətini
dəyişir (`apply_zoom`). 50 000 hüceyrəli modeldə tam yenidən çəkiliş
saniyələr çəkərdi və hərəkət kəsikli olardı.

## Yaxınlaşdırma

`Yaxınlaşdırma` sürgüsü ekran böyütməsidir:

| Dəyər | Nəticə |
|---|---|
| **100 %** (defolt) | tam model çərçivəyə sığır |
| 150–250 % | modelə yaxından baxış, detallar aydınlaşır |
| 25–75 % | uzaqdan ümumi görünüş |

Yaxınlaşdırma **yalnız görüntüyə** təsir edir: ox hədləri, dərinlik
işarələri və çəkilən həndəsə dəyişmir. Test bunu qoruyur.

matplotlib-in 3D oxu defolt olaraq geniş boş kənar buraxır, ona görə
100 %-ə baza əmsalı (`BASE_FIT = 1.35`) tətbiq olunur — əks halda model
çərçivənin yalnız yarısını tutardı.

Şaquli miqyas idarəsi (`Z×`) çıxarıldı: model həmişə **həqiqi nisbətdə**
göstərilir. Hündür qutuda plan sıxılırdı və nəzarət qarışıqlıq yaradırdı.

## Məhdudiyyətlər

- **Şəffaflıq yoxdur.** matplotlib-in 3D-də şəffaflıq sıralaması
  düzgün işləmir, ona görə hüceyrələr tam qeyri-şəffafdır. Daxili
  strukturu görmək üçün kəsim həddi və təbəqə aralığı işlədilir.
- **Çox böyük modellər.** 100 000+ hüceyrədə çəkmə bir neçə saniyə
  çəkir. Bu halda təbəqə aralığını daraltmaq lazımdır.
- **Fırlatma sürəti.** Hər fırlatmada səhnə yenidən çəkilir; böyük
  modellərdə bu hiss olunur.

`pyvista` quraşdırılsa, gələcəkdə ikinci renderer kimi əlavə edilə
bilər — `VolumeRenderer` interfeysi buna hazırdır.## Siçanla idarəetmə

CMG Builder ilə eyni təyinat:

| Düymə | Nə edir |
|---|---|
| **Sol** basıb sürüşdür | Modeli fırladır |
| **Orta** basıb sürüşdür | Səhnəni sürüşdürür (Pan) |
| **Sağ** basıb sürüşdür | Yaxınlaşdırır |
| **Çarx** | Yaxınlaşdırır / uzaqlaşdırır |

### Alət paneli niyə çıxarıldı

3D tabında matplotlib naviqasiya paneli **zərərlidir**. `Axes3D._on_move`
başlanğıcda yoxlayır:

    get_navigate_mode() is not None  ->  return

Yəni panelin "pan" və ya "zoom" rejimi aktiv olan kimi siçanla
fırlatmaq **tamamilə dayanır**. İstifadəçi düyməni basır və modelin
niyə fırlanmadığını anlamır.

Ona görə 3D tabında yalnız "Şəkli saxla…" düyməsi qalıb. 2D xəritədə
panel faydalıdır və orada saxlanılıb. Test bunu qoruyur.

### Çarxın sürəti

Çarx hər addımda səhnəni yenidən çəkmir — yalnız kamera nisbətini
dəyişir (`apply_zoom`). 50 000 hüceyrəli modeldə tam yenidən çəkiliş
saniyələr çəkərdi və hərəkət kəsikli olardı.

## Yaxınlaşdırma

`Yaxınlaşdırma` (25–500 %, defolt **100 %**) kameranı modelə
yaxınlaşdırır. 100 % — model tam kadra sığır; artırdıqca detallara
yaxından baxılır.

Yalnız kameraya təsir edir: koordinatlar, ox hədləri və çəkilən üzlər
dəyişmir. Test bunu qoruyur.

Şaquli mübaliğə (`Z×`) idarəetməsi **çıxarıldı** — model həmişə həqiqi
nisbətdə göstərilir. Yaxından baxmaq lazım olanda yaxınlaşdırma
işlədilir.

## Məhdudiyyətlər

- **Şəffaflıq yoxdur.** matplotlib-in 3D-də şəffaflıq sıralaması
  düzgün işləmir, ona görə hüceyrələr tam qeyri-şəffafdır. Daxili
  strukturu görmək üçün kəsim həddi və təbəqə aralığı işlədilir.
- **Çox böyük modellər.** 100 000+ hüceyrədə çəkmə bir neçə saniyə
  çəkir. Bu halda təbəqə aralığını daraltmaq lazımdır.
- **Fırlatma sürəti.** Hər fırlatmada səhnə yenidən çəkilir; böyük
  modellərdə bu hiss olunur.

`pyvista` quraşdırılsa, gələcəkdə ikinci renderer kimi əlavə edilə
bilər — `VolumeRenderer` interfeysi buna hazırdır.

---

# VTK motoru (ResInsight tipli) — v53

İstifadəçi ResInsight-ın 3D görünüşünü bəyəndi və eynisini istədi.
ResInsight VTK üzərində qurulub (C++); eyni kitabxananın Python
versiyası mövcuddur.

## Niyə matplotlib kifayət etmir

matplotlib 3D üçün nəzərdə tutulmayıb — 2D çəkici olub 3D-ni "təqlid"
edir:

- **Dərinlik sıralaması dəqiq deyil** (z-buffer yoxdur, mərkəz
  nöqtəsinə görə təxmini sıralama). Bu problemin nəticəsi artıq
  sənədləşib: fault görüntüsündə dolğu əvəzinə kontur işlətməyə məcbur
  olduq, çünki dolğu "cırıq" görünürdü (bax `FAULTS.md`).
- **Hər fırlatmada bütün üzlər yenidən sıralanır və çəkilir** — böyük
  gridlərdə dözülməz ləng.
- **İşıqlandırma primitivdir** — kölgə əl ilə hesablanır.

VTK isə əsl 3D qrafika kitabxanasıdır (OpenGL üzərində) — ResInsight,
ParaView, Petrel onu işlədir.

## Dizayn: ƏVƏZ ETMİR, YANINDA YAŞAYIR

`vtk_volume.py` mövcud `volume.py`-ı əvəz ETMİR. İstifadəçi "3D
görüntü" tabında motoru seçir:

    VTK (sürətli)  ·  matplotlib

VTK **məcburi asılılıq deyil**: quraşdırılmayıbsa (`available()` False)
siyahıda görünmür və proqram tam olduğu kimi işləyir. Bu ayrılıq
sayəsində `pip install vtk` etməyən istifadəçi heç nə itirmir.

Üstəlik iki qat müdafiə var: (1) VTK-nın Qt körpüsü qurula bilməsə
(OpenGL sürücüsü, uzaq masaüstü), motor siyahıdan avtomatik silinir;
(2) çəkiliş zamanı xəta olsa, proqram matplotlib-ə qayıdır və səbəbi
loqa yazır — istifadəçi boş ekranla qalmır.

## Performans fərqinin mənbəyi

matplotlib motorunda **hər yeniləmədə** bütün həndəsə yenidən
hesablanırdı. VTK səhnəsi BİR DƏFƏ qurulur (`VtkReservoirScene`),
sonra yalnız dəyərlər yenilənir (`update_values()`). Zaman slider-ini
sürüşdürəndə fərq dərhal hiss olunur.

Kəsim həddi/K-təbəqə filtri də hüceyrələri SİLMİR — VTK-nın
`BlankCell` mexanizmi ilə "gizlədir", həndəsə toxunulmaz qalır.

## Rəng xəritəsi

`_build_lookup_table()` HƏM ad (str), HƏM matplotlib `Colormap`
obyekti qəbul edir. İnterfeys `MapRenderer._select_volume()`-dən
obyekt alır — bu halda matplotlib-in ÖZ rəngləri birbaşa nümunələnir,
beləliklə iki motor arasında keçiddə rəng sıçrayışı olmur.

## Yol boyu tapılan iki səhv

1. **Kənar düyünlərdə çıxıntı** — ilk versiyada düyün koordinatları
   `min(i, nx-1)` ilə sıxılırdı, bu, modelin kənarında hüceyrə eninə
   bərabər artıq çıxıntı yaradırdı (ilk render sınağında göründü).
   Vektorlaşdırılıb düzəldildi — həm doğru, həm sürətli.
2. **`unhashable type: LinearSegmentedColormap`** — interfeys rəng
   ADI yox, matplotlib OBYEKTİ ötürür; ilk versiya yalnız ad
   gözləyirdi.

## v54 — bütün funksiyalar köçürüldü

İstifadəçi istədi: "matplotlib-də hansı funksiyalar varsa tam olaraq
VTK-ya köçür". Köçürülənlər:

| Funksiya | Vəziyyət |
|---|---|
| Xassə seçimi, rəng xəritəsi, legend | ✅ |
| Kəsim həddi (dəyər filtri) | ✅ (`BlankCell`) |
| K-təbəqə diapazonu | ✅ (`BlankCell`) |
| Şəffaflıq | ✅ |
| Kənarlar (aç/bağla) | ✅ |
| Baxış bucaqları (4 hazır) | ✅ |
| Yaxınlaşdırma | ✅ |
| İşıq sürgüsü | ✅ (əsl material xassələri) |
| Quyular (lülə + perforasiya + ad) | ✅ |
| Faultlar | ✅ |
| Zaman slider-i | ✅ |

### Quyular — matplotlib-dən daha yaxşı

matplotlib-də lülə sadəcə qalın XƏTT idi. VTK-da əsl SİLİNDR
(`vtkCylinderSource`), perforasiyalar isə kürələr — model fırlananda
həqiqi 3D obyekt kimi davranır, dərinlik hissi düzgün verilir.

Ad üçün `vtkBillboardTextActor3D` işlədilir: mətn model fırlananda
HƏMİŞƏ kameraya baxır (adi 3D mətn fırlanıb oxunmaz olardı).

**Tapılan səhv:** ilk versiyada lülə tam hüceyrələrin içində
gizlənirdi — yalnız adlar görünürdü. Lülə indi modelin səthindən
yuxarıda başlayır (ümumi qalınlığın 25 %-i qədər), beləliklə istənilən
qalınlıqda nisbətli görünür.

### Faultlar — VTK-da ƏSL DOLĞU

Bu, VTK-ya keçidin ən aydın qazancıdır. matplotlib motorunda fault
üçün dolğu İŞLƏTMƏK MÜMKÜN OLMAMIŞDI — dərinlik sıralaması dəqiq
olmadığı üçün yarı-şəffaf müstəvi "cırıq" görünürdü, ona görə orada
yalnız kontur çəkilir (bax `FAULTS.md`).

VTK-da əsl z-buffer var: fault indi yarı-şəffaf müstəvi kimi düzgün
görünür. Şəffaflıq çarpandan asılıdır (sealing fault daha
qeyri-şəffaf) — istifadəçi cədvələ baxmadan hansı faultun axını nə
qədər bloklandığını görür.

### İşıq

matplotlib-də kölgə ƏL İLƏ hesablanırdı (üzün istiqamətinə görə
rəng qaraldılırdı). VTK-da bu, materialın əsl işıq xassələridir —
`Ambient`/`Diffuse`/`Specular` nisbəti sürgü ilə tənzimlənir.

## v55 — koordinat şəbəkəsi və istiqamət oxu

İstifadəçi ResInsight-dakı iki elementi göstərib istədi.

### Koordinat şəbəkəsi (`vtkCubeAxesActor`)

Modelin ətrafında X/Y/Z ölçü oxları, rəqəmlərlə. Bu, modelə MİQYAS
HİSSİ verir — onsuz model "havada asılı" görünürdü.

**Dərinlik oxunda iki səhv tapıldı və düzəldildi:**

1. **Avtomatik üstlü miqyaslama** — VTK defolt olaraq etiketləri
   "×10³" ilə sıxır. 2000–2032 m aralığı bu formatda "2 2 2 2" kimi
   görünürdü (bütün fərqlər itirdi). `SetLabelScaling(False, ...)`
   ilə söndürüldü.
2. **Format** — `SetZLabelFormat("%.0f")`, çünki VTK-nın defoltu
   ("%-#6.3g") 3 əhəmiyyətli rəqəmlə bu aralığı ayırd edə bilmir.

Əlavə olaraq: nazik rezervuarlarda (qalınlıq / areal ölçü < 0.25)
dərinlik etiketləri avtomatik kiçildilir, kiçik bölgülər söndürülür —
sıxlıq problemini azaldır.

**Şaquli mübaliğə** etiketlərə TƏSİR ETMİR: mübaliğə yalnız görüntünü
uzadır, ox aralığı əsl dərinliyə bölünərək verilir, istifadəçi həmişə
real metr görür (test bunu qoruyur).

### İstiqamət oxu (`vtkOrientationMarkerWidget`)

Sağ aşağı küncdə kiçik X/Y/Z işarəsi — model fırlananda o da fırlanır.

Bu, adi aktyor DEYİL (öz kiçik görüntü sahəsi var) və İNTERAKTOR
tələb edir — ona görə səhnə qurulanda yox, interfeys interaktoru
hazırlayandan sonra bağlanır. Offscreen render zamanı (test, şəkil
saxlama) interaktor olmadığı üçün səssizcə atlanılır.

## Qalan iş

- Kəsim müstəviləri (ResInsight-ın güclü tərəfi — VTK-da `vtkPlane`
  ilə asandır)
- Zaman üzrə avtomatik animasiya

## v56 — koordinat şəbəkəsi və istiqamət oxu

İstifadəçi ResInsight-ın iki elementini istədi (şəkil üzərində
işarələdi): modelin ətrafındakı **ölçü şəbəkəsi** və sağ aşağı
küncdəki **istiqamət oxu**.

### Koordinat şəbəkəsi (`vtkCubeAxesActor`)

Modelə MİQYAS HİSSİ verir — onsuz model "havada asılı" görünürdü.
X/Y oxları metrlə, dərinlik oxu ayrıca.

**Dərinlik etiketləri üçün üç yanaşma sınandı:**

1. Kiçik font — etiketlər hələ də üst-üstə düşürdü
2. `SetScreenSize` — VTK-nın etiket SAYINA təsir etmədi
   (`vtkCubeAxesActor`-da etiket sayını birbaşa idarə edən metod
   YOXDUR, yoxlanılıb)
3. **NİSBİ dərinlik** — işlədi

Səbəb: mütləq dərinlik (2000, 2004, 2008…) 4 rəqəmlidir, aralıq isə
dardır (onlarla metr) — etiketlər sıxılıb oxunmaz olur. Nisbi
dərinlik (0, 5, 10… m — tavandan aşağı) rəqəmləri 1-2 rəqəmə endirir.
Mütləq dərinlik oxun BAŞLIĞINDA saxlanılır ("tavan 2000 m-dən"), ona
görə məlumat itmir.

Şaquli mübaliğə oxa TƏSİR ETMİR — rəqəmlər həmişə əsl metrləri
göstərir (test bunu qoruyur).

### İstiqamət oxu (`vtkOrientationMarkerWidget`)

Sağ aşağı küncdə X/Y/Z işarəsi, model ilə birlikdə fırlanır.

Bu, adi aktyor DEYİL — öz kiçik görüntü sahəsi var və İNTERAKTOR
tələb edir. Ona görə səhnə qurulanda yox, interfeys interaktoru
hazırlayandan sonra qoşulur; offscreen render zamanı (test, şəkil
saxlama) sadəcə atlanılır.

## v57 — oxlar sadələşdirildi (istifadəçi rəyi)

İstifadəçi v56-nın nəticəsini bəyənmədi və şəkil göndərdi: yazılar
NƏHƏNG idi (bütün ekranı tuturdu), etiketlər üst-üstə düşürdü, "ə"
hərfi itirdi ("Dərinlik" → "Dinlik").

### Üç düzəliş

1. **Başlıqlar tamamilə SİLİNDİ.** İstifadəçi belə istədi: *"yalnız
   rəqəmlər qalsın, sağ aşağıdan onsuz bilirik hansı X, hansı Y
   oxudur"*. Bu, eyni zamanda "ə" hərfi probleminin KÖKÜNÜ kəsdi —
   VTK-nın defolt fontu Azərbaycan hərflərini dəstəkləmir.

2. **Etiket ölçüsü kiçildildi** (`SetScreenSize` 9-11, font 10).
   VTK mətni MODEL ölçüsünə görə miqyaslayır — böyük modeldə
   rəqəmlər nəhəng olurdu.

3. **Nazik modeldə dərinlik etiketləri gizlədilir.** Qalınlıq areal
   ölçünün 12 %-indən azdırsa, dərinlik oxu ekranda o qədər qısa olur
   ki, rəqəmlər bir-birinin üstünə yığılır. Hədd ölçülüb: 41×41×5
   (820 m areal, 50 m qalınlıq, nisbət 0.06) hələ sıx görünürdü.

Nəticə iki halda yoxlanıldı: nazik model (41×41×5) — tam təmiz,
yalnız X/Y rəqəmləri; qalın model (14×12×6) — dərinlik oxu (0–140 m)
aydın oxunur.

## v58 — VTK xəta pəncərəsi düzəldildi (v57-nin yan təsiri)

v57-də başlıqları gizlətmək üçün onlara BOŞ SƏTİR (`""`) verilmişdi.
Nəticə: model düzgün çəkilirdi, LAKİN VTK-nın daxili `vtkVectorText`
filtri boş mətni qəbul etmir — hər kadrda `"Text is not set!"` xətası
atırdı və Windows-da fasiləsiz `vtkOutputWindow` pəncərəsi açılırdı
(istifadəçi şəkildə göstərdi).

**Düzəliş:** boşluq simvolu (`" "`) — VTK üçün keçərli mətndir,
ekranda isə heç nə görünmür.

**Ölçüldü:** boş sətir 3 kadrda 3 xəta verir, boşluq simvolu 0 xəta.

Bu səhvin xüsusiyyəti: kod "işləyirdi", şəkil düzgün çıxırdı, testlər
keçirdi — problem yalnız istifadəçinin ekranında görünürdü. Ona görə
`test_rendering_produces_no_vtk_errors` testi əlavə olundu: render
zamanı VTK-nın xəta çıxışı tutulur və boş olması yoxlanılır.
