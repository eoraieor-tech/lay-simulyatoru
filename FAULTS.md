# Fault transmissivliyi (B3)

## Fikir

Fault iki hüceyrə blokunu ayıran müstəvidir. Fiziki mənası: fay
zonasında gil yaxılması (gouge) axını qismən və ya tam bloklaya
bilər. Modeldə bu, sadəcə həmin müstəvidəki üzlərin transmissivliyinin
bir çarpanla azaldılması kimi ifadə olunur — Eclipse-in
`FAULTS`/`MULTFLT` cütlüyü ilə eyni fikir.

    çarpan = 1.0   →  fay şəffafdır (adi hüceyrə sərhədi kimi)
    çarpan = 0.0   →  fay tam bağlıdır (axın yoxdur, sealing)
    0 < çarpan < 1 →  qismən keçirici

## Ad — "Fay" yox, "Fault"

İnterfeysdə və bütün istifadəçiyə görünən mətnlərdə ingiliscə
`Fault` termini işlədilir (Eclipse/CMG-dəki adlandırma ilə uyğun) —
`fay` yalnız kod şərhlərində qalır.

## İki subyekt, bir bişmiş nəticə

Geoloq və mühəndis fərqli məlumatla işləyir:

- **Geoloq** `GeologicalModel.Fault`-u yaradır — hansı I/J/K
  müstəvisində, hansı diapazonda. `polyline`/`throw`/`dip` yalnız
  təsviridir, xəritədə göstərmək üçün; hesablamaya təsir etmir.
- **Mühəndis** simulyasiya üçün yalnız transmissivlik çarpanını
  dəyişir, geoloji detalları bilməli deyil.

`FaultReference` bu ikisinin **bişmiş** nəticəsidir: modelin qurulma
anında `GeologicalModel.fault_references()` vasitəsilə həndəsə
köçürülür. Diskretizasiya birbaşa bunu oxuyur.

## Niyə müstəvi, əyri səth yox

`Fault.polyline` ümumi 3D əyri ola bilərdi, amma hazırkı `CellGeometry`
yalnız bərabər ölçülü bloklar saxlayır (bax `ARCHITECTURE.md`, 5.1).
Əyri fay-hüceyrə kəsişməsini hesablamaq corner-point həndəsə tələb
edir ki, bu hələ yoxdur.

Ona görə fay **grid oxlarına düz bucaqlı müstəvi** kimi təyin olunur —
bu, həm cari həndəsə ilə tam uyğundur, həm də demək olar bütün
kommersiya simulyatorlarının (Eclipse daxil olmaqla, `FAULTS` açar
sözü) əslində işlətdiyi yanaşmadır.

## Tərif

    I1 I2  J1 J2  K1 K2  ÜZ     (Eclipse sintaksisi, 1-based)

Daxili təmsil 0-based və sadədir: `axis` (I/J/K), `plane_index`
(sərhəd — `plane_index` və `plane_index+1` arasında), `range_a`,
`range_b` (digər iki ox üzrə diapazon, `None` = bütün grid).

## Tətbiq

`TwoPointFluxDiscretization._apply_fault_multipliers()` — hər fay
üçün: hansı bağlantılar (`Connections`) onun müstəvisinə düşür,
tapılır, transmissivlik çarpanla vurulur.

Bir üz birdən çox faya düşərsə, çarpanlar **vurulur** — iki qismən
keçirici fay üst-üstə düşəndə axının daha da azalması fiziki
cəhətdən doğrudur.

## Ölçülmüş fiziki doğruluq

20×10 grid, i=10/11 sərhədində tam en boyu fay, `RATE` idarə olunan
vurucu bir tərəfdə, istismarçı digər tərəfdə:

| Fay | RF (400 gün) |
|---|---|
| Yoxdur | 29.7 % |
| Qismən (çarpan 0.05) | 20.4 % |
| Tam bağlı (sealing) | **0.9 %** |

Tam bağlı fay istismarçını təzyiq dəstəyindən demək olar tamamilə
kəsir — gözlənilən nəticə, çünki grid iki əlaqəsiz hissəyə bölünür.

## Diaqnostika

- Grid-dən kənar müstəvi/diapazon → **XƏTA** (model işə salınmır)
- Eyni adlı fay təkrarlanır → **XƏBƏRDARLIQ** (çarpanlar vurulacaq)

## Fayl formatları

**CSV:**

    name,axis,plane_index,a_low,a_high,b_low,b_high,multiplier,sealing
    F1,I,10,0,40,,,0.1,0

Diapazon sütunları boş buraxılarsa, həmin ox üzrə bütün grid əhatə
olunur.

**Eclipse (`FAULTS` + `MULTFLT`):**

    FAULTS
      'F1'  11 11  1 41  1 5  'I' /
    /
    MULTFLT
      'F1'  0.1 /
    /

`MULTFLT` verilməyibsə çarpan 1.0 (şəffaf, yalnız qeydiyyat) qəbul
edilir.

## İnterfeys

`4 · FAULTS` bölməsi: CSV yüklə, Eclipse FAULTS yüklə, ya da əl ilə
əlavə et (pəncərə ilə ad/ox/müstəvi/çarpan/sealing). Cədvəldə bütün
fault-lar görünür.

## 3D görüntüdə

Hər fault öz tam diapazonunun HƏDD qutusunu (bounding box) əhatə
edən **qırmızı konturla** çəkilir — dolğu yox, çünki matplotlib-in
3D dərinlik-sıralaması iri, yarı-şəffaf poliqonu hüceyrə üzləri ilə
kəsişəndə "cırıq"/qara zolaq artefaktı yaradır (mplot3d-nin məlum
məhdudiyyəti). Kontur bu problemə məruz qalmır.

Xəttin qalınlığı çarpandan asılıdır: tam bağlı (sealing) fault qalın,
şəffaf fault nazik xətlə çəkilir — istifadəçi cədvələ baxmadan hansı
faultun axını nə qədər məhdudlaşdırdığını görür.

`3 D görüntü` tabında **"Faults"** qutusu ilə göstərilib-gizlədilə
bilər (quyular kimi).
