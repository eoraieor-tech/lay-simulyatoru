# OPM Flow nəticələrinin idxalı — strateji dönüş

> **Status qeydi (yeniləndi: 1 oktyabr 2026, Seans 51).** Aşağıda təsvir
> olunan «fizikanı OPM Flow-a həvalə etmək» qərarı (v69) **10 sentyabr 2026-da
> ləğv olundu**: sahibkar «öz fizikamız» seçdi, üç fazalı (qaz) həlledici
> bərpa olunub servisə və UI-yə qaytarıldı (`ISH_HESABATI.md` → Seans 2, 6;
> `ROADMAP.md` 0.10). Yəni hazırkı mühərrik iki VƏ üç fazalıdır.
>
> OPM idxalının özü (`imex2d/io/opm_import.py`) müstəqil imkan kimi qalır və
> istəyə bağlı `resdata` paketini tələb edir (`requirements.txt`-də yoxdur;
> paket olmayanda `tests/test_opm_import.py` atlanır). OPM Flow indi başqa
> rolda işlədilir — SPE1CASE2 üçün **etalon** kimi; onun summary faylları
> `resdata`-sız öz oxuyucumuzla oxunur (`io/eclipse_summary.py`, bax `SPE1.md`).

## Fikir

A7-nin öz üç fazalı Nyuton həlledicisi açıq bir davamlılıq problemi
daşıyırdı (bax `A7_PLAN.md`) — quyu öz BHP sərhədinə çox yaxınlaşanda
bəzən yığılmırdı. Həmin həlledici v69-da kod bazasından çıxarıldı
(sonradan qaytarıldı — yuxarıdakı status qeydinə bax). Aşağıdakı OPM
idxalı o tarixi qərarın nəticəsi olaraq qalıb və müstəqil şəkildə
faydalıdır; "qaz fazasının OPM-dən idxalı" ideyası (aşağıda, "Növbəti
addımlar") heç vaxt tətbiq olunmayıb və hazırda planlaşdırılmır.

Bunun əvəzinə strateji qərar verildi: **fizikanı** real, sınanmış bir
simulyatora (OPM Flow — açıq mənbəli, Eclipse formatına uyğun) həvalə
edirik, öz proqramımızın güclü tərəfini (3D görüntü, analiz) OPM-in
NƏTİCƏLƏRİNİ göstərmək üçün işlədirik.

Planlaşdırılan iki mərhələ:
1. **Bu mərhələ:** OPM Flow-un 3D görüntüsünü öz görüntümüzlə əvəz et
2. **Növbəti:** uğurlu olsa, OPM Flow-dan qaz fazasını (real, sınanmış
   fizika ilə) əlavə et

## Texniki yanaşma

`resdata` kitabxanası (libecl-in müasir varisi) Eclipse binar
formatını oxuyur — OPM Flow-un çıxışı da bu formatdadır (`.EGRID`
grid, `.UNRST` hesabat addımları).

### Təhlükəsizlik qeydi — VACİB

resdata-nın "yüksək səviyyəli" rahatlıq metodları (`iget_restart_sim_time`
və s.) natamam/qeyri-standart fayllarda **native SIGABRT** ilə çökür —
bunu Python-un `try/except`-i TUTA BİLMİR (proses bütöv dayanır).
Kəşf yolu ilə tapıldı: sintetik test faylı ilə sınayarkən bu native
çökmə birbaşa müşahidə olundu.

Həll: `imex2d/io/opm_import.py` YALNIZ aşağı səviyyəli, təhlükəsiz
metodları işlədir (`num_named_kw()`, `iget_named_kw()`) və zamanı
`DOUBHEAD`-in 0-cı sahəsindən (elapsed sim time, Eclipse
sənədləşməsinə görə bütün versiyalarda sabit) birbaşa oxuyur.

### Model uyğunlaşdırması

`build_display_model()` OPM halını bizim öz `ReservoirModel`-imizə
çevirir — **ayrıca adapter sinfi yaratmaq əvəzinə** birbaşa öz
`CartesianGrid`/`CellGeometry` siniflərimizdən istifadə edir (bunlar
artıq düzbucaqlı həndəsəni tam əhatə edir). Bu sayədə mövcud
`VolumeRenderer` heç bir dəyişiklik olmadan işləyir.

### Sadələşdirmə

Yalnız **düzbucaqlı** (dx/dy/dz sabit) torlar dəstəklənir.
Corner-point (mürəkkəb həndəsə) torlar hələ dəstəklənmir — açıq
xəbərdarlıq verilir, səssizcə yanlış göstərmək əvəzinə.

## UI inteqrasiyası

`Layihə → OPM Flow nəticəsini yüklə (.EGRID+.UNRST)…`

Diqqətəlayiq dizayn qərarı: heç bir yeni widget yaradılmadı. Mövcud
"3D görüntü" tabındakı bütün idarəetmə (vaxt slider-i, xassə seçimi,
K-təbəqə filtri) artıq `self.reservoir_model`/`self.result.snapshots`
üzərində işləyir — idxal əməliyyatı sadəcə bunları doldurur, sonra
`update_volume()` çağırılır. Tam təkrar istifadə.

## Doğrulama

12 test (`test_opm_import.py`, pytest — `tmp_path` fixture-i lazımdır):
sintetik round-trip Eclipse halı (grid + 3 hesabat addımı, PRESSURE/
SWAT/SGAS) yaradılır, oxunur, dəyərlər yoxlanılır, VƏ öz
`VolumeRenderer`-imizlə çəkilməsi (çökmədən) təsdiqlənir — bu, əsl
inteqrasiya məqsədini (OPM göstəricisini əvəz etmək) sınayır.

## Növbəti addımlar

- Corner-point (mürəkkəb həndəsə) tor dəstəyi
- Quyu mövqelərinin idxalı (hazırda göstərilmir)
- Zaman üzrə animasiya (OPM-in bəyənilən xüsusiyyəti)
- Rəng xəritəsi/legend-in OPM-ə bənzədilməsi (VTK-ya keçid — əvvəlki plan)
- Uğurlu olduqdan sonra: OPM-dən qaz fazası nəticələrinin idxalı
