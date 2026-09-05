# Eclipse mübadiləsi (B5)

Uzunmüddətli alət təcrid olunmuş qala bilməz: geoloq Petrel-də model
qurur, mühəndis CMG və ya Eclipse-də işlədir. Bu modul həmin zəncirə
qoşulmanı təmin edir.

`Layihə` menyusunda:

| Bənd | Nə edir |
|---|---|
| **GRDECL grid oxu…** | Xarici grid faylını geoloji model kimi yükləyir |
| **Eclipse deck yaz (.DATA)…** | Modeli Eclipse formatında ixrac edir |

## Oxuma

Dəstəklənən açar sözlər:

| Qrup | Açar sözlər |
|---|---|
| Ölçü | `SPECGRID`, `DIMENS` |
| Həndəsə | `DX`, `DY`, `DZ`, `TOPS` · `COORD`, `ZCORN` (approksimasiya) |
| Xassələr | `PORO`, `PERMX`, `PERMY`, `PERMZ`, `NTG` |
| Regionlar | `SATNUM`, `PVTNUM`, `EQLNUM`, `FIPNUM` |
| Digər | `ACTNUM` (**tətbiq olunur** — bax aşağı), `MULTX/Y/Z` |

### Təkrar sintaksisi

    PORO
      8405*0.22 /

`n*value` — dəyər n dəfə. Bu, real fayllarda mütləqdir: 200 000
hüceyrəli modeldə hər dəyəri ayrıca yazmaq faylı yüzlərlə meqabayt
edərdi. Yazıcı da eyni sıxılmanı tətbiq edir.

### Tapılan incə səhv

İlk versiyada bölmə başlıqları (`RUNSPEC`, `GRID`, `PROPS`) adi açar
söz kimi qəbul olunurdu. Onlar isə Eclipse-də `/` TƏLƏB ETMİR — nəticədə
oxucu növbəti `/` işarəsinə qədər hər şeyi udurdu və **ilk massiv
(adətən `DX`) itirdi**.

Model səssizcə defolt 100 m hüceyrə ölçüsü ilə qurulurdu. Dövrə testi
bunu üzə çıxardı: RF 37.4 % əvəzinə 5.6 %. İndi `test_section_headers_
do_not_swallow_the_next_array` bunu qoruyur.

## Dürüst məhdudiyyətlər

Bunlar gizlədilmir — hər biri istifadəçiyə xəbərdarlıq kimi göstərilir:

| Məhdudiyyət | Davranış |
|---|---|
| **Corner-point həndəsə** | Oxunur, bərabər bloka **approksimasiya** olunur. Orta `DX/DY/DZ` hesablanır, istifadəçi xəbərdarlıq alır |
| **Dəyişkən hüceyrə ölçüsü** | Orta qiymət götürülür + xəbərdarlıq |
| **`ACTNUM`** | **Tam dəstəklənir.** Qeyri-aktiv hüceyrə simulyasiyadan çıxarılır: `PV = 0`, qonşuluq bağlantısı qurulmur, xətti sistem `n_active` naməlumla həll olunur, həmin hüceyrədəki perforasiya söndürülür (`WI = 0`). Bax `imex2d/domain/grid.py::ActiveMap` və `tests/test_actnum.py` |
| **`INCLUDE`** | Dəstəklənmir; həmin fayl oxunmur — xəbərdarlıq verilir |

Səbəb: `CellGeometry` hazırda yalnız bərabər ölçülü bloklar saxlayır
(bax `ARCHITECTURE.md`, 5.1). Corner-point dəstəyi həmin sinfin yeni
implementasiyasını tələb edir — diskretizasiya kodu dəyişməyəcək,
çünki o, `face_areas()` və `face_half_distances()` interfeysindən
istifadə edir.

## Yazma

`.DATA` deck-i altı bölmədən ibarətdir: `RUNSPEC`, `GRID`, `PROPS`,
`SOLUTION`, `SUMMARY`, `SCHEDULE`.

Yazılanlar: grid həndəsəsi, `PORO`/`PERMX`/`PERMY`/`PERMZ`/`NTG`,
`SATNUM`, `SWOF` cədvəli (Corey əyriləri + kapilyar təzyiq), PVT
(`PVDO`/`PVTW` və ya `PVCDO`), `EQUIL`, quyular (`WELSPECS`,
`COMPDAT`, `WCONPROD`, `WCONINJE`), `TSTEP`.

**Vahidlər: METRIC** — modelin öz sistemi ilə eynidir, ona görə heç
bir çevirmə aparılmır və çevirmə səhvi riski yoxdur.

Deck avtomatik yaradılır və tam işlək olmaya bilər — hər simulyatorun
öz xüsusiyyətləri var. Struktur düzgündür və əl ilə tamamlanmağa
hazırdır.

## Dövrə testi

`test_round_trip_reproduces_the_simulation_result` — modeli yazır,
geri oxuyur, hər ikisini işə salır və nəticələri müqayisə edir:

| Göstərici | Nəticə |
|---|---|
| Həndəsə (`dx`, `dy`, `dz`, dərinlik) | dəqiq |
| `PORO` | 1e-5 daxilində |
| `PERMX` | nisbi 1e-4 daxilində |
| OOIP | nisbi 1e-4 daxilində |
| Recovery Factor | 0.01 % daxilində |

Fərq yalnız fayla yazılan onluq dəqiqlikdən gəlir.

## ACTNUM — qeyri-aktiv hüceyrələr

`ACTNUM` deck-də varsa, o, **grid topologiyasının bir hissəsi** olur
(`CartesianGrid.actnum` → `ActiveMap`) və idxal ANINDA simulyasiya
modelini reduksiya edir — dörd qatda eyni vaxtda:

| Qat | Davranış |
|---|---|
| Topologiya (`grid.build_connections`) | Aktiv↔qeyri-aktiv üz **ümumiyyətlə qurulmur**. `T = 0` yazmaq kifayət DEYİL — sıfır element seyrək matrisin strukturunda qalır və hüceyrəni sistemə naməlum kimi gətirir |
| Həcm (`ReservoirModel.pore_volume`) | Qeyri-aktiv hüceyrədə `PV = 0` → OOIP, akkumulyasiya, material balansı və CFL addımı avtomatik düzgün çıxır |
| Xətti sistem (IMPES + Nyuton) | Matris `n_active` (Nyutonda `2·n_active`) ölçüsündədir; qeyri-aktiv hüceyrənin NAMƏLUMU YOXDUR |
| Quyular (`PeacemanWellModel`) | Qeyri-aktiv hüceyrədəki perforasiya bağlantı siyahısına salınmır (`WI = 0`) + diaqnostika xəbərdarlığı. Quyunun BÜTÜN perforasiyaları qeyri-aktivdirsə — XƏTA |

Massivlərin **saxlama** formatı qlobal (`ncell`) qalır: 3D görüntü,
hesabat, serializasiya və fayl formatlarının hamısı qlobal indeksləmə
üzərində qurulub. Reduksiya YALNIZ xətti sistem sərhədində baş verir.

```python
grid.n_active                  # simulyasiyadakı naməlum sayı
grid.active.global_to_active   # (ncell,)    qeyri-aktiv üçün -1
grid.active.active_to_global   # (n_active,)
```

Kənar hallar:

* bütün hüceyrələr qeyri-aktivdirsə idxal `GrdeclError` ilə dayanır;
* `ACTNUM`-un ölçüsü grid ilə uyğun gəlmirsə xəbərdarlıq verilir və
  massiv NƏZƏRƏ ALINMIR (səssiz sürüşmə yoxdur);
* `NaN` (deck-in `n*` defoltu) **aktiv** sayılır — oxunmamış hüceyrəni
  səssizcə yox etmək təhlükəlidir;
* MPFA-O bu modelləri HƏLƏ həll edə bilmir və bunu AÇIQ bildirir
  (`unsupported_features`); TPFA tam dəstəkləyir.

Testlər: `tests/test_actnum.py` (o cümlədən qeyri-aktiv "doldurucu"
hüceyrələrin nəticəni DƏYİŞMƏDİYİNİ sübut edən bərabərlik testi).
