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
| Digər | `ACTNUM` (xəbərdarlıqla), `MULTX/Y/Z` |

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
| **`ACTNUM`** | Qeyri-aktiv hüceyrələr **aktiv sayılır**; həcm hesabı böyük çıxır — xəbərdarlıq verilir |
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
