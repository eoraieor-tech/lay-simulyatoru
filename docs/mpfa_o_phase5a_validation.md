# Phase 5A — MPFA-O validasiya hesabatı

Tapşırıq §34-də tələb olunan A–P bölmələri + §35 yekun siyahısı.
Riyazi spesifikasiya: [`mpfa_o_phase5a.md`](mpfa_o_phase5a.md).
Bütün rəqəmlər `tests/test_mpfa_o.py` və onun köməkçiləri ilə
təkrar-istehsal olunandır.

---

## A. Dəyişən fayllar

**YENİ**

| fayl | məzmun |
|---|---|
| `docs/mpfa_o_phase5a.md` | riyazi spesifikasiya (§1–§17) |
| `docs/mpfa_o_phase5a_validation.md` | bu hesabat |
| `imex2d/discretization/__init__.py` | paket girişi |
| `imex2d/discretization/mpfa_o_interaction.py` | bölgələr, sub-üzlər, kəsilməzlik nöqtələri |
| `imex2d/discretization/mpfa_o_local_system.py` | lokal sistem, sub-cell qradiyentləri, diaqnostika |
| `imex2d/discretization/mpfa_o.py` | qlobal əmsal yığımı, `MPFAODiscretization` |
| `tests/test_mpfa_o.py` | 67 test |

**DƏYİŞDİRİLMİŞ (yalnız ƏLAVƏ, davranış pozulmadan)**

| fayl | dəyişiklik |
|---|---|
| `imex2d/domain/general_grid_geometry.py` | `face_index(cell, local_name)` O(1) axtarışı + `connection_faces()`; mövcud API/nəticələr toxunulmayıb |
| `imex2d/interfaces/discretization.py` | NON-ABSTRAKT `supports_multipoint_stencil()` (defolt `False`) |
| `ARCHITECTURE.md` | §5.16 MPFA-O bölməsi |

**TOXUNULMAYIB**: `imex2d/simulation/discretization.py` (TPFA),
`implicit/residual.py`, `implicit/jacobian.py`, `implicit/newton.py`,
PVT/SCAL/quyu/UI/serializasiya — heç biri.

---

## B. Riyazi formulyasiya

Sub-cell `(c,v)` daxilində xətti təzyiq və 3 kəsilməzlik şərti:

```
p(x) = p_c + g_(c,v)·(x − x_c)
D_(c,v) g_(c,v) = π_S − p_c·1        D sətirləri: (x_σ − x_c)ᵀ
x_σ = x_v + η(x_F − x_v)              η = 1 (defolt)
```

Yarım-axın (TAM tenzor, `nᵀKn` skalyarlaşdırma YOXDUR):

```
q_(c,σ) = −Γ a_σ^(c)ᵀ K_c D_(c,v)⁻¹ (π_S − p_c·1)
W_(c)   = Γ · A_c · K_c · D_(c,v)⁻¹        (3×3)
```

Lokal sistem və axın bərpası:

```
C π_unk = D p + E π_bnd
q_sub   = (F C⁻¹ D + G) p + (F C⁻¹ E + H) π_bnd
        =      T_cell    p +      T_bnd    π_bnd
q_F     = Σ_{v ∈ 4 künc} q_(F,v)
```

## C. Qarşılıqlı təsir bölgələri

Hər grid TƏPƏSİ üçün bir bölgə. Struktur topologiyadan
DETERMİNİSTİK qurulur (`(di,dj,dk)` oktantı → yerli təpə indeksi →
həmin təpəni saxlayan 3 yerli üz → `GeneralGridGeometry.face_index`).
**Koordinat müqayisəsi/floating-point axtarış YOXDUR** → `O(N)`.

Ölçülmüş struktur (3×3×2 grid, 48 bölgə):

| bölgə tipi | hüceyrə | sub-üz | naməlum π (Dirichlet) |
|---|---|---|---|
| daxili | 8 | 12 | 12 |
| sərhəd (üz) | 4 | 8..10 | 4..6 |
| sərhəd (künc) | 1 | 3 | 0 |

Sub-üz `σ=(F,v)` poliqonu `[v, mid(v,v₊), x_F, mid(v₋,v)]`; tiling
xassəsi `Σ_q a_(F,u_q) = A_F n_F` **testlə** yoxlanılır
(`test_reference_case_sub_face_area_vectors_sum_to_the_face_area_vector`,
skew gridd-də `rtol=1e-10`).

## D. Lokal sistemlər

* **naməlumlar**: daxili sub-üzlərin `π_σ`-ları (Neumann bağlanışında
  sərhəd sub-üzləri də);
* **sətirlər**: hər naməlum üçün bir bağlanış tənliyi — daxili sub-üzdə
  `q_o + q_n = 0`, `NEUMANN_ZERO`-da `q_o = 0`;
* **sütunlar**: `C` → naməlum π, `D` → bölgə hüceyrə təzyiqləri,
  `E` → xaricdən verilən sərhəd π;
* **axın bərpası**: `F/G/H` sətirləri OWNER tərəfindən qurulur → sub-üz
  səviyyəsində antisimmetriya İDENTİKDİR.

Bütün matrislər (`C,D,E,F,G,H,T_cell,T_bnd`) və hər sub-cell-in
`D_(c,v)`, `D⁻¹`, `W` matrisləri `MPFAOLocalSystem`-də AÇIQ
atributlardır; `describe()` tam dump verir.

### Əl ilə yoxlanıla bilən istinad halı (§30)

2×2×2 vahid kub, `K = I`, `Γ = 1`, `η = 1`, mərkəzi node `(1,1,1)`:

```
hüceyrə mərkəzləri : (0.5|1.5, 0.5|1.5, 0.5|1.5)
σ0                 : üz 8, owner 0, neighbor 4
  A_σ = 0.25        a_σ = (0, 0, 0.25)      x_σ = (0.5, 0.5, 1.0)
D_(c=0,v) = [[0, 0, 0.5], [0, 0.5, 0], [0.5, 0, 0]]
W_(c=0)   = 0.5 · I
C: 12×12    D: 12×8    T_cell: 12×8
T_cell[0] = [+0.25, 0, 0, 0, −0.25, 0, 0, 0]
```

Əl hesabı: `T_σ = Γ A_σ / (h_o/k + h_n/k) = 1·0.25/(0.5+0.5) = 0.25` ✔

`p = x` sahəsində: `q_σ = 0.25·(p_o − p_n)`; konservasiya qalığı
`3.3e-16`.

## E. Tenzorun işlənməsi — 6 komponentin HAMISI təsir edir

Skew grid, baza `K = diag(500, 50, 10)`; hər komponent ayrıca dəyişdirilib
(diaqonal ×1.1, off-diaqonal 0 → 15) və `T_cell`-in NİSBİ dəyişməsi:

| komponent | nisbi əmsal dəyişməsi | nəticə |
|---|---|---|
| Kxx | 9.818e-02 | TƏSİR EDİR |
| Kyy | 1.062e-02 | TƏSİR EDİR |
| Kzz | 7.959e-03 | TƏSİR EDİR |
| **Kxy** | **2.508e-02** | **TƏSİR EDİR** |
| **Kxz** | **2.855e-02** | **TƏSİR EDİR** |
| **Kyz** | **2.879e-02** | **TƏSİR EDİR** |

Test: `test_g_off_diagonal_components_change_the_coefficients`.

## F. Çoxnöqtəli sübut

3×3×2 skew grid, `K = R diag(500,50,10) Rᵀ` (37°, ox `(0.3,−0.5,0.8)`):

```
daxili üz stensil ölçüləri: {8: 20 üz, 12: 12 üz, 18: 1 üz}
```

Ən geniş üz (üz 28, owner 4, neighbor 13) — **18 hüceyrə**:

```
T[28, 4] = +5.645548     T[28,13] = −5.645548
T[28, 3] = +2.899346     T[28,14] = −2.899346
T[28, 1] = +2.404583     T[28,16] = −2.404583
T[28, 0] = +0.963491     T[28,17] = −0.963491
... (qalan 10 hüceyrə)
```

**Anti-pseudo sübutu**: `p_owner = p_neighbor` qurulduqda
`q[28] = 121.75` (axın miqyası 4203.9). İki-nöqtəli sxem burada
MƏCBURİ `0` verərdi. Test:
`test_anti_pseudo_mpfa_flux_is_not_a_two_point_relation`.

## G. Konservasiya

`max_σ |q_{o,σ} + q_{n,σ}|` — hər tərəf ÖZ `g_(c,v)` qradiyentindən
MÜSTƏQİL hesablanır (tavtologiya deyil):

| konfiqurasiya | konservasiya qalığı | hüceyrə üzrə `max abs div` |
|---|---|---|
| Kartezian + izotrop | 1.705e-13 | 8.669e-13 |
| Skew + izotrop | 2.842e-13 | 8.669e-13 |
| Skew + fırlanmış tam tenzor | 1.776e-13 | 6.253e-13 |
| Warped + fırlanmış tam tenzor | 3.126e-13 | 6.395e-13 |

Təsadüfi (fiziki olmayan) təzyiq və sərhəd dəyərləri ilə də eyni
səviyyə (`test_i_local_conservation_of_sub_face_fluxes`).

## H. Manufactured solution: ədədi vs analitik

`p = 2x − 3y + 1.25z + 11`, sabit `K` → `u = −K∇p` analitik olaraq
məlumdur. Sərhəd `π` dəyərləri analitik həldən verilir (Dirichlet).

| konfiqurasiya | `max|q_ədədi − q_analitik|` | nisbi |
|---|---|---|
| Kartezian + izotrop | 5.826e-13 | 4.56e-15 |
| Skew + izotrop | 4.405e-13 | 3.37e-15 |
| Skew + fırlanmış tam tenzor | 3.268e-13 | 2.10e-15 |
| Warped + fırlanmış tam tenzor | 7.674e-13 | 4.97e-15 |

Yəni MPFA-O xətti təzyiq sahələrini **maşın dəqiqliyində** bərpa edir —
qeyri-ortoqonal həndəsə və sıfırdan fərqli off-diaqonal tenzorla da.

## I. Ortoqonal limit: MPFA vs TPFA

4×3×1 Kartezian, `K = 200 mD · I`, təsadüfi təzyiq sahəsi:

```
max |q_MPFA − q_TPFA|                  = 9.663e-13   (nisbi 8.73e-15)
max |T_MPFA[üz, owner] − T_TPFA[üz]|   = 3.553e-15
```

Stensil səviyyəsində də: bütün daxili üzlərdə `|S_f| = 2`, əmsallar
`+T` / `−T`, və daxili üzlərin `T_bnd` sətri sıfırdır. MPFA DAXİLƏN
TPFA ÇAĞIRMIR — bu, `η=1` seçimi ilə formulyasiyanın öz limitidir
(spesifikasiya §16).

## J. Fırlanmış anizotropluq

* `test_g_rotated_tensor_flux_follows_the_rotation` — 0°/15°/37°/72°/115°
  fırlanma bucaqlarında analitik axınla uyğunluq `rtol=1e-11`.
* `test_n_rotation_invariance_of_physical_flux` — həndəsə + `K` + qradiyent
  BİRLİKDƏ 41° fırlananda axın dəyişmir (`rtol=1e-10`).
* `test_n_arbitrary_perturbation_is_distinguished_from_rotation` —
  fırlanma OLMAYAN deformasiya axını DƏYİŞİR (yəni əvvəlki test boş
  tavtologiya deyil).
* `test_o/p` — miqyas (`q ~ s`, s = 0.25/4/100) və sürüşdürmə
  (invariant) fiziki cəhətdən düzgün.

## K. Şərtlənmə

`K = R diag(500,50,10) Rᵀ`:

| grid | bölgə | `κ(C)` min | `κ(C)` maks | maks `κ(D_sub-cell)` | sinqulyar | pis-şərtlənmiş |
|---|---|---|---|---|---|---|
| Kartezian | 48 | 1.00 | 25.07 | 2.00 | 0 | 0 |
| Skew 0.35 | 48 | 1.00 | 22.64 | 2.38 | 0 | 0 |
| Skew 0.80 | 48 | 1.00 | 21.52 | 2.97 | 0 | 0 |
| Warped | 36 | 1.00 | 22.35 | 2.33 | 0 | 0 |

Sinqulyarlıq idarəsi:

* degenerativ (yastılanmış) hüceyrə → `validate_interaction_regions`
  AÇIQ `ValueError` verir (sub-üz sahəsi 0);
* sinqulyar `D_(c,v)` və ya `C` → `MPFAOSingularSystemError` +
  `diagnostics` (`region_id`, `cell`, `condition_number`, `rank`, `det`);
* `κ > condition_warning_threshold` → `ill_conditioned` bayrağı,
  hesablama davam edir.

Klipləmə / ε-əlavəsi / simmetrikləşdirmə / TPFA-ya keçid **YOXDUR**
(kodda da, testdə də təsdiqlənib).

## L. Sərhəd

| bağlanış | sərhəd DOF (3×3×2) | `T_bnd` | davranış |
|---|---|---|---|
| `DIRICHLET` (defolt) | 168 | (75, 168) | `π_bnd` XARİCİ giriş; fiziki BC İCAD EDİLMİR |
| `NEUMANN_ZERO` (opt-in) | 0 | (75, 0) | sərhəd üzlərində `max|q| = 3.427e-13` (daxili miqyas 246.9) |

Sərhəd bölgəsi `is_boundary_region` bayrağı ilə AÇIQ fərqləndirilir;
künc node-da bölgə 1 hüceyrə + 3 sub-üz + 0 naməlumdur və düzgün
işlənir.

## M. Performans və yaddaş

| grid | ncell | üz | bölgə | vaxt | µs/hüceyrə | `T_cell` nnz | nnz/sətir | yaddaş | sıx olsaydı |
|---|---|---|---|---|---|---|---|---|---|
| 4×4×2 | 32 | 128 | 75 | 0.125 s | 3901 | 1040 | 8.1 | 0.03 MB | 0.03 MB |
| 8×8×2 | 128 | 480 | 243 | 0.450 s | 3513 | 4752 | 9.9 | 0.10 MB | 0.5 MB |
| 12×12×3 | 432 | 1512 | 676 | 1.375 s | 3183 | 18360 | 12.1 | 0.32 MB | 5.2 MB |
| 16×16×4 | 1024 | 3456 | 1445 | 3.472 s | 3390 | 46368 | 13.4 | 0.74 MB | 28.3 MB |
| 24×24×4 | 2304 | 7680 | 3125 | 7.227 s | 3137 | 106400 | 13.9 | 1.65 MB | 141.6 MB |

* vaxt/hüceyrə **SABİT** (~3.1–3.9 ms) → `O(N)`, kvadratik davranış yoxdur;
* `T_cell`/`T_bnd` `scipy.sparse.csr` — sətir başına ≤18 əmsal → yaddaş `O(N)`
  (sıx saxlama `O(N²)` olardı: 2304 hüceyrədə 1.65 MB ↔ 141.6 MB);
* mütləq sürət OPTİMALLAŞDIRILMAYIB (bölgə üzrə Python dövrü) — §31
  "prematurely optimize etmə" göstərişinə uyğun. Vektorlaşdırma Phase 5D
  işidir.

## N. Testlər

```
hədəflənmiş (tests/test_mpfa_o.py):
    67 passed, 0 failed

tam dəst (pytest):
    1088 passed, 0 failed, 1 skipped
```

Yeganə skip — `tests/test_opm_import.py` ("resdata quraşdırılmayıb"),
MPFA-dan ƏVVƏL DƏ mövcud olan, opsional asılılıq skipidir.

Kateqoriya örtüyü (§28): A (API) ✔, B (sabit təzyiq) ✔, C (xətti təzyiq) ✔,
D (izotrop ortoqonal) ✔, E (skew) ✔, F (diaqonal anizotrop) ✔,
G (fırlanmış anizotrop) ✔, H (çoxnöqtəli stensil) ✔, I (lokal konservasiya) ✔,
J (sərhəd bölgəsi) ✔, K (sinqulyar sistem) ✔, L (etibarsız tenzor) ✔,
M (NaN/Inf) ✔, N (fırlanma invariantlığı) ✔, O (miqyas) ✔,
P (sürüşdürmə) ✔, Q (təsadüfi skew, 8 seed) ✔, R (TPFA reqressiyası) ✔,
§29 (anti-pseudo) ✔, §30 (istinad halı) ✔, §31 (performans/yaddaş) ✔.

## O. TPFA reqressiyası

* `imex2d/simulation/discretization.py` **HEÇ DƏYİŞMƏYİB** (git: fayl
  dəyişikliklər siyahısında yoxdur);
* `test_homogeneous_transmissibility_matches_analytic_value`,
  `test_transmissibility_uses_harmonic_mean`, fay testləri, `test_regression.py`
  (five-spot etalon: addım sayı / OOIP / RF) — hamısı dəyişmədən keçir;
* `default_flux_discretization()` HƏLƏ DƏ `TwoPointFluxDiscretization`
  qaytarır; MPFA AÇIQ seçim tələb edir;
* `IFluxDiscretization`-a əlavə olunan `supports_multipoint_stencil()`
  NON-ABSTRAKTDIR, defoltu `False` — mövcud implementasiyalar
  dəyişməmişdir;
* `test_r_mpfa_build_does_not_mutate_the_model` — MPFA qurulması modeli
  və sonradan hesablanan TPFA transmissivliyini dəyişmir.

## P. Qalan iş

**Phase 5A — TAMAMLANDI**: bölgələr, lokal sistemlər, tam tenzor,
çoxnöqtəli əmsallar, sərhəd strukturu, diaqnostika, validasiya.

**Phase 5B — TƏLƏB OLUNUR**: `ResidualAssembler`/`JacobianAssembler`
inteqrasiyası; `JacobianAssembler._build_pattern`-in 2-hüceyrəli blok
fərziyyəsinin N-hüceyrəli stensilə genişləndirilməsi; mobilite/upstream
çəkiləndirmə; sərhəd şərti qatının `π_bnd` təchizatı.

**Phase 5C — TƏLƏB OLUNUR**: mobilite/PVT analitik törəmələri
(`∂q/∂p` artıq `T_cell`-dir, amma qeyri-xətti hissə yoxdur).

**Phase 5D — TƏLƏB OLUNUR**: corner-point/struktursuz bölgə qurucusu;
fay (fault) transmissivlik çarpanlarının MPFA-da işlənməsi (hazırda
AÇIQ xəbərdarlıq verilir); bölgə dövrünün vektorlaşdırılması;
qeyri-aktiv (ACTNUM) hüceyrələr.

---

## §35 Yekun siyahı

| # | sual | cavab | sübut |
|---|---|---|---|
| 1 | Riyazi cəhətdən HƏQİQİ MPFA-O? | **BƏLİ** | `docs/mpfa_o_phase5a.md` §1–§16, `mpfa_o_local_system.py` |
| 2 | Bölgələr AÇIQ implement edilib? | **BƏLİ** | `MPFAOInteractionRegion`, `test_e/j` |
| 3 | Lokal çoxnöqtəli sistemlər? | **BƏLİ** | `C π = D p + E π_b`, `test_h` |
| 4 | Tam tenzor K lokal formulyasiyaya girir? | **BƏLİ** | `W = Γ A K D⁻¹`, `test_g` |
| 5 | Kxy/Kxz/Kyz HƏQİQƏTƏN işlənir? | **BƏLİ** | E bölməsi cədvəli |
| 6 | Üz axını >2 təzyiqdən asılı ola bilir? | **BƏLİ** | 18-hüceyrəli stensil, F bölməsi |
| 7 | Lokal konservasiya AÇIQ test edilib? | **BƏLİ** | G bölməsi, `test_i` (müstəqil yarım-axınlar) |
| 8 | Ortoqonal izotrop limit validasiya edilib? | **BƏLİ** | I bölməsi, nisbi 8.7e-15 |
| 9 | Manufactured solution validasiya edilib? | **BƏLİ** | H bölməsi |
| 10 | Fırlanmış anizotropluq validasiya edilib? | **BƏLİ** | J bölməsi |
| 11 | Skew/qeyri-ortoqonal grid validasiya edilib? | **BƏLİ** | C/E/Q testləri |
| 12 | Sinqulyar lokal sistemlər aşkarlanır? | **BƏLİ** | K bölməsi, `MPFAOSingularSystemError` |
| 13 | Etibarsız tenzor rədd edilir? | **BƏLİ** | `test_l`, `test_m` |
| 14 | Sərhəd bölgələri AÇIQ təmsil olunur? | **BƏLİ** | L bölməsi, `MPFAOBoundaryClosure` |
| 15 | TPFA TAM DƏYİŞMƏZDİR? | **BƏLİ** | O bölməsi |
| 16 | MPFA nüvəsi qeyri-xətti fizikadan müstəqildir? | **BƏLİ** | `build_mpfa_o_coefficients(grid, geometry, K, Γ)` |
| 17 | Sonlu-fərq qısayolu YOXDUR? | **BƏLİ** | pertürbasiya kodu mövcud deyil |
| 18 | MPFA daxilində gizli TPFA fallback YOXDUR? | **BƏLİ** | TPFA idxalı yoxdur; `test_anti_pseudo_*` |
| 19 | Empirik düzəliş əmsalı YOXDUR? | **BƏLİ** | yalnız `Γ` (vahid sabiti, TPFA ilə eyni) |
| 20 | Phase 5B üçün uyğundurmu? | **BƏLİ** | `T_cell` seyrək `(nface×ncell)`, `supports_multipoint_stencil()` giriş nöqtəsi |
