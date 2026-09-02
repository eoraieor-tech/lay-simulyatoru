# Phase 5B-1 — MPFA-O qlobal residual və seyrək stensil inteqrasiyası

Bu sənəd **implementasiyadan ƏVVƏL** yazılıb (tapşırıq §2). Phase 5A
nüvəsi (bax [`mpfa_o_phase5a.md`](mpfa_o_phase5a.md)) DƏYİŞDİRİLMİR —
onun ətrafına inteqrasiya qatı qurulur.

> **Miqyas xəbərdarlığı**: bu faza QALIQ (residual) səviyyəsindədir.
> **Qeyri-xətti Jacobian və Nyuton inteqrasiyası BU FAZADA YOXDUR**
> (§2.5). MPFA rejimində `JacobianAssembler`/`FullyImplicitEngine`/
> `ImpesEngine` AÇIQ xəta verir — saxta uyğunluq YARADILMIR (§24).

---

## 1. Mövcud arxitektura: TPFA qalığa NECƏ çatır

```
FullyImplicitEngine.__init__
    flux_discretization.build(model) → DiscretizedGrid
        .connections            (nconn: cell_a, cell_b, axis)
        .transmissibility       (nconn,)
        .pore_volume/.cell_volume
        .compute_flux(d_phi) = transmissibility · d_phi
    ↓
ResidualAssembler(model, grid, wells, relperm, pvt, capillary)
    .fluid_state(state)   → μ, B, λ = kr/μ, Pc      (HÜCEYRƏ üzrə)
    .potentials(state, f) → ΔΦ_w, ΔΦ_o             (ÜZ üzrə fərq)
    .face_fluxes()        → grid.compute_flux(ΔΦ) · (λ/B)_upstream
    .net_influx()         → np.add.at ilə cell_a(−)/cell_b(+)
    .well_rates()         → Peaceman (MPFA-dan KƏNAR)
    .accumulation()       → PV(p)·S/B
    .residual()           → R = (acc_new − acc_old)/Δt − influx − q_well
    ↓
JacobianAssembler  →  NewtonSolver  →  AdaptiveTimeStepper
```

Mövcud potensial konvensiyası (`residual.py::potentials`,
`impes_engine.py::_face_potential_terms` — İKİSİ DƏ EYNİ):

```
ΔΦ_o(a→b) = (p_a − p_b) − ρ̄_o · g · (D_a − D_b) · PA_TO_BAR
ΔΦ_w(a→b) = (p_a − p_b) − (Pc_a − Pc_b) − ρ̄_w · g · (D_a − D_b) · PA_TO_BAR
ρ̄_α      = ½ (ρ_α,a + ρ_α,b),      ρ_α = ρ_α,səth / B_α
```

**İşarə konvensiyası (mövcud)**: `flux > 0` ⟺ axın `cell_a → cell_b`;
`net_influx` `cell_a`-ya `−flux`, `cell_b`-yə `+flux` yazır;
`R = akkumulyasiya/Δt − influx − quyu`, `R = 0` tarazlıqdır.

**Kritik məhdudiyyətlər (auditlə təsdiqləndi)**:

| yer | fərziyyə |
|---|---|
| `ResidualAssembler.__init__` | `grid.transmissibility` MÖVCUDDUR |
| `ResidualAssembler.face_fluxes` | axın = `f(ÜZ-başına ΔΦ)` |
| `JacobianAssembler._build_pattern` | hər üz üçün DƏQİQ 2 hüceyrəli blok |
| `JacobianAssembler._flux` | `∂ΔΦ/∂p_a=+1`, `∂ΔΦ/∂p_b=−1` |
| `ImpesEngine.__init__` | `_discretization.transmissibility` MÖVCUDDUR |

Birinci ikisi bu fazada ümumiləşdirilir; son üçü MPFA rejimində AÇIQ
RƏDD edilir (Phase 5B-2).

---

## 2. Yeni MPFA yolu

```
ReservoirState (p, Sw)
    ↓  fluid_state()  (MÜŞTƏRƏK — dəyişmir)
    ↓
HÜCEYRƏ üzrə faza potensialı  Φ_α,c        ← §4 (YENİ)
    ↓
MPFA HƏNDƏSİ operatoru  T_conn  (nconn × ncell, csr)   ← §3, DÖVLƏTDƏN ASILI DEYİL
    ↓
baza Darsi axını   q_pot,f = Σ_j T_conn[f,j] Φ_α,j
    ↓
upstream çəkiləndirmə   sign(q_pot,f)                   ← §6
    ↓
faza axını   q_α,f = q_pot,f · (λ_α/B_α)_upstream
    ↓
net_influx (DƏYİŞMƏYƏN `np.add.at` məntiqi)
    ↓
R_i = (acc_new − acc_old)/Δt − influx_i − q_well,i      (DƏYİŞMİR)
```

**Müştərək qalan** (§14): akkumulyasiya, PVT, nisbi keçiricilik,
kapilyar təzyiq, quyu üzvləri, material balansı, konvergensiya ölçüsü,
zaman addımı. **Yalnız məkan diskretizasiyası** fərqlənir.

---

## 3. Qlobal stensil: lokaldan qlobala

Phase 5A `MPFAOCoefficients.T_cell` `(nface × ncell)` csr verir —
qarşılıqlı təsir bölgələrindəki `T_cell^(v)` bloklarının fiziki üzlər
üzrə CƏMİ (bir üzün 4 küncündən 4 pay; COO təkrarlarının cəmlənməsi ilə).
**Bu əmsallar OLDUĞU KİMİ saxlanılır** (§3 — heç bir yenidən hesablama,
heç bir yığma).

Qlobal inteqrasiya YALNIZ SƏTİR SEÇİMİDİR:

```
T_conn = T_face[connection_faces, :]        (nconn × ncell) csr
```

`connection_faces` — `GeneralGridGeometry.connection_faces()`, hər
`Connections` girişi üçün QLOBAL üz indeksi (Phase 5A-da qurulub).
Bu, DETERMİNİSTİK BİYEKSİYADIR:

```
daxili fiziki üz  ⟷  DƏQİQ BİR Connections girişi
sərhəd fiziki üz  →  Connections-da YOXDUR (axını NEUMANN_ZERO ilə ≡ 0)
```

Üz axını:

$$q_f=\sum_{j\in S_f} T_{f,j}\,\Phi_j,\qquad |S_f|>2$$

Hüceyrə-hüceyrə (divergensiya) operatoru — sabit mobilitə halında
xətti operator və HƏQİQİ bağlantı naxışı:

$$A = D\,T_{conn},\qquad
D\in\mathbb R^{ncell\times nconn},\;
D_{a,f}=-1,\; D_{b,f}=+1$$

`A` sətri `i` — `R_i`-nin HANSI hüceyrə təzyiqlərindən asılı olduğunu
göstərir; bu, `Connections` qonşuluğundan GENİŞDİR (diaqonal
bağlantılar) və `global_stencil_pattern()` ilə açılır (§22).

### İkiqat sayılmanın qarşısı (§8)

| mərhələ | sahiblik konvensiyası |
|---|---|
| bölgə sub-üzü `σ=(F,v)` | `q_σ` YALNIZ owner tərəfindən təyin olunur |
| sub-üz → fiziki üz | 4 künc payı `T_face[F,:]`-də CƏMLƏNİR (Phase 5A) |
| fiziki üz → `Connections` | biyeksiya (`connection_faces`), TƏK giriş |
| `Connections` → hüceyrə | `np.add.at`: `cell_a −q`, `cell_b +q` (DƏYİŞMİR) |

Hər fiziki üz qlobal balansa **DƏQİQ BİR DƏFƏ** daxil olur. Sərhəd
üzləri `Connections`-da yoxdur və `NEUMANN_ZERO` bağlanışı ilə axını
sıfırdır — nə buraxılır, nə də iki dəfə sayılır. Testlər:
`test_p_no_double_counting_*`.

---

## 4. Təzyiq / potensial (§12)

**BU FAZA: TAM POTENSİAL** (`pressure-only` DEYİL) — mövcud qalıq
cazibə və kapilyar üzvlərini TƏLƏB EDİR, onları susdurmaq fiziki
cəhətdən yanlış olardı.

MPFA operatoru SKALYAR sahəyə tətbiq olunur, ona görə HÜCEYRƏ üzrə
potensial lazımdır (üz-başına fərq yox):

```
Φ_o,c = p_c          − ρ_o,c · g · D_c · PA_TO_BAR
Φ_w,c = p_c − Pc_c   − ρ_w,c · g · D_c · PA_TO_BAR
ρ_α,c = ρ_α,səth / B_α,c        (EYNİ düstur, `residual.py`-dən)
```

`g`, `PA_TO_BAR`, `ρ_səth`, `B`, `Pc` — HAMISI mövcud koddan;
**İKİNCİ cazibə konvensiyası İCAD EDİLMİR** (§12).

### Bilinən və SƏNƏDLƏŞDİRİLMİŞ fərq

TPFA üz-başına ORTALANMIŞ sıxlıq (`ρ̄ = ½(ρ_a+ρ_b)`) işlədir. Çoxnöqtəli
stensildə "həmin iki hüceyrə" anlayışı YOXDUR, ona görə HÜCEYRƏ sıxlığı
işlədilir. İki-nöqtəli limitdə:

```
MPFA : T[(p_a−p_b) − (Pc_a−Pc_b) − g·PA_TO_BAR·(ρ_a D_a − ρ_b D_b)]
TPFA : T[(p_a−p_b) − (Pc_a−Pc_b) − g·PA_TO_BAR· ρ̄ ·(D_a − D_b)]
```

Bunlar **EYNİDİR** o zaman ki:
* cazibə yoxdur (`D_a = D_b`, düz lay) — **VƏ YA**
* faza sıxlığı hüceyrələr arasında bərabərdir (`ρ_a = ρ_b`).

Fərq `O(Δρ · ΔD)` tərtibindədir. Kapilyar üzv **DƏQİQ EYNİDİR**
(`Pc` onsuz da hüceyrə kəmiyyətidir). Bu, qüsur deyil — çoxnöqtəli
sxem üçün standart və yeganə ardıcıl seçimdir; testlər hər üç halı
(bərabər, fərqli sıxlıq, cazibəsiz) AÇIQ ölçür
(`test_h_*`, `test_o_gravity_*`).

---

## 5. Qalıq yığımı (§14)

`ResidualAssembler` DAXİLİNDƏ DAR bir budaqlanma:

```python
self._multipoint = grid.supports_multipoint_stencil()   # duck-typed, defolt False

def face_fluxes(self, state, fluid):
    if self._multipoint:
        return self._multipoint_face_fluxes(state, fluid)
    ...  # MÖVCUD TPFA yolu — BİR SƏTİR BELƏ DƏYİŞMİR
```

`net_influx`, `well_rates`, `accumulation`, `residual`,
`scaled_residual_norm`, `material_balance_error` — **TOXUNULMUR**.

Qalıq (dəyişmir):

$$R_{\alpha,i}=\frac{(PV\,S_\alpha/B_\alpha)^{n+1}_i-(PV\,S_\alpha/B_\alpha)^{n}_i}{\Delta t}
-\underbrace{\sum_{f\ni i}\pm q_{\alpha,f}}_{\text{influx}}-q^{well}_{\alpha,i}$$

---

## 6. Mobilitə və upstream (§9/§10/§11)

**Ayrılıq (§11)**: `T_conn` YALNIZ `{həndəsə, topologiya, K, Γ}`-dan
asılıdır → BİR DƏFƏ qurulur (`MPFAODiscretization.build`) və hər
qalıq qiymətləndirməsində TƏKRAR İSTİFADƏ olunur. Lokal sistemlər
Nyuton iterasiyasında YENİDƏN QURULMUR.

**Mobilitə**: mövcud `IRelativePermeabilityProvider`/`IPVTProvider`
işlədilir — İKİNCİ relperm çərçivəsi YARADILMIR. Nəqledici mobilitə
`M_α = λ_α/B_α = kr_α/(μ_α B_α)` — `face_fluxes`-dəki EYNİ ifadə.

**Upstream (§10)**: istiqamət `p_owner − p_neighbor`-dan DEYİL, HƏQİQİ
çoxnöqtəli axından çıxarılır:

```
q_pot,f = Σ_j T_conn[f,j] Φ_α,j
upstream(f) = cell_a  if q_pot,f ≥ 0  else  cell_b
q_α,f = q_pot,f · M_α[upstream(f)]
```

İki-nöqtəli limitdə `q_pot = T·ΔΦ`, `T > 0` → `sign(q_pot) = sign(ΔΦ)`
→ TPFA-nın `where(ΔΦ ≥ 0, a, b)` qaydası ilə **BİRƏBİR** eyni
(`≥` daxil olmaqla). Test: `test_o_upwind_matches_tpfa_in_two_point_limit`.

---

## 7. Sərhəd (§13)

Mövcud simulyator QAPALI rezervuar həll edir: `Connections` YALNIZ
daxili üzləri saxlayır, sərhəd üzü ANLAYIŞI residual qatında YOXDUR
→ implisit **axınsız (no-flow)** sərhəd.

Ona görə qalıq yolu Phase 5A-nın `MPFAOBoundaryClosure.NEUMANN_ZERO`
bağlanışını işlədir — bu, "fiziki BC icad etmək" DEYİL, mövcud
simulyatorun ARTIQ tətbiq etdiyi şərtin MPFA qarşılığıdır.

| BC növü | 5B-1 statusu |
|---|---|
| no-flow (Neumann-0) | **DƏSTƏKLƏNİR** — qalıq yolunun defoltu |
| Dirichlet təzyiq | **struktur HAZIRDIR** (`T_bnd`, `boundary_dofs`), residual qatı hələ dəyər ötürmür — Phase 5B-2 |
| ümumi Neumann (sıfırdan fərqli axın) | **YOXDUR** — AÇIQ bildirilir |

---

## 8. Dəstəklənməyən xüsusiyyətlərin RƏDDİ (§26/§27)

`MPFADiscretizedGrid.unsupported_features` siyahısı daşıyır. Qalıq
qatına qoşulanda siyahı BOŞ DEYİLSƏ `NotImplementedError` verilir —
səssiz, fiziki cəhətdən YANILDICI nəticə QADAĞANDIR.

* **Fay (fault)**: Phase 5A MPFA-da fay çarpanlarını tətbiq ETMİR.
  TPFA çarpanlarını MPFA-ya "yapışdırmaq" QADAĞANDIR (§26/§32) →
  faylı modeldə MPFA qalığı **RƏDD EDİLİR**.
* **NaN/Inf təzyiq**: hansı hüceyrələr olduğu göstərilərək rədd edilir.
* **Ölçü uyğunsuzluğu**: operator/vektor ölçüləri yoxlanılır.
* **Etibarsız tenzor/həndəsə**: Phase 5A-dakı yoxlamalar
  (`MPFAOTensorError`, `MPFAOSingularSystemError`) DƏYİŞMƏDƏN işləyir.

---

## 9. Seyrəklik və performans (§5/§6/§28)

* `T_face`, `T_conn`, `A` — HAMISI `scipy.sparse.csr_matrix`.
  `np.zeros((ncell, ncell))` HEÇ BİR YERDƏ YOXDUR.
* Yığım COO üçlükləri ilə (Phase 5A-da) → `tocsr()`; CSR girişləri
  dövr daxilində DƏYİŞDİRİLMİR.
* `T_conn` sətir seçimidir (`T_face[rows]`) — `O(nnz)`.
* `A = D @ T_conn` — seyrək-seyrək hasili, sətir başına məhdud
  (`≤ 18 × 2`) girişlə → `O(N)`.
* Bölgələr/lokal əmsallar/qlobal naxış **KEŞLƏNİR** (dövlətdən asılı
  deyil); qalıq qiymətləndirməsi yalnız `T_conn @ Φ` + mobilitə
  vurmasıdır.

---

## 10. Sinaq matrisi

`tests/test_mpfa_o_global_assembly.py` — A–T (§29) + §16/§17/§20/§30/§31.

---

## 11. Phase 5B-2 / 5C üçün qalan (§2.5/§23)

**Bu fazada Jacobian YENİDƏN QURULMUR.** Sabit mobilitəli halda
`A = D T_conn` XƏTTİ MPFA operatorudur, AMMA:

```
xətti MPFA operatoru   ≠   tam qeyri-xətti rezervuar Jacobian-ı
```

Çatışmayan: `∂M_α/∂p`, `∂M_α/∂S_w`, `∂Φ/∂S_w` (kapilyar), quyu
törəmələri çoxnöqtəli stensillə, `JacobianAssembler._build_pattern`-in
N-hüceyrəli ümumiləşdirilməsi, CPR ön-şərtçisinin yeni naxışa
uyğunlaşdırılması. Sonlu-fərq qısayolu **QADAĞANDIR** (§23).

MPFA rejimində `JacobianAssembler`, `FullyImplicitEngine`, `ImpesEngine`
AÇIQ `NotImplementedError` verir və bu sənədə istinad edir.
