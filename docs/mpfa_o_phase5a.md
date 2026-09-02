# Phase 5A — TRUE MPFA-O riyazi nüvəsi (dizayn sənədi)

Bu sənəd **implementasiyadan ƏVVƏL** yazılıb (tapşırıq §3) və
`imex2d/discretization/` altındakı MPFA-O nüvəsinin DƏQİQ riyazi
formulyasiyasını təyin edir. Kod bu sənədə əməl edir; sənəd "ümumi
təsvir" deyil, **spesifikasiyadır**.

> **Miqyas xəbərdarlığı (dürüstlük)**: bu faza YALNIZ *lokal
> diskretizasiya nüvəsidir*. Residual/Jacobian/Nyuton/quyu/PVT
> inteqrasiyası BU FAZADA YOXDUR (bax §17 "Qalan iş"). TPFA
> (`TwoPointFluxDiscretization`) TOXUNULMAZ qalır.

---

## 1. Bu implementasiyada MPFA-O nə deməkdir

MPFA-O ("O-method", Aavatsmark 2002) — **lokal, çoxnöqtəli, konservativ**
sonlu-həcm diskretizasiyasıdır. Darsi qanunu

$$\mathbf u = -\,\mathbf K\,\nabla p ,\qquad
\mathbf K=\begin{bmatrix}K_{xx}&K_{xy}&K_{xz}\\K_{xy}&K_{yy}&K_{yz}\\
K_{xz}&K_{yz}&K_{zz}\end{bmatrix},\quad \mathbf K=\mathbf K^{\mathsf T}\succ 0$$

**tam tenzor** formasında saxlanılır. Metod:

1. Grid təpəsi (vertex) ətrafında **qarşılıqlı təsir bölgəsi**
   (interaction region) qurur;
2. bölgədəki hər hüceyrə payında (**sub-control volume**) təzyiqi
   **xətti** approksimasiya edir;
3. sub-üzlərdə **təzyiq kəsilməzliyini** (konstruksiyaca) və
   **axın kəsilməzliyini** (lokal xətti sistemlə) tələb edir;
4. alınan lokal sistemi həll edib **çoxnöqtəli transmissivlik
   əmsallarını** çıxarır.

MPFA-O **DEYİL** (və bu implementasiyada QADAĞANDIR):

* TPFA + düzəliş əmsalı;
* $K_n=\mathbf n^{\mathsf T}\mathbf K\mathbf n$ ilə iki-nöqtəli hesab;
* ən kiçik kvadratlar qradiyenti + iki-nöqtəli axın;
* təzyiqi pertürbasiya edib sonlu fərqlə əmsal çıxarmaq.

Bu qadağalar `tests/test_mpfa_o.py`-də **testlə** yoxlanılır (§29 anti-
pseudo-MPFA testi), sadəcə niyyət kimi qalmır.

---

## 2. Qarşılıqlı təsir bölgəsi (interaction region)

Hər **grid təpəsi** (node) $v$ üçün bir bölgə $\Omega_v$ təyin olunur.

$$\mathcal C(v)=\{\,c \;:\; v \text{ hüceyrə } c\text{-nin təpəsidir}\,\},
\qquad m_v=|\mathcal C(v)|\le 8$$

Struktur (Kartezian topologiyalı) grid üçün node indeksi
$(I,J,K)$, $0\le I\le n_x$, $0\le J\le n_y$, $0\le K\le n_z$, və

$$\mathcal C(I,J,K)=\{(i,j,k)\;:\;i\in\{I{-}1,I\}\cap[0,n_x{-}1],\;
j\in\{J{-}1,J\}\cap[0,n_y{-}1],\;k\in\{K{-}1,K\}\cap[0,n_z{-}1]\}$$

Bölgə qurulması **TAM TOPOLOJİDİR** — heç bir koordinat axtarışı /
floating-point uyğunlaşdırma YOXDUR (eyni fəlsəfə: `GeneralGridGeometry.
_build_faces`). Buna görə mürəkkəblik $O(N)$-dir (§31).

**Sub-üzlər.** $\sigma=(F,v)$ cütü — $F$ üzünün $v$ təpəsinə bitişik
dörddəbir hissəsi. $F$ üzü $(u_0,u_1,u_2,u_3)$ təpələri ilə (owner-dan
KƏNARA baxan sağ-əl sırası) verilibsə və $u_q=v$ isə:

$$\sigma=(F,v)\;\text{poliqonu}=\bigl[\;u_q,\;\tfrac{u_q+u_{q+1}}2,\;
\mathbf x_F,\;\tfrac{u_{q-1}+u_q}2\;\bigr]$$

$\mathbf x_F$ — üzün sahə-çəkili mərkəzi. Bu 4 poliqon $F$-i **tam
örtür** (müstəvi üz üçün DƏQİQ), ona görə

$$\sum_{q=0}^{3}\mathbf a_{(F,u_q)}=\mathbf a_F \equiv A_F\,\mathbf n_F$$

Hər hüceyrənin hər təpəsində **DƏQİQ 3** üzü var (hekzahedr) →
$|S(c,v)|=3$ HƏMİŞƏ, sərhəddə də. Daxili node üçün
$m_v=8$, sub-üz sayı $n_v=12$ ($8\times3/2$).

## 3. Sub-control volume (sub-cell)

$\omega_{c,v}\subset c$ — $c$ hüceyrəsinin $v$-yə ən yaxın payı
(hekzahedr üçün "səkkizdəbir"). Bu faza $\omega_{c,v}$-nin **həcmini
istifadə etmir** (akkumulyasiya Phase 5B-dədir) — yalnız onun **3 sub-
üzü** $S(c,v)$ və hüceyrə mərkəzi $\mathbf x_c$ istifadə olunur.

## 4. Bölgə daxilində təzyiq approksimasiyası

Hər sub-cell-də təzyiq **xəttidir**:

$$p^{(c,v)}(\mathbf x)=p_c+\mathbf g_{c,v}\!\cdot(\mathbf x-\mathbf x_c),
\qquad \mathbf g_{c,v}\in\mathbb R^{3}$$

$\mathbf g$ 3 naməlum daşıyır → **3 şərt** lazımdır, və sub-cell-in
məhz **3 sub-üzü** var. Hər sub-üzdə bir **kəsilməzlik nöqtəsi**
(continuity point):

$$\boxed{\;\mathbf x_\sigma=\mathbf x_v+\eta\,(\mathbf x_F-\mathbf x_v),
\qquad \eta\in(0,1]\;}$$

və şərt:

$$p^{(c,v)}(\mathbf x_\sigma)=\pi_\sigma
\qquad \forall\,\sigma\in S(c,v)$$

$\pi_\sigma$ — sub-üzün (skalyar) **kəsilməzlik təzyiqi**, bölgənin
köməkçi naməlumu. Təzyiq kəsilməzliyi bununla **konstruksiyaca**
təmin olunur: $\sigma$-nı bölüşən hər iki hüceyrə EYNİ $\pi_\sigma$-nı
işlədir.

$D_{c,v}\in\mathbb R^{3\times3}$ sətirləri $(\mathbf x_\sigma-\mathbf x_c)^{\mathsf T}$,
$\sigma\in S(c,v)$ (sabit, sənədləşdirilmiş sıra ilə):

$$D_{c,v}\,\mathbf g_{c,v}=\boldsymbol\pi_{S(c,v)}-p_c\mathbf 1
\;\;\Longrightarrow\;\;
\mathbf g_{c,v}=D_{c,v}^{-1}\bigl(\boldsymbol\pi_{S(c,v)}-p_c\mathbf 1\bigr)$$

### $\eta$ parametri

* $\eta=1$ (**DEFOLT**): $\mathbf x_\sigma=\mathbf x_F$ — üz mərkəzi.
  Bu seçim K-ortoqonal grid-də metodu **DƏQİQ** TPFA-ya reduksiya edir
  (§16), çünki Kartezian hüceyrədə $D_{c,v}$ **diaqonal** olur.
* $\eta=\tfrac12$: təpə ilə üz mərkəzi arasında orta nöqtə. Xətti
  sahələr üçün YENƏ DƏQİQDİR, amma Kartezian grid-də stensil
  **2-nöqtəli DEYİL** (4-nöqtəli) — bu, metodun HƏQİQİ,
  sənədləşdirilmiş xassəsidir, qüsur deyil; `test_mpfa_o.py`-də hər
  ikisi yoxlanılır.

## 5. Tam tenzor $\mathbf K$ axına necə daxil olur

Sub-cell axın sıxlığı:

$$\mathbf u^{(c,v)}=-\,\mathbf K_c\,\mathbf g_{c,v}
=-\,\mathbf K_c D_{c,v}^{-1}\bigl(\boldsymbol\pi_{S}-p_c\mathbf 1\bigr)$$

$\mathbf K_c$ **TAM 3×3 SİMMETRİK MATRİS** kimi daxil olur — heç bir
mərhələdə $\mathbf n^{\mathsf T}\mathbf K\mathbf n$ skalyarına
"yığılmır" (bu QADAĞANDIR, §4). $K_{xy},K_{xz},K_{yz}$ komponentləri
$\mathbf K_c D^{-1}$ hasilində birbaşa iştirak edir.

## 6. Yarım-axınlar (half-fluxes) və işarə konvensiyası

$\sigma=(F,v)$ sub-üzü üçün **istiqamətləndirilmiş sahə vektoru**
$\mathbf a_\sigma$ — HƏMİŞƏ $F$-in **owner**-indən KƏNARA. Hüceyrə
$c$-yə görə:

$$\mathbf a_\sigma^{(c)}=
\begin{cases}+\mathbf a_\sigma,& c=\mathrm{owner}(F)\\
-\mathbf a_\sigma,& c=\mathrm{neighbor}(F)\end{cases}$$

$c$-dən $\sigma$ vasitəsilə **çıxan** axın:

$$q_{c,\sigma}
=-\,\Gamma\;\mathbf a_\sigma^{(c)\mathsf T}\,\mathbf K_c\,\mathbf g_{c,v}
=-\,\Gamma\;\mathbf a_\sigma^{(c)\mathsf T}\mathbf K_c D_{c,v}^{-1}
\bigl(\boldsymbol\pi_{S(c,v)}-p_c\mathbf 1\bigr)$$

$\Gamma$ — vahid sisteminin Darsi sabiti (`UnitSystem.darcy_constant`,
TPFA ilə **EYNİ**). Çəki sətri:

$$\mathbf w_{c,\sigma}=\Gamma\,\mathbf a_\sigma^{(c)\mathsf T}
\mathbf K_c D_{c,v}^{-1}\in\mathbb R^{1\times3}
\;\Longrightarrow\;
q_{c,\sigma}=-\!\!\sum_{\tau\in S(c,v)}\!\! w_{c,\sigma,\tau}\,\pi_\tau
\;+\;\Bigl(\sum_\tau w_{c,\sigma,\tau}\Bigr)p_c$$

**İşarə konvensiyası (bütün kodda vahid)**: `flux[f] > 0` ⟺ axın
`owner(f)`-dan `neighbor(f)`-ə (yəni `face_normals[f]` istiqamətində).

## 7. Lokal sistem

Naməlumlar (DAXİLİ sub-üzlərin $\pi$-ləri, və Neumann bağlanışında
sərhəd sub-üzləri də):

$$\mathbf x_{\text{local}}=\boldsymbol\pi_{\text{unk}}\in\mathbb R^{n_u}$$

Tənliklər — **hər naməlum $\pi$ üçün BİR sətir**:

| sətir tipi | tənlik | mənası |
|---|---|---|
| daxili sub-üz $\sigma$ | $q_{o,\sigma}+q_{n,\sigma}=0$ | axın kəsilməzliyi |
| sərhəd, `NEUMANN_ZERO` | $q_{o,\sigma}=0$ | axınsız bağlanış |
| sərhəd, `DIRICHLET` | — (naməlum DEYİL) | $\pi_\sigma$ XARİCDƏN verilir |

Matris formasında:

$$\boxed{\;C\,\boldsymbol\pi_{\text{unk}}
= D\,\mathbf p_{\mathcal C(v)} + E\,\boldsymbol\pi_{\text{bnd}}\;}$$

* **sətir** $r$ ↔ naməlum sub-üz $\sigma_r$-in bağlanış tənliyi;
* **sütun** $t$ ↔ naməlum $\pi_{\sigma_t}$;
* $D$ sütunları ↔ bölgədəki hüceyrə təzyiqləri $p_{c_0..c_{m-1}}$;
* $E$ sütunları ↔ XARİCDƏN verilən sərhəd $\pi$-ləri (Dirichlet).

Hər əmsala **hansı həndəsə** ($\mathbf a_\sigma$, $D_{c,v}$) və
**hansı tenzor** ($\mathbf K_c$) daxil olduğu §6 düsturundan birbaşa
oxunur; kod `MPFAOLocalSystem`-də bu matrisləri **açıq atribut** kimi
saxlayır (§8 tələbi — "opaque helper" QADAĞANDIR).

## 8. Axın bərpası və çoxnöqtəli əmsallar

Sub-üz axını (owner-dan çıxan, §6 işarəsi):

$$q_\sigma=F\boldsymbol\pi_{\text{unk}}+G\,\mathbf p+H\boldsymbol\pi_{\text{bnd}}$$

$\boldsymbol\pi_{\text{unk}}$-i əvəz edərək:

$$\boxed{\;\mathbf q_{\text{sub}}
=\underbrace{\bigl(FC^{-1}D+G\bigr)}_{T^{\text{cell}}}\mathbf p
+\underbrace{\bigl(FC^{-1}E+H\bigr)}_{T^{\text{bnd}}}\boldsymbol\pi_{\text{bnd}}\;}$$

$T^{\text{cell}}\in\mathbb R^{n_\sigma\times m_v}$ — **bölgə
transmissivlik matrisi**. Tam üz axını sub-üz paylarının cəmidir:

$$q_F=\sum_{v\in\mathcal V(F)}q_{(F,v)}
=\sum_{c\in\mathcal S_F}T_{F,c}\,p_c ,\qquad
|\mathcal V(F)|=4 \text{ (hekzahedr üzü)}$$

Qlobal `T` matrisi `(nface × ncell)` — hər üzün 4 küncündən
gələn payların **cəmi**.

## 9. Lokal konservasiya

$q_\sigma$ **BİR DƏFƏ** (owner tərəfindən) təyin olunur; neighbor
axını riyazi olaraq $-q_\sigma$-dır → sub-üz səviyyəsində konservasiya
**identikdir**. Bundan asılı OLMAYARAQ, lokal sistemin **həqiqətən
həll olunduğunu** yoxlamaq üçün müstəqil qalıq hesablanır:

$$\varepsilon_\sigma=\bigl|q_{o,\sigma}+q_{n,\sigma}\bigr|,\qquad
\varepsilon_{\max}=\max_\sigma \varepsilon_\sigma$$

burada $q_{o,\sigma}$ və $q_{n,\sigma}$ **ayrı-ayrılıqda**, hər
tərəfin ÖZ $\mathbf g_{c,v}$-sindən hesablanır (§15 tələbi). Üz
səviyyəsində: $q_{i,F}+q_{j,F}=0$ (§9).

## 10. Sərhəd bölgələri

`MPFAOBoundaryClosure` **AÇIQ** enum-dur:

* `DIRICHLET` (defolt) — sərhəd $\pi_\sigma$ **NAMƏLUM DEYİL**, xarici
  giriş kimi qəbul edilir; nüvə $T^{\text{bnd}}$ əmsallarını qaytarır.
  Bu, "fiziki sərhəd şərti icad etmək" DEYİL (§20) — sadəcə gələcək
  BC qatının dəyər ötürməsi üçün **riyazi struktur**.
* `NEUMANN_ZERO` — $q_{o,\sigma}=0$ (axınsız). Bu, HƏQİQİ fiziki
  şərtdir, ona görə **opt-in**, defolt DEYİL, və ayrıca test edilir.

Bölgə `is_boundary_region` bayrağı daşıyır. Nüvə heç vaxt "hər
bölgənin 8 hüceyrəsi var" fərz ETMİR.

## 11. Degenerativ / sinqulyar hallar

| yoxlama | reaksiya |
|---|---|
| $\mathbf K$-də NaN/Inf | `MPFAOTensorError` — **rədd** |
| $\lambda_{\min}(\mathbf K)\le0$ (qeyri-SPD) | `MPFAOTensorError` — **rədd** |
| $D_{c,v}$ sinqulyar (degenerativ sub-cell) | `MPFAOSingularSystemError` |
| $C$ sinqulyar / rank azlığı | `MPFAOSingularSystemError` |
| $\kappa(C)>$ `condition_warning_threshold` | diaqnostikada **xəbərdarlıq**, hesablama davam edir |

**QADAĞAN** (§18/§19): eigenvalue klipləmə, diaqonala $\varepsilon$
əlavəsi, simmetrikləşdirmə, səssiz TPFA-ya keçid, səssiz
requlyarizasiya. Heç biri kodda YOXDUR.

Diaqnostika (`MPFAOLocalSystem.diagnostics()`): `region_id`,
`condition_number`, `rank`, `determinant`, `singular`,
`n_unknowns`, `n_cells`, `is_boundary_region`.

## 12. İndeksləmə konvensiyaları

* **hüceyrə**: qlobal `cell` (0..ncell-1); bölgə daxilində `local`
  indeks `region.cell_local[cell]`.
* **üz**: `GeneralGridGeometry` qlobal `face_index`; owner/neighbor
  `face_owner`/`face_neighbor` (`-1` = sərhəd).
* **node**: struktur grid-də `node_id = (K*(ny+1) + J)*(nx+1) + I`.
* **sub-üz**: `(face_index, node_id)` cütü — qlobal unikal.
* **sərhəd DOF**: sərhəd sub-üzləri qlobal sırada nömrələnir
  (`boundary_dof[(face, node)]`).
* $S(c,v)$ sırası: `HEX_FACE_VERTEX_INDICES` açar sırası
  (`Z-,Z+,Y-,Y+,X-,X+`) ilə **deterministik**.

## 13. 2D / 3D

Nüvə **TAM 3D**-dir — tək kod yolu, $\mathbb R^3$-də $3\times3$
$D_{c,v}$ və tam $3\times3$ $\mathbf K$. **AYRI 2D alqoritmi YOXDUR**
(§26 qadağası: "XY/XZ/YZ müstəvilərinə 2D alqoritm tətbiq etmək"
saxta 3D-dir).

"2D" məsələ ($n_z=1$) bu nüvədə **tək layı olan 3D grid** kimi
həll olunur: üst/alt üzlər sərhəddir və §10-a görə açıq bağlanış
alır. `NEUMANN_ZERO` bağlanışı ilə bu, fiziki cəhətdən **məhdud
(confined) lay** deməkdir və $K_{xz}=K_{yz}=0$ olduqda klassik 2D
MPFA-O ilə üst-üstə düşür.

## 14. Gözlənilən stensil ölçüsü

| konfiqurasiya | bir sub-üz | tam üz |
|---|---|---|
| daxili node, 3D | $\le 8$ hüceyrə | $\le 18$ hüceyrə |
| daxili node, tək lay ($n_z=1$) | $\le 4$ | $\le 6$ |
| Kartezian + K-ortoqonal $\mathbf K$ + $\eta=1$ | **2** | **2** |

Sonuncu sətir §16-nın nəticəsidir və **testlə** təsdiqlənir.

## 15. Geometriyadan asılı əvvəlcədən hesablama

$T^{\text{cell}},T^{\text{bnd}}$ YALNIZ `{həndəsə, topologiya,
K, Γ}`-dan asılıdır — **təzyiq/doyma/mobilite/PVT-dən ASILI DEYİL**.
Ona görə `MPFAOCoefficients` bir dəfə qurulur və keşlənir; Phase 5B-də
mobilite vurulması **bölgələri yenidən qurmadan** mümkün olacaq (§24).

## 16. TPFA-ya reduksiya şərtləri

$\eta=1$, ortoqonal Kartezian hüceyrə, K-ortoqonal tenzor
($\mathbf K\mathbf n_F \parallel \mathbf n_F$ hər üz üçün — məs.
$\mathbf K=k\mathbf I$ və ya oxa-uyğun diaqonal) olduqda:

$D_{c,v}=\operatorname{diag}(\pm h_x/2,\pm h_y/2,\pm h_z/2)$ →
$\mathbf a_\sigma^{\mathsf T}\mathbf K_c D_{c,v}^{-1}$ YALNIZ $\sigma$-nın
öz sütununda sıfırdan fərqlidir → hər kəsilməzlik tənliyi TƏK $\pi_\sigma$
saxlayır → $\pi_\sigma=(k_o p_o/h_o+k_n p_n/h_n)/(k_o/h_o+k_n/h_n)$ və

$$q_\sigma=\Gamma\frac{A_\sigma}{h_o/k_o+h_n/k_n}\,(p_o-p_n)$$

— yəni **eyni harmonik-orta TPFA düsturu** (`TwoPointFluxDiscretization.
build`), sub-üz sahələri üzrə cəmləndikdə $A_F$ verir. Bu, MPFA-nın
DAXİLƏN TPFA çağırması ilə DEYİL, formulyasiyanın öz limitidir (§13).

## 17. Qalan iş (bu fazada YOXDUR)

* **5B**: residual/Jacobian inteqrasiyası, mobilite/upstream çəkiləndirmə,
  `JacobianAssembler._build_pattern` çoxnöqtəli stensilə genişləndirilməsi.
* **5C**: analitik törəmələr (mobilite/PVT törəmələri).
* **5D**: struktursuz/corner-point bölgə qurucusu, fayların (fault)
  MPFA-da işlənməsi, paralelləşdirmə.
