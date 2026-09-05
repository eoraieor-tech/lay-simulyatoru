# IMEX-2D — Arxitektura sənədi

Refaktorinq sənədi. Yeni mühəndislik hesablaması əlavə edilməyib;
mövcud fizika olduğu kimi saxlanılıb və qatlara bölünüb.

## 0. Refaktorinqin doğruluğunun sübutu

Eyni giriş məlumatı ilə köhnə və yeni kod eyni nəticəni verir:

| Göstərici | Köhnə (`core.py`) | Yeni (qatlı arxitektura) |
|---|---|---|
| Recovery Factor | 16.840 % | 16.840 % |
| Zaman addımı sayı | 4314 | 4314 |
| OOIP | 1 029 064.3 m³ | 1 029 064.3 m³ |

Bu, refaktorinqin **davranışı qorunan** (behaviour-preserving) olduğunu göstərir.

---

## 1. Köhnə arxitekturanın problemləri

### 1.1 Model anlayışı ümumiyyətlə yox idi
`Simulator2D.__init__` altı ayrı konfiqurasiya obyekti alıb massivləri
**özü qururdu**. Model heç yerdə obyekt kimi mövcud deyildi — o, yalnız
simulyatorun daxili vəziyyəti idi. Nəticə: modeli yadda saxlamaq,
müqayisə etmək, iki simulyasiyada təkrar işlətmək mümkün deyildi.

### 1.2 Geologiya ilə simulyasiya bir-birinə qarışmışdı
Keçiricilik sahəsi `MainWindow.perm_field()` metodunda yaranırdı — yəni
geoloji modelləşdirmə **interfeys kodunun içində** idi. Streamlit-dəki
3D lay modelini bu proqrama bağlamaq mümkün deyildi.

### 1.3 UI biznes məntiqini daşıyırdı
Model yoxlaması `QMessageBox` çağırışları ilə `run_simulation()` içində
idi. Skriptdən istifadə edəndə heç bir yoxlama işə düşmürdü.

### 1.4 Rendering hesablama ilə eyni sinifdə idi
`update_map`, `update_results` metodları həm məlumat seçirdi, həm çəkirdi.
Başsız (headless) hesabat və ya avtomatik test mümkün deyildi.

### 1.5 Sabitlər hər yerdə səpələnmişdi
`rtol=1e-8`, `drop_tol=1e-4`, 25 addımlıq ön-şərtçi yeniləməsi, 1.15
artım əmsalı — mühərrikin içində hardcode edilmişdi.

### 1.6 Genişlənmə nöqtəsi yox idi
PVT əlavə etmək üçün mühərrikin özünü açıb dəyişmək lazım gələcəkdi
(Open/Closed prinsipinin pozulması).

---

## 2. Yeni qat strukturu

```
imex2d/
  domain/         Model obyektləri. Heç bir daxili qatdan asılı deyil.
  interfaces/     Provider müqavilələri (ABC). Yalnız domain-dən asılıdır.
  application/    Layihə, iş axını, servislər, konfiqurasiya.
  simulation/     Hesablama mühərriki və köməkçiləri.
  rendering/      Çəkmə. Qt tanımır.
  ui/             Qt interfeysi. Yalnız orkestrasiya.
app.py            Composition root — asılılıqlar burada bağlanır.
```

Asılılıq istiqaməti həmişə içəriyə doğrudur:

```
ui ──> rendering ──┐
 │                 ├──> application ──> simulation ──> interfaces ──> domain
 └─────────────────┘
```

`domain` heç kimdən asılı deyil. `simulation` konkret provider siniflərini
deyil, `interfaces`-dəki abstraksiyaları tanıyır.

---

## 3. UML — sinif diaqramı (mətn formatı)

```
┌─ domain ─────────────────────────────────────────────────────────┐
│                                                                   │
│  CartesianGrid           CellGeometry            PropertyMap      │
│  + nx, ny, nz            + dx, dy, dz            + name           │
│  + ncell                 + volumes()             + values         │
│  + index(i,j,k)          + face_areas()          + unit           │
│  + build_connections()   + face_half_distances() + uniform()      │
│           │                      │                     │         │
│           └──────────┬───────────┘                     │         │
│                      ▼                                 ▼         │
│              GeologicalModel                    RockProperties    │
│              + grid, geometry                   + porosity        │
│              + property_maps                    + permx, permy    │
│              + regions                          + compressibility │
│              + horizons: Horizon[]                                │
│              + faults: Fault[]                  FluidProperties   │
│              + validate()                       (PLACEHOLDER)     │
│                      │                                            │
│                      │ ReservoirModelBuilder.build()              │
│                      ▼                                            │
│  ╔══════════════════════════════════════════════════════════════╗ │
│  ║  ReservoirModel   « yeganə həqiqət mənbəyi »                 ║ │
│  ╠══════════════════════════════════════════════════════════════╣ │
│  ║  + grid: CartesianGrid                                       ║ │
│  ║  + geometry: CellGeometry                                    ║ │
│  ║  + rock: RockProperties                                      ║ │
│  ║  + fluids: FluidProperties           (placeholder)           ║ │
│  ║  + property_maps: Dict[str, PropertyMap]                     ║ │
│  ║  + regions: RegionSet                                        ║ │
│  ║  + fault_references: FaultReference[]                        ║ │
│  ║  + horizon_references: HorizonReference[]                    ║ │
│  ║  + wells: Well[]                                             ║ │
│  ║  + initial_conditions: InitialConditions  (placeholder)      ║ │
│  ║  + scal_parameters: CoreyParameters                          ║ │
│  ║  + units: UnitSystem                                         ║ │
│  ╟──────────────────────────────────────────────────────────────╢ │
│  ║  + connections()      + pore_volume()                        ║ │
│  ║  + active_wells()     + validate()      + summary()          ║ │
│  ╚══════════════════════════════════════════════════════════════╝ │
└───────────────────────────────────────────────────────────────────┘

┌─ interfaces (yalnız müqavilə, implementasiya YOXDUR) ────────────┐
│  «interface» IPVTProvider                                        │
│      + oil_fvf(p)          + oil_viscosity(p)                    │
│      + water_fvf(p)        + water_viscosity(p)                  │
│      + total_compressibility(p, sw)                              │
│                                                                   │
│  «interface» IRelativePermeabilityProvider                       │
│      + krw(sw, region)     + kro(sw, region)                     │
│      + saturation_limits() + endpoint_water_mobility()           │
│      + max_fractional_flow_derivative()                          │
│                                                                   │
│  «interface» ICapillaryPressureProvider                          │
│      + pcow(sw, region)    + dpcow_dsw(sw, region)               │
│                                                                   │
│  «interface» IInitializationProvider                             │
│      + initialize(model) -> InitialState                         │
│                                                                   │
│  «interface» ILinearSolver     «interface» IProgressReporter     │
│      + solve(A, b, x0)             + report(fraction, msg)       │
│                                                                   │
│  «interface» ISimulationEngine                                   │
│      + run(reporter) -> SimulationResult                         │
└───────────────────────────────────────────────────────────────────┘

┌─ simulation ─────────────────────────────────────────────────────┐
│  «interface» IFluxDiscretization      PeacemanWellModel           │
│      + build(model) -> DiscretizedGrid   + build_connections(model)│
│           ▲                               │                       │
│           │ implements                    ▼                       │
│  TwoPointFluxDiscretization (TPFA,   WellConnection[]             │
│  DEFOLT) ──── future MPFAODiscretization (HƏLƏ YOXDUR)            │
│           │                                                        │
│           ▼                                                        │
│      DiscretizedGrid                                               │
│        + connections, pore_volume, cell_volume                     │
│        + compute_flux(d_phi)  ← ResidualAssembler BUNU çağırır,   │
│          TPFA-nın özünü DEYİL (bax §5.13)                          │
│                                                                    │
│  ImpesEngine ──implements──> ISimulationEngine                    │
│      - model: ReservoirModel        « inject »                    │
│      - config: SimulationConfig     « inject »                    │
│      - relperm: IRelativePermeabilityProvider   « inject »        │
│      - linear_solver: ILinearSolver             « inject »        │
│      - pvt: IPVTProvider = None                 « inject »        │
│      - capillary: ICapillaryPressureProvider = None               │
│      - initialization: IInitializationProvider = None             │
│      + run(reporter) -> SimulationResult                          │
│                                                                    │
│  ScipyCgIluSolver ──implements──> ILinearSolver                   │
│  CoreyRelativePermeabilityAdapter ──implements──>                 │
│                              IRelativePermeabilityProvider        │
└───────────────────────────────────────────────────────────────────┘

┌─ application ────────────────────────────────────────────────────┐
│  Project                          SimulationConfig                │
│   + geological_models{}            + end_time                     │
│   + reservoir_models{}             + time_stepping                │
│   + runs{}                         + linear_solver                │
│   + new_run()                      + output                       │
│                                    + validate()                   │
│  ReservoirModelBuilder                                            │
│   + build(geo, wells, fluids, scal, initial) -> ReservoirModel    │
│                                                                    │
│  SimulationService                                                │
│   - relperm_provider, linear_solver, pvt_provider, ...            │
│   + create_engine(model, config)                                  │
│   + run(model, config, reporter) -> SimulationResult              │
│   + run_in_project(project, model_name, config)                   │
└───────────────────────────────────────────────────────────────────┘

┌─ rendering (Qt YOXDUR) ──────────────────────────────────────────┐
│  MapRenderer · ProductionCurveRenderer · ScalRenderer             │
│  ValidationRenderer                                               │
│      + draw(axes, ...)   — hazır Axes qəbul edir, çəkir           │
└───────────────────────────────────────────────────────────────────┘

┌─ ui ─────────────────────────────────────────────────────────────┐
│  MainWindow                                                       │
│   - project, service, geology_builder, model_builder  « inject »  │
│   + rebuild_model()   — panellərdən domain obyektləri toplayır    │
│                                                                    │
│  GridGeometryPanel · RockFluidPanel · ScalPanel ·                 │
│  WellPanel · NumericalPanel                                       │
│      + values() -> domain / config obyekti                        │
│                                                                    │
│  SimulationWorker (QThread)                                       │
│  QtProgressReporter ──implements──> IProgressReporter             │
└───────────────────────────────────────────────────────────────────┘
```

---

## 4. UML — iş axını (ardıcıllıq diaqramı)

```
İstifadəçi   MainWindow   GeologyBuilder  ModelBuilder  SimulationService  ImpesEngine
    │            │              │               │              │               │
    │ parametr   │              │               │              │               │
    ├───────────>│              │               │              │               │
    │            │ build()      │               │              │               │
    │            ├─────────────>│               │              │               │
    │            │<─────────────┤ GeologicalModel               │               │
    │            │                              │              │               │
    │            │ build(geo, wells, fluids, scal, initial)     │               │
    │            ├─────────────────────────────>│              │               │
    │            │<─────────────────────────────┤ ReservoirModel               │
    │            │                                             │               │
    │ "İşə sal"  │                                             │               │
    ├───────────>│ create_engine(model, config)                │               │
    │            ├────────────────────────────────────────────>│               │
    │            │                                             │ validate()    │
    │            │                                             ├──> model      │
    │            │                                             │ ImpesEngine(  │
    │            │                                             │   model,      │
    │            │                                             │   config,     │
    │            │                                             │   providers)  │
    │            │                                             ├──────────────>│
    │            │<────────────────────────────────────────────┤               │
    │            │ SimulationWorker.start()                                    │
    │            ├────────────────────────────────────────────────────────────>│
    │            │              progress(fraction, message)                    │
    │            │<────────────────────────────────────────────────────────────┤
    │            │              SimulationResult                               │
    │            │<────────────────────────────────────────────────────────────┤
    │            │ renderer.draw(axes, result)                                 │
    │<───────────┤                                                             │
```

---

## 5. Dizayn qərarları və səbəbləri

### 5.1 Grid topologiyası ilə həndəsə ayrıldı
`CartesianGrid` yalnız hüceyrələrin qonşuluğunu bilir; ölçülər
`CellGeometry`-dədir. **Səbəb:** corner-point və ya qeyri-struktur grid
gələndə yalnız həndəsə sinfi əvəz olunacaq. Diskretizasiya kodu
`face_areas()` və `face_half_distances()` metodlarından istifadə edir,
`dx`-i birbaşa oxumur.

> **Bu qərar özünü doğrultdu** (Phase 5D, §5.18): `CornerPointGeometry`
> məhz həmin müqaviləni tətbiq edən bir alt sinif kimi əlavə olundu və
> nə topologiya (`CartesianGrid`), nə də diskretizasiya kodu dəyişdi.

### 5.2 Geoloji model ilə rezervuar modeli ayrıldı
| Geoloji model | Rezervuar modeli |
|---|---|
| Grid, həndəsə, xassə xəritələri | + flüidlər, quyular, SCAL, ilkin şərtlər |
| Horizontlar, faylar (həndəsə ilə) | + onlara **istinadlar** |
| Regionlar | + vahid sistemi |
| Geoloqun məhsulu | Rezervuar mühəndisinin məhsulu |

**Səbəb:** bir geoloji modeldən bir neçə rezervuar modeli (fərqli quyu
sxemi, fərqli SCAL) qurmaq mümkün olmalıdır. Həmçinin geoloji model
xarici mənbədən (sənin 3D Lay Modeli tətbiqindən, GRDECL faylından)
gələ bilər — o zaman yalnız yeni "import" sinfi yazılacaq.

### 5.3 Fay və horizont HƏNDƏSƏSİ deyil, İSTİNADI saxlanılır
`FaultReference` yalnız ad, mənbə identifikatoru və transmissivlik
çarpanı daşıyır. **Səbəb:** fayın həndəsəsi geoloji məlumatdır və
simulyasiyaya yalnız üz transmissivliyi kimi təsir edir. Bu ayrılıq
məlumatın iki yerdə saxlanmasının qarşısını alır.

### 5.4 ReservoirModel yeganə həqiqət mənbəyidir
Məhərrik konstruktorda hazır model alır və heç nə qurmur. **Səbəb:**
- eyni model müxtəlif konfiqurasiyalarla təkrar işlədilə bilər
- model yadda saxlanıla, müqayisə edilə, faylla mübadilə oluna bilər
- test yazmaq mümkün olur (model qur → mühərrikə ver → nəticəni yoxla)

### 5.5 Provider-lər üçün dependency injection
Məhərrik `IRelativePermeabilityProvider` interfeysindən asılıdır,
`CoreyParameters` sinfindən yox. **Səbəb (Dependency Inversion):**
region əsaslı SCAL cədvəlləri, histerezis və ya maşın öyrənmə əsaslı
model gələndə mühərrikin bir sətri belə dəyişməyəcək.

PVT, kapilyar təzyiq və initialization provider-ləri həmin vaxt
**inject edilə bilər, lakin implementasiya edilməmişdi**. Bu, artıq
CARİ deyil: hər üçü sonradan implementasiya edildi —
`BlackOilPVTProvider` (`imex2d/simulation/pvt/black_oil.py`),
`BrooksCoreyCapillaryProvider` (`imex2d/simulation/capillary.py`) və
onun cədvəl-əsaslı qardaşı `TableCapillaryPressureProvider`
(`imex2d/simulation/scal_tables_provider.py`),
`EquilibriumInitializationProvider`
(`imex2d/simulation/initialization/equilibrium.py`). Məhərrik onlar
`None` olanda YENƏ DƏ modelin statik dəyərləri ilə işləyir (geriyə
uyğunluq qorunur) — Dependency Inversion qərarının özü dəyişməyib,
sadəcə interfeyslər artıq boş deyil.

### 5.6 Corey adapteri müvəqqətidir
`CoreyRelativePermeabilityAdapter` yeni modul deyil — mövcud düsturları
interfeysə bağlayan nazik təbəqədir. Real (cədvəl-əsaslı) SCAL modulu
(`TableRelativePermeabilityProvider`, bax `imex2d/simulation/
scal_tables_provider.py`) artıq yazılıb, AMMA Corey adapteri **silinmədi**
— sadə modellər üçün hələ də faydalıdır və bütün mövcud testlər ondan
asılıdır. Model `scal_tables` təyin edilibsə cədvəl provideri, olmasa
Corey adapteri işlədir (ikisi eyni interfeysi paylaşır, bax `SCAL.md`).
Alternativ (mühərrikin `CoreyParameters`-i birbaşa oxuması) interfeysi
mənasız edərdi.

### 5.7 Bütün yoxlamalar modelə köçürüldü
`ReservoirModel.validate()` və `SimulationConfig.validate()` siyahı
qaytarır, `QMessageBox` çağırmır. **Səbəb:** eyni yoxlama həm interfeysdə,
həm skriptdə, həm testdə işləməlidir. UI yalnız siyahını göstərir.

### 5.8 Xətti həlledici ayrıca sinif oldu
`ILinearSolver` + `ScipyCgIluSolver`. **Səbəb:** AMG, PETSc və ya GPU
həlledicisi əlavə edəndə mühərrikə toxunulmayacaq. Ön-şərtçinin
ömrü də artıq həlledicinin daxili işidir, mühərrikin yox.

### 5.9 Konfiqurasiya obyektləri
`LinearSolverConfig`, `TimeSteppingConfig`, `OutputConfig`,
`SimulationConfig`. **Səbəb:** əvvəl `rtol=1e-8` mühərrikin içində,
`cfl=0.45` başqa dataclass-da, `20 gün` isə UI spinbox-unda idi.
İndi hər sabitin bir sahibi var.

### 5.10 Rendering Qt-dən ayrıldı
Renderer-lər hazır `matplotlib.Axes` qəbul edir. **Səbəb:** eyni kod
PDF hesabatda, Jupyter dəftərində və avtomatik testdə işləyir. Yuxarıdakı
smoke testdə bütün renderer-lər `Agg` backend ilə, Qt olmadan işlədi.

### 5.11 UI panelləri domain obyekti qaytarır
`ScalPanel.values() -> CoreyParameters`, `NumericalPanel.simulation_config()
-> SimulationConfig`. **Səbəb:** panel ilə mühərrik arasında "sərbəst
ədədlər" ötürülmür. Yeni giriş forması (məsələn fayl oxuyucusu) eyni
obyektləri qursa, qalan hər şey işləyir.

### 5.12 Composition root — `app.py`
Bütün `new` çağırışları burada. **Səbəb:** PVT modulu yazılanda dəyişəcək
yeganə fayl budur:
```python
service = SimulationService(
    relperm_provider=TableBasedRelPerm(...),
    pvt_provider=BlackOilPVT(...),      # ← yalnız bu sətir əlavə olunacaq
)
```

### 5.13 Axın diskretizasiyası: TPFA (indiki) + gələcək MPFA-O
("Numerical Discretization Architecture Preparation" fazası — bax
`imex2d/simulation/discretization.py` modul docstring-i tam detal üçün.)

```text
Flow Solver (FullyImplicitEngine / ImpesEngine)
    ↓
Residual Assembly (ResidualAssembler)
    ↓  yalnız grid.connections / grid.pore_volume / grid.compute_flux()-a güvənir
Flux Discretization Interface (IFluxDiscretization)
    ↓                                              ↓
TPFA (TwoPointFluxDiscretization, İNDİKİ)     Future MPFA-O (HƏLƏ YOXDUR)
```

**Nə edildi**: `IFluxDiscretization` (bax `imex2d/interfaces/discretization.py`)
— tək abstrakt metod `build(model) -> DiscretizedGrid`.
`TwoPointFluxDiscretization` bunu tətbiq edir, riyaziyyatı BİRƏBİR
eynidir (harmonik orta transmissivlik). `DiscretizedGrid`-ə YENİ bir
metod (`compute_flux(d_phi)`) əlavə olundu — `ResidualAssembler.
face_fluxes()` artıq `self.transmissibility * ΔΦ`-ni BİRBAŞA yazmır,
`self.grid.compute_flux(ΔΦ)` çağırır. Mobilite/upstream çəkiləndirmə
YENƏ DƏ `ResidualAssembler`-dədir, çünki bu, diskretizasiya sxemindən
ASILI OLMAYAN fizikadır. `FullyImplicitEngine`/`ImpesEngine`
konstruktoruna `flux_discretization: Optional[IFluxDiscretization] = None`
əlavə olundu — verilməyəndə `default_flux_discretization()` (= TPFA)
işlədilir, DEFOLT DAVRANIŞ dəyişmir.

**Tenzor permeabilitə (Phase 2 — "Full Tensor Permeability Implementation")**:
`RockProperties.permeability_tensor: Optional[PermeabilityTensor] = None`
— Kxx/Kyy/Kzz (məcburi) + Kxy/Kxz/Kyz (istəyə görə), hər biri `PropertyMap`
(hüceyrə-hüceyrə dəyişə bilər, bax audit §8). Yalnız 6 MÜSTƏQİL komponent
saxlanılır — simmetrik cütlər (`Kyx=Kxy` və s.) `as_matrices()`-də hər
çağırışda RİYAZİ bərpa olunur, ayrıca TƏKRAR saxlanmır.

*Doğrulama* (`PermeabilityTensor.validate()`): NaN/Inf hər komponentdə
tutulur; müsbət-müəyyənlik YALNIZ diaqonalın müsbətliyinə görə YOX,
`np.linalg.eigvalsh`-lə hesablanan HƏQİQİ məxsusi qiymətlərə görə
yoxlanılır (λ_min(K) > 0) — vektorlaşdırılıb (bütün hüceyrələr TƏK LAPACK
çağırışında, O(N)), Python dövrü yoxdur. Güclü anizotropluq (10⁴-ə qədər
sınanıb) süni kəsilmə OLMADAN qəbul edilir; etibarsız tenzor "təmir"
edilmir, sadəcə RƏDD edilir (`RockProperties.validate()` bunu ötürür).

*Fırlanma* (`rotate(R)`): `K_rot = R·K·Rᵀ`, `R`-in ORTOQONALLIĞI (`R·Rᵀ=I`)
yoxlanılır — əks halda `ValueError` (qeyri-ortoqonal "fırlanma" məxsusi
qiymətləri poza bilər). Testlərlə təsdiqlənib: izotrop tenzor izotrop
qalır, məxsusi qiymətlər saxlanılır, müsbət-müəyyənlik qorunur.

*Vahid çevirməsi* (`convert_units(from, to)`): mövcud `unit_conversions.
convert`-i BÜTÜN 6 komponentə EYNİ amillə tətbiq edir. **Məhdudiyyət**:
CSV/GRDECL idxal sərhədində tenzor komponentlərini AVTOMATİK çevirən
boru xətti HƏLƏ YOXDUR (çünki tenzor K-nı geologiya boru xəttindən
oxuyan mexanizmin ÖZÜ də yoxdur, bax §13) — `convert_units()` istifadəçi
kodunun əlində olan, hazır, düzgün alətdir, avtomatik ÇAĞIRILAN yer YOXDUR.

*TPFA münasibəti*: off-diaqonal komponent aşkarlansa, TPFA onu SƏSSİZCƏ
diaqonala yumşaltmır (`has_off_diagonal()`) — `DiscretizedGrid.warnings`-ə
HANSI komponentin (Kxy/Kxz/Kyz) və neçə hüceyrənin təsirləndiyini AÇIQ
göstərən xəbərdarlıq yazır, öz nəticəsini (transmissivliyini) DƏYİŞDİRMƏDƏN
(yalnız diaqonaldan istifadəyə davam edir).

*Serializasiya*: `imex2d/application/serialization.py` — bütün 6 komponent
ayrıca saxlanılır/bərpa olunur (`_permeability_tensor_to_dict`/`_from_dict`),
`None` (köhnə/tenzor-siz model) geriyə-uyğun qalır.

**DOĞRU İFADƏ (audit §21/Phase 2 §L)**: "tam tenzor permeabilite TƏMSİL
OLUNUR VƏ DOĞRULANIR, amma tam tenzor axın diskretizasiyası MPFA TƏLƏB
EDİR" — "simulyator tərəfindən tam dəstəklənir" YANLIŞDIR, İSTİFADƏ
EDİLMİR.

**Həndəsə hazırlığı**: `CellGeometry`-yə üç YENİ, əlavə metod —
`cell_centroid()`, `face_centroid(conn)`, `face_normal(conn)` (bax
`imex2d/interfaces/geometry.py::IGridGeometry` — sənədləşdirilmiş
müqavilə, `CellGeometry` bunu formal İRSƏN ALMIR, çünki bu kod
bazasında domain dataclass-ları interfeys-inject edilən strategiyalar
deyil — yalnız metod adları uyğunlaşdırılıb). Mövcud `face_areas`/
`face_half_distances`/`volumes` DƏYİŞMƏYİB, TPFA bunlardan istifadəyə
davam edir.

**Jacobian inteqrasiya nöqtəsi (HƏLƏ DƏYİŞDİRİLMƏYİB)**: `JacobianAssembler.
_build_pattern()` (`implicit/jacobian.py`) hər üz üçün DƏQİQ 2-hüceyrəli
blok fərz edir; `_flux()` `transmissibility`-ni birbaşa oxuyub tək-cüt
∂ΔΦ/∂p=±1 fərziyyəsi ilə törəmə qurur. MPFA-O gələndə hər ikisi
N-hüceyrəli stensilə ümumiləşdirilməli olacaq — bax `discretization.py`
modul docstring-i, "Jacobian inteqrasiya nöqtəsi" bölməsi, tam sətir
istinadları ilə. **Bu fayl bu fazada dəyişdirilməyib.**

**Bilərəkdən EDİLMƏYƏNLƏR**: MPFA-O-nun özü (MPFA-L də), native
corner-point həndəsə, tenzor K-nın TPFA-da (və ya hər hansı digər axın
həlledicisində) FAKTİKİ istifadəsi, geologiya boru xəttindən (SGS/SIS/
interpolyasiya) tenzor K komponentlərinin AVTOMATİK doldurulması (co-
kriging və s.) — hamısı gələcək faza.

### 5.14 Ümumi (qeyri-ortoqonal) çoxüzlü həndəsə nüvəsi (Phase 3)
Yeni modul: `imex2d/domain/polyhedral_geometry.py::HexahedralCell`/`Face`
— HÜCEYRƏ-BAŞINA (per-cell), potensial qeyri-ortoqonal 8-təpəli hüceyrə
həndəsəsi. `CellGeometry`-dən (Kartezian, vektorlaşdırılmış, DƏYİŞMƏYİB)
FƏRQLİ, ONU ƏVƏZ ETMİR — bir səviyyə AŞAĞIDA yerləşən SAF riyazi nüvədir:

```text
Geometry NÜVƏSİ (HexahedralCell/Face — saf riyaziyyat, Connections-dan
                  MÜSTƏQİL, bax audit §18/§21)
    ↓
Grid Geometry (CellGeometry — Kartezian; gələcək CornerPointGeometry
               bu nüvəni İSTİFADƏ EDƏ BİLƏR, HƏLƏ YAZILMAYIB)
    ↓
Topology (Connections — HANSI hüceyrələr bağlıdır)
    ↓
Discretization (TPFA indiki, MPFA-O gələcək)
    ↓
Flow
```

**Həcm/mərkəz**: tetraedr-parçalanması (8 təpənin ortası olan daxili
istinad nöqtəsindən hər üzün üçbucaqlarına) — standart, qabarıq (VƏ
yüngül qeyri-qabarıq) çoxüzlülər üçün RİYAZİ CƏHƏTDƏN DÜZGÜN üsul.
Sıfır/mənfi həcm (degenerativ/tərs-yönümlü hüceyrə) AÇIQ RƏDD edilir.

**Üz sahəsi/mərkəzi/normalı**: fan-triangulyasiya (`(v0,v1,v2),
(v0,v2,v3)`) + sahə-çəkili orta. Tam müstəvi üz üçün DƏQİQ; əyri
(warped) üz üçün YALNIZ TƏXMİNİ (bax `Face.is_planar()`) — bu, SÜKUTLA
GİZLƏDİLMİR, `is_planar()` vasitəsilə AÇIQ yoxlanıla bilər.

**Fırlanma/qeyri-ortoqonallıq**: hüceyrə fırlandıqda həcm/kənar
uzunluğu/üz sahəsi QORUNUR, normal DÜZGÜN fırlanır (test edilib).
`non_orthogonality_angle(d_ij, n_f)` — YALNIZ DİAQNOSTİKA, TPFA-nı
DƏYİŞDİRMİR (bax audit §12).

**Doğrulama**: NaN/Inf, sıfır/mənfi həcm, sıfır/mənfi sahə, sıfır-
uzunluqlu normal — hamısı AÇIQ RƏDD edilir, HEÇ NƏ "təmir" olunmur.

**Serializasiya**: BU FAZADA YOXDUR — çünki `HexahedralCell`/`Face`
HƏLƏ HEÇ BİR SAXLANILAN modelə (`GeologicalModel`/`ReservoirModel`)
BAĞLANMAYIB (bax audit §25: yalnız istifadə olunan sahələr saxlanılır).
Gələcək `CornerPointGeometry` modelə bağlananda, onun ÖZ serializasiyası
YALNIZ təpələri saxlamalıdır (sahə/həcm/mərkəz/normal HƏMİŞƏ YENİDƏN
hesablana bilər, TƏKRAR saxlanmamalıdır).

**Sərhəd üzləri**: bu modul "sərhəd VS daxili" TƏSNİFATINI ETMİR (bu,
topologiya sualıdır, bax §18) — hər hüceyrənin HƏR üzü (sərhəddə olsun
ya olmasın) etibarlı sahə/mərkəz/normal daşıyır; hansı üzün sərhəd
olduğunu XARİCİ topologiya (`Connections`) müəyyən edir.

**Bilərəkdən EDİLMƏYƏNLƏR** (Phase 3-də): MPFA-O-nun bu nüvədən
istifadəsi, ümumi (Kartezian-olmayan) GRID-səviyyəli həndəsə sinfi,
native corner-point idxal, ACTNUM/qeyri-aktiv hüceyrə inteqrasiyası,
qeyri-qabarıq (non-convex) hüceyrələr üçün formal doğrulama, faylların
(COORD/ZCORN) bu nüvəyə bağlanması.

> **Sonrakı fazaların yenilədiyi bəndlər**: GRID-səviyyəli həndəsə
> §5.15-də (Phase 4), MPFA-O-nun istifadəsi §5.16/§5.17-də (Phase 5A/5B),
> native corner-point idxal və `COORD/ZCORN`-un bu nüvəyə bağlanması isə
> §5.18-də (Phase 5D) YAZILIB. ACTNUM-un MPFA yolunda dəstəyi və
> qeyri-qabarıq hüceyrələrin formal doğrulaması HƏLƏ gələcək fazadır.

### 5.15 Grid-səviyyəli ümumi həndəsə (Phase 4)
Yeni modul: `imex2d/domain/general_grid_geometry.py::GeneralGridGeometry` —
§5.14-dəki `HexahedralCell`/`Face` NÜVƏSİNİ ÇOXLU hüceyrəyə tətbiq edir,
`Connections` (topologiya) ilə DETERMİNİSTİK inteqrasiya ilə:

```text
Geometry NÜVƏSİ (HexahedralCell/Face — TƏK hüceyrə)
    ↓
Grid Topology (Connections — HANSI hüceyrələr bağlıdır)
    ↓
General Grid Geometry (BU FAZA — ÇOXLU hüceyrə, paylaşılan üzlər,
                        sərhəd üzləri, vektorlaşdırılmış massivlər)
    ↓
Flux Discretization (TPFA indiki, MPFA-O gələcək)
```

**Üz deduplikasiyası** (audit §5/§6): floating-point koordinat
müqayisəsinə güvənmir — `Connections.axis`-dan (0=X,1=Y,2=Z) istifadə
edərək DETERMİNİSTİK: `cell_a`-nın müsbət-ox yerli üzü === `cell_b`-nin
mənfi-ox yerli üzü (EYNİ fiziki üz, TƏK `Face` obyekti saxlanılır). O(N)
— heç bir geometrik axtarış yoxdur.

**Sərhəd üzləri**: HƏR hüceyrənin 6 yerli üzündən `Connections`-da
əhatə OLUNMAYANLAR avtomatik sərhəd üzü kimi qeydə alınır (`neighbor=
None`) — sərhəd/daxili təsnifatı ÜÇÜN xüsusi giriş TƏLƏB OLUNMUR.

**Doğrulanmış invariantlar** (bax `tests/test_general_grid_geometry.py`,
22 test): fırlanma/sürüşmə/miqyaslanma altında həcm/mərkəz/sahə/normal
düzgün dəyişir; daxili üzlərdə owner/neighbor sahə-mərkəz-normal
uyğunluğu; qapanma münasibəti `Σ(A_f·n_f)≈0` maşın-dəqiqliyində;
qonşuluq/hüceyrə-üz xəritələməsi `Connections`-la tam uyğun; NaN/Inf/
sıfır-həcm/sıfır-sahə AÇIQ rədd edilir; qurma N-də XƏTTİ miqyaslanır
(ölçülüb, test edilib).

**Serializasiya**: BU FAZADA DA YOXDUR — `GeneralGridGeometry` heç bir
saxlanılan modelə (`ReservoirModel`/`Project`) BAĞLANMAYIB (eyni
qərar/səbəb, bax §5.14-ün serializasiya qeydi).

**TPFA reqressiyası**: dəyişməz — bu modul heç bir mövcud sinif
tərəfindən ÇAĞIRILMIR, `test_five_spot_reference_case_reproduces_
legacy_results` (dəqiq addım sayı/OOIP/RF) dəyişmədən keçir.

**Bilərəkdən EDİLMƏYƏNLƏR** (Phase 4-də): MPFA-O-nun bu qatdan
istifadəsi, native corner-point idxal (COORD/ZCORN → bu qat), fay
transmissivlik çarpanlarının bu qata daxil edilməsi (bu, DİSKRETİZASİYA
işidir, HƏNDƏSƏ yox).

> **Sonrakı fazalar**: MPFA-O bu qatdan Phase 5A-da (§5.16) istifadə
> etməyə başladı; `COORD/ZCORN` → bu qat bağlantısı Phase 5D-də (§5.18)
> `hexahedral_vertices()` dispetçeri ilə qurulub. Fay çarpanları HƏLƏ
> gələcək fazadır.

### 5.16 MPFA-O riyazi nüvəsi (Phase 5A)

Yeni paket: `imex2d/discretization/` — **TAM tenzorlu, çoxnöqtəli**
MPFA-O (O-method) LOKAL diskretizasiya nüvəsi. TAM riyazi spesifikasiya:
**`docs/mpfa_o_phase5a.md`** (bu bölmə YALNIZ arxitektura xülasəsidir).

```text
Geometry  (domain/general_grid_geometry.py::GeneralGridGeometry)
    ↓
Interaction Region  (mpfa_o_interaction.py::MPFAOInteractionRegion)
    ↓
Lokal MPFA-O riyaziyyatı  (mpfa_o_local_system.py::MPFAOLocalSystem)
    ↓
Axın əmsalları  (mpfa_o.py::MPFAOCoefficients / MPFAODiscretization)
```

**Riyaziyyat.** Hər grid TƏPƏSİ ətrafında bölgə qurulur; bölgədəki hər
hüceyrə payında (sub-control volume) təzyiq XƏTTİDİR, `p = p_c + g·(x−x_c)`.
Qradiyent `g` həmin payın DƏQİQ 3 sub-üzündəki kəsilməzlik nöqtəsi
şərtindən çıxır: `D g = π − p_c·1`. Yarım-axın
`q = −Γ a_σᵀ K_c D⁻¹ (π − p_c·1)` — burada `K_c` **TAM 3×3 simmetrik
matrisdir**, heç bir mərhələdə `nᵀKn` skalyarına yığılMIR. Sub-üzlərdə
axın kəsilməzliyi lokal sistem verir (`C π = D p + E π_bnd`), həlli isə
ÇOXNÖQTƏLİ transmissivlik matrisini: `T_cell = F C⁻¹ D + G`.

**Tenzor K necə daxil olur**: `RockProperties.permeability_tensor`
(`PermeabilityTensor.as_matrices()`) — 6 komponentin HAMISI
(Kxx,Kyy,Kzz,Kxy,Kxz,Kyz). Tenzor verilməyibsə diaqonal
`diag(PERMX,PERMY,PERMZ)` qurulur və bu, AÇIQ xəbərdarlıqla bildirilir.
Etibarsız tenzor (NaN/Inf, qeyri-simmetrik, qeyri-SPD) RƏDD EDİLİR —
klipləmə/ε-əlavəsi/simmetrikləşdirmə YOXDUR.

**Stensil**: daxili üz üçün 3D-də ≤18 hüceyrə. K-ortoqonal Kartezian
gridd-də (η=1) formulyasiya ÖZÜ 2 nöqtəyə yığılır və TPFA düsturunu
`1e-11` nisbi dəqiqliklə bərpa edir — MPFA DAXİLƏN TPFA ÇAĞIRMIR.

**Sərhəd**: `MPFAOBoundaryClosure.DIRICHLET` (defolt) sərhəd `π`-lərini
XARİCİ giriş kimi saxlayır və `T_bnd` əmsallarını verir — fiziki BC
İCAD EDİLMİR, yalnız gələcək BC qatı üçün struktur. `NEUMANN_ZERO`
(axınsız) AÇIQ opt-in seçimdir.

**Diaqnostika**: hər bölgə üçün `condition_number`, `rank`,
`determinant`, `singular`, `ill_conditioned`, `max_sub_cell_condition`
(`MPFAORegionDiagnostics`); sinqulyar sistem AÇIQ `MPFAOSingular
SystemError` verir — səssiz requlyarizasiya və ya TPFA-ya keçid YOXDUR.

**Yaddaş/performans**: bölgə qurulması TAM TOPOLOJİDİR (koordinat
axtarışı yox) → `O(N)`; lokal həll bölgə başına SABİT ölçülüdür (≤12×12);
`T_cell`/`T_bnd` `scipy.sparse.csr_matrix`-dir (sətir başına ≤18 əmsal)
— sıx saxlama `O(N²)` yaddaş olardı.

**Status (DÜRÜST təsnifat)**:

| | vəziyyət |
|---|---|
| 2D (tək lay, məhdud 3D lay kimi) | **implement + validasiya edilib** |
| 3D (tam formulyasiya) | **implement + validasiya edilib** |
| Xətti/manufactured həll dəqiqliyi | **validasiya edilib** (≈1e-13, nisbi ≈5e-15) |
| Lokal konservasiya | **validasiya edilib** (≈3e-13) |
| Ortoqonal izotrop limit ≡ TPFA | **validasiya edilib** (nisbi 9e-15) |
| Fırlanmış tam tenzor / skew grid | **validasiya edilib** |
| Residual/Jacobian/Nyuton inteqrasiyası | **YOXDUR** (Phase 5B) |
| Mobilite/upstream/PVT | **YOXDUR** (Phase 5B) |
| Analitik törəmələr | **YOXDUR** (Phase 5C) |
| Corner-point/struktursuz bölgə qurucusu | **YOXDUR** (Phase 5D) |
| Fay (fault) çarpanları MPFA-da | **YOXDUR** (Phase 5D, AÇIQ xəbərdarlıq verilir) |

**TPFA reqressiyası**: `TwoPointFluxDiscretization` HEÇ DƏYİŞMƏYİB —
düsturu, nəticələri, `default_flux_discretization()` DEFOLTU eynidir.
`IFluxDiscretization`-a YALNIZ NON-ABSTRAKT, defoltu `False` olan
`supports_multipoint_stencil()` əlavə edilib (geriyə uyğunluq pozulmur).
`MPFADiscretizedGrid.compute_flux(d_phi)` QƏSDƏN `NotImplementedError`
verir — çoxnöqtəli axını tək üzün `ΔΦ`-sindən uydurmaq metodu TPFA-ya
çevirmək olardı.

**Testlər**: `tests/test_mpfa_o.py` — 67 test (A–R kateqoriyaları,
anti-pseudo-MPFA sübutu, əl ilə yoxlanıla bilən istinad halı,
performans/yaddaş reqressiyası).

### 5.17 MPFA-O qlobal residual/stensil inteqrasiyası (Phase 5B-1)

Phase 5A LOKAL nüvəni verirdi; bu faza onu **qalıq (residual)
səviyyəsinə** bağlayır. Spesifikasiya: **`docs/mpfa_o_phase5b1.md`**.

```text
Phase 5A: T_face (nface × ncell, csr)
    ↓  sətir seçimi (connection_faces — DETERMİNİSTİK biyeksiya)
T_conn (nconn × ncell)          ← mpfa_global.py::MPFAGlobalOperator
    ↓  Φ_α = p − Pc − ρ_α g D   (HÜCEYRƏ üzrə, §4)
q_pot,f = Σ_j T_conn[f,j] Φ_α,j
    ↓  upstream = sign(q_pot)   (HƏQİQİ çoxnöqtəli axından, §6)
q_α,f = q_pot,f · (λ_α/B_α)_upstream
    ↓  net_influx (DƏYİŞMƏYƏN np.add.at)
R_i = (akkumulyasiya)/Δt − influx − quyu     (DƏYİŞMİR)
```

**Yeni fayl**: `imex2d/discretization/mpfa_global.py::MPFAGlobalOperator`
— `T_conn`, divergensiya `D`, `A = D·T_conn` (SABİT mobilitəli XƏTTİ
operator), `global_stencil_pattern()`, `conservation_report()`.
Hamısı `scipy.sparse`; `np.zeros((ncell, ncell))` HEÇ BİR YERDƏ YOXDUR.

**`ResidualAssembler`-də DƏYİŞİKLİK — DAR budaqlanma**: `_multipoint`
bayrağı (duck-typed `supports_multipoint_stencil()`) + `cell_potentials()`
+ `_multipoint_face_fluxes()`. `net_influx`, `well_rates`,
`accumulation`, `residual`, `material_balance_error` TOXUNULMAYIB —
akkumulyasiya/PVT/relperm/quyu/zaman addımı MÜŞTƏRƏKDİR.

**İkiqat sayma (§8)**: sub-üz payları Phase 5A-da fiziki üzdə cəmlənir;
fiziki üz ↔ `Connections` girişi BİYEKSİYADIR; sərhəd üzləri
`NEUMANN_ZERO`-da sıfır axın daşıyır. Hər üz balansa DƏQİQ BİR DƏFƏ
daxil olur — testlə (üz-səviyyəli və əlaqə-səviyyəli balansın
müqayisəsi) yoxlanılır.

**Potensial konvensiyası (§4)**: TAM POTENSİAL (pressure-only DEYİL).
Sabitlər/sıxlıq/Pc mövcud koddan gəlir — İKİNCİ cazibə konvensiyası
YOXDUR. YEGANƏ (SƏNƏDLƏŞDİRİLMİŞ) fərq: TPFA üz-ortalanmış sıxlıq
(`½(ρ_a+ρ_b)`), MPFA isə hüceyrə sıxlığı işlədir (çoxnöqtəli stensildə
"həmin iki hüceyrə" yoxdur). Cazibəsiz VƏ YA bərabər sıxlıqda ikisi
BİRƏBİR eynidir; fərq `O(Δρ·ΔD)`-dir və testlə ÖLÇÜLÜR.

**Status (DÜRÜST təsnifat)**:

| | vəziyyət |
|---|---|
| Qlobal seyrək MPFA stensili (`T_conn`, `A`, naxış) | **implement + validasiya edilib** |
| Qalıq (residual) inteqrasiyası | **implement + validasiya edilib** |
| Mobilitə/upstream (çoxnöqtəli axının işarəsindən) | **implement + validasiya edilib** |
| Daxili üz + qlobal konservasiya | **validasiya edilib** (≈1e-12 / ≈6e-14) |
| Ortoqonal izotrop limitdə MPFA qalığı ≡ TPFA qalığı | **validasiya edilib** (nisbi 1.5e-16) |
| no-flow (Neumann-0) sərhəd | **dəstəklənir** |
| Dirichlet sərhəd | **struktur hazırdır**, qalıq qatı π ötürmür → RƏDD EDİLİR (5B-2) |
| Qeyri-xətti analitik Jacobian | **YOXDUR** (5B-2) — MPFA ilə `JacobianAssembler` AÇIQ imtina edir |
| Nyuton/`FullyImplicitEngine`/`ImpesEngine` MPFA ilə | **YOXDUR** (5B-2) — AÇIQ imtina |
| MPFA-da fay (fault) | **YOXDUR** (5D) — model faylıdırsa qalıq RƏDD EDİLİR |
| corner-point həndəsə (COORD/ZCORN → 8 təpə, dəqiq həcm/sahə/normal) | **implement edilib** (5D, bax §5.18) — MPFA-O `build()` təpələri birbaşa istehlak edir |
| ACTNUM (MPFA yolunda) | **YOXDUR** (5D) — model qeyri-aktiv hüceyrə daşıyırsa qalıq AÇIQ RƏDD EDİLİR; TPFA yolu ACTNUM-u tam dəstəkləyir |

**Simulyator HƏLƏ MPFA ilə qeyri-xətti hasilat proqnozu VERƏ BİLMİR** —
Phase 5B-1 YALNIZ qalıq qiymətləndirməsidir.

**TPFA reqressiyası**: `TwoPointFluxDiscretization` və TPFA-nın
qalıq/Jacobian/Newton yolu HEÇ DƏYİŞMƏYİB;
`default_flux_discretization()` HƏLƏ DƏ TPFA qaytarır.

**Testlər**: `tests/test_mpfa_o_global_assembly.py` — 57 test
(A–T kateqoriyaları, lokal→qlobal uyğunluq, qlobal anti-pseudo-MPFA,
ikiqat-sayma detektoru, əl ilə yoxlanıla bilən qlobal matris).

### 5.18 Native corner-point həndəsə (Phase 5D)

`imex2d/domain/corner_point_geometry.py` — §5.14-də "gələcək faza"
kimi qeyd olunan **native corner-point idxalı** artıq yazılıb.
`COORD` (pillar oxları) və `ZCORN` (künc dərinlikləri) daha bərabər
bloka approksimasiya EDİLMİR; hər hüceyrə üçün 8 təpə açıq qurulur və
bütün ölçülər həmin təpələrdən hesablanır.

    COORD + ZCORN
        ↓  pillar boyunca xətti interpolyasiya  (corner_point_nodes)
    (ncell, 8, 3) təpələr
        ↓                                    ↓
    CornerPointGeometry                 GeneralGridGeometry
    (CellGeometry alt sinfi —           (üz-əsaslı, owner/neighbor —
     TPFA/hesabat/görüntü)               MPFA-O, §5.15)

**Niyə `CellGeometry`-nin ALT SİNFİ**: §5.1-dəki qərar ("topologiya ilə
həndəsə ayrıldı, diskretizasiya `face_areas()`/`face_half_distances()`
işlədir, `dx`-i birbaşa oxumur") məhz bu anı nəzərdə tuturdu. Ona görə
yeni sinif həmin müqaviləni tətbiq edir və model zəncirinin qalanı
(`GeologicalModel` → `ReservoirModelBuilder` → `ReservoirModel` → TPFA)
BİR SƏTİR belə dəyişmir.

**Vektorlaşdırma**: §5.14 nüvəsi (`HexahedralCell`/`Face`) hüceyrə
başına Python obyektidir. `corner_point_geometry.py` EYNİ riyaziyyatın
tam vektorlaşdırılmış (bütün grid üçün tək numpy keçidi) formasıdır;
ikisinin ədədi uyğunluğu testlə qıfıllanıb, ona görə iki kod yolu
ayrıla bilməz.

**Üz parçalanmasında düzəliş**: `Face._triangles()` `v0`-dan sadə fan
əvəzinə MƏRKƏZ-fan işlədir. Səbəb: sadə fan üzün bir diaqonalını seçir,
paylaşılan üz isə owner və neighbor tərəfindən FƏRQLİ təpədən başlanır
(`X+` = `(1,2,6,5)` vs `X-` = `(3,0,4,7)`) — ƏYRİ üzdə iki hüceyrə eyni
fiziki üzə fərqli sahə/mərkəz verərdi, yəni hüceyrələr arasında boşluq/
örtüşmə (qeyri-konservativ diskretizasiya). Mərkəz-fan simmetrikdir.
MÜSTƏVİ üzlərdə nəticə eynidir, ona görə Kartezian nəticələr dəyişmir.

**Diaqnostika** (heç biri sükutla keçmir): sol-əlli deck bir dəfəlik
qlobal sarğı düzəlişi ilə, dejenerativ pillar / pinch-out / mənfi həcmli
hüceyrələr isə sayılıb `DiagnosticReport`-a yazılır.

**Testlər**: `tests/test_corner_point_geometry.py` — 44 test
(rekonstruksiya, dəqiq həcm, əyri üz sahəsi, həqiqi normal, MPFA-O/TPFA
inteqrasiyası, Kartezian geriyə uyğunluğu, yönüm/dejenerativ hallar).

---

## 6. Növbəti modulun necə qoşulacağı

Bu bölmə refaktorinqin sübutu kimi PVT modulunun əlavəsini nümunə
göstərirdi — **bu iş artıq edilib** (`imex2d/simulation/pvt/
black_oil.py::BlackOilPVTProvider`, bax 5.5). Faktiki qoşulma
planlaşdırılan üç addımla üst-üstə düşdü: yeni provider sinfi yazıldı,
mühərrik onu `None` olanda köhnə sabit dəyərlərlə, veriləndə isə
provider vasitəsilə işlədir, `app.py`/`SimulationService`-də inject
edilir. `domain`, `rendering`, `ui` qatlarına toxunulmadı — refaktorinqin
məqsədi olan qat ayrılığı bu əlavə ilə TƏSDİQLƏNDİ.

---

## 7. İşə salmaq

```bash
pip install PyQt5 matplotlib numpy scipy
# OPM Flow nəticələrinin idxalı üçün (istəyə bağlı, bax OPM_IMPORT.md):
pip install resdata pytest
# ResInsight tipli sürətli 3D görüntü üçün (istəyə bağlı, bax VISUALIZATION.md):
pip install vtk
python app.py
```

Skript rejimi (interfeys olmadan):

```python
from imex2d.application.scenarios import SyntheticGeologicalModelBuilder, five_spot
from imex2d.application.model_builder import ReservoirModelBuilder
from imex2d.application.config import SimulationConfig
from imex2d.application.simulation_service import SimulationService
from imex2d.domain.scal import CoreyParameters
from imex2d.simulation.scal_adapter import CoreyRelativePermeabilityAdapter

geology = SyntheticGeologicalModelBuilder().build(41, 41, 20, 20, 10, 0.22, 150.0)
scal = CoreyParameters()
model = ReservoirModelBuilder().build(geology, five_spot(geology.grid), scal=scal)
service = SimulationService(CoreyRelativePermeabilityAdapter(scal))
result = service.run(model, SimulationConfig(end_time=1500))
print(result.final_recovery_factor)
```
