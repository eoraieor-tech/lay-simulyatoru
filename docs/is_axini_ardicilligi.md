# İş axını: hansı düstur hansı ardıcıllıqla işləyir

**Yazılıb:** 25 sentyabr 2026 (Seans 47) · **Cütü:** [nezeri_esaslar.md](nezeri_esaslar.md)

[İnteraktiv xəritənin](teqdimat/README.md) («İş axını» tabı) hər bloku üçün
**icra ardıcıllığı**: proqram həmin blokda hansı düsturu hansı sıra ilə
hesablayır. Ardıcıllıq yaddaşdan yazılmayıb — kodun icra yolundan
çıxarılıb (`geology_service.build`, `implicit/engine.run`,
`implicit/newton.solve`, `implicit/residual.residual`,
`impes_engine._solve_pressure`).

`nezeri_esaslar.md` sualı «bu düstur haradan gəlir» idi; bu sənədin sualı
**«bu düstur nə vaxt, nədən sonra, nəyin içində hesablanır»**.

---

## 0 · Ümumi ardıcıllıq — bir baxışda

```
GİRİŞ            quyu CSV · grid · faylar · PVT · SCAL · quyular · müşahidə
  ↓ (bir dəfə)
GEOLOGİYA        QC → çevirmə → variogram → anizotropluq → qonşuluq
                 → kriging/SGS/SIS → geri çevirmə → çarpaz-doğrulama
  ↓
MODEL            GeologicalModel → ReservoirModel (+PVT cədvəli, SCAL, quyular)
  ↓ (bir dəfə, mühərrik qurulanda)
ÖN-HESABLAMA     T_ij (TPFA/MPFA-O) · PV · WI (Peaceman) · ilkin p, Sw, Sg
  ↓
ZAMAN DÖVRƏSİ    ┌─ hər addım: rate payları → Nyuton → THP → səth debiti
                 │              → BHP limiti → qəbul → sıralar → adaptiv Δt
                 │   ┌─ hər Nyuton iterasiyası: PVT → kr → Pc → Φ → upstream
                 │   │   → axın → akkumulyasiya → quyu → R → CNV/MB → J → δ
                 └───┴─────────────────────────────────────────────────────
  ↓ (post-proses)
NƏTİCƏ           THP traversi → sıralar/anlar → qrafik · ixrac · doğrulama
                 → (istəsən) history matching: bütün dövrə N dəfə təkrarlanır
```

Ən vacib iki ardıcıllıq — §4 (bir Nyuton iterasiyasının içi) və
§5 (IMPES addımı). Qalan hər şey onların ətrafındadır.

---

## 1 · Geologiya bloku — «quyudan xəritəyə»

`WellBasedGeologicalModelBuilder.build()`-un FAKTİKİ sırası. Sıra
təsadüfi deyil: hər addım özündən əvvəlkinin nəticəsini tələb edir.

| # | Addım | Düstur / əməliyyat | Nə üçün məhz burada |
|---|---|---|---|
| 1 | Giriş yoxlaması | `dataset.validate()` | Yararsız data variogramı səssizcə korlayardı |
| 2 | Keyfiyyət (QC) | təkrar koordinat siyasəti, kənar dəyər balı, paylanma təhlili | Kənar dəyər `sill`-i şişirdir → bütün xəritə dəyişir |
| 3 | Grid + **müvəqqəti** həndəsə | `CartesianGrid`, düz qutu | Struktur hələ bilinmir |
| 4 | **Struktur səthlər** (TOP/BOTTOM) | AREAL interpolyasiya (yalnız X,Y; Z = 0) | **Sıra kritikdir:** qalan xassələr hüceyrə mərkəzinin Z-ini işlədir, ona görə həqiqi həndəsə ƏVVƏL qurulmalıdır |
| 5 | Həqiqi həndəsə | `_build_structural_geometry` → `model.geometry` | Bundan sonra `targets` yenidən hesablanır |
| 6 | **Kateqorik** sahələr (fasiya) | SIS və ya deterministik indikator kriging: `I_k(x) ∈ {0,1}` → kriging → ehtimal → normallaşdırma → nümunə | Kəsilməzlərdən ƏVVƏL, çünki fasiya-şərtli SGS ona istinad edir |
| 7 | Lay mövcudluğu | `availability` — hansı xassə hansı K-təbəqəsində var | Məlumatsız laya uydurma dəyər yazılmasın |
| 8 | **Kəsilməz** xassələr | aşağıdakı 8 addımlıq zəncir (hər xassə üçün ayrıca) | — |
| 9 | Çarpaz-doğrulama | leave-one-out / k-fold → RMSE, bias → etibar kateqoriyası | Nəticə hazırdır, indi ona qiymət verilir |

### 1.1 Bir kəsilməz xassənin içindəki zəncir (8 addım)

```
1. STRATEGİYA        ad → PropertyStrategy (PERMX → log-kriging, PORO → logit…)
2. ÇEVİRMƏ           Y = ln(Z + offset)          (keçiricilik)
                     Y = logit((Z − a)/(b − a))  (hədli: PORO, NTG)
                     Y = norm.ppf((rank − 0.5)/n) (SGS üçün normal-score)
3. DENEYSEL VARİOGRAM   γ̂(h) = 1/(2N(h)) · Σ (Y(xᵢ) − Y(xᵢ+h))²
4. MODEL FİT         γ(h) = nugget + sill·g(h/a) → least_squares ilə
                     (sferik / eksponensial / qauss; ən yaxşısı seçilir)
5. DOĞRULAMA         yararsız parametr Kriging-ə ÇATMIR (ValueError)
6. ANİZOTROPLUQ      x'' = S·R_dip·R₀·x  →  d_ani = ‖x''ᵢ − x''ⱼ‖
7. QONŞULUQ          cKDTree ilə anizotrop radiusda N ən yaxın nöqtə
8. YERLİ KRİGİNG     [Γ 1; 1ᵀ 0]·[w; μ] = [γ₀; 1]
                     Ŷ(x₀) = Σ wᵢYᵢ        σ²_K = Σ wᵢγ₀ᵢ + μ
```

Sonra **geri çevirmə** (`exp`, `logit⁻¹`, tərs normal-score) və fiziki
hədd tətbiq olunur.

**SGS seçilibsə** 8-ci addımdan sonra bir addım ARTIQ olur — və məhz bu,
SGS-i Kriging-dən ayırır:

```
8a. NÜMUNƏ    Y*(x₀) ~ N( Ŷ(x₀), σ²_K(x₀) )      ← kriging qiyməti YAZILMIR
8b. ƏLAVƏ     Y*(x₀) kondisioner çoxluğa qoşulur  ← növbəti hüceyrə onu görür
```

Ardıcıllıq təsadüfi yol (random path) ilə gedir, ona görə hər toxum (seed)
ayrı realizasiya verir.

---

## 2 · Model bloku — ReservoirModel yığılır

```
GeologicalModel (grid, həndəsə, xassə xəritələri, regionlar)
   + quyular (perforasiya, rejim, limit, lülə həndəsəsi)
   + PVT     (korrelyasiya cədvəli VƏ YA deck cədvəli)
   + SCAL    (Corey parametrləri VƏ YA SWOF/SGOF cədvəli)
   + ilkin şərtlər (datum p, kontaktlar)
   ↓
ReservoirModel.validate()  →  xəta siyahısı; natamam model mühərriyə KEÇMİR
```

**PVT cədvəli məhz burada qurulur** (mühərrik korrelyasiyanı heç vaxt
görmür). Daxili sıra:

```
1. Pb        verilməyibsə: Pb = 0.6·p_max  (və ya Standing-in tərsi ilə Rsb-dən)
2. ANCHOR    Rsb = Standing(Pb) → Bob = VB(Rsb) → μob = BR(μod, Rsb)
3. ŞAXƏ      p < Pb:  Rs = Standing(p), Bo = VB(Rs(p)), μo = BR(μod, Rs(p))
             p ≥ Pb:  Rs = Rsb,  Bo = Bob·exp(co(Pb − p)),  μo = μob·(p/Pb)^0.278
4. SU        Bw = McCain(p,T),  μw = Meehan(T, duzluluq)
5. QAZ       Tpc,Ppc = Sutton(γg) → Z = Beggs-Brill(ppr,Tpr)
             → Bg = (Z·T/p)/(Zsc·Tsc/psc),  μg = Lee-Gonzalez-Eakin(ρg)
6. CƏDVƏL    PVTTable — mühərrik bundan sonra YALNIZ np.interp işlədir
```

---

## 3 · Mühərrik qurulanda BİR DƏFƏ hesablananlar

Zaman dövrəsinə girməzdən əvvəl, `create_engine` daxilində:

| # | Kəmiyyət | Düstur | Qeyd |
|---|---|---|---|
| 1 | Üz transmissivliyi | `T = C·A / (d_a/k_a + d_b/k_b)` | harmonik orta; MPFA-O seçilibsə əvəzinə `T_cell[F,c]` matrisi qurulur |
| 2 | Fay çarpanı | `T ← T · multiplier` (sealing → 0) | üzə tətbiq olunur |
| 3 | Məsamə həcmi | `PV = V·φ` (NTG varsa ona vurulur) | `PV_ref` — sıxılma sonra tətbiq olunur |
| 4 | Quyu indeksi | `WI = 2π·C·√(K1K2)·h / (ln(re/rw) + S)`, `re` — anizotrop Peaceman | hər perforasiya üçün; ACTNUM = 0 olan hüceyrə SİYAHIYA DÜŞMÜR |
| 5 | İlkin təzyiq | `p(z) = p_datum + ρ_lay·g·(z − z_datum)/1e5` | `ρ_lay = ρ_səth/B`; qaz varsa canlı neft: `ρo = (ρo_səth + Rs·ρg_səth)/Bo` |
| 6 | İlkin doyumluluq | OWC-dən yuxarı: `Pc(D) = (ρw − ρo)g(D_owc − D)/1e5` → `Sw = Pc⁻¹` | Pc provider yoxdursa kəskin kontakt; GOC varsa qaz papağı |
| 7 | Provider seçimi | Corey ⟶ cədvəl ⟶ (qaz varsa) Stone II sarğısı | qaz fazası PVT-də varsa üç fazalı mühərrik AVTOMATİK seçilir |

Bu yeddisi zaman dövrəsində **bir də hesablanmır** (istisna: `PV(p)` süxur
sıxılması ilə hər iterasiyada yenilənir — §4, addım 4).

---

## 4 · Fully implicit: bir zaman addımı və bir Nyuton iterasiyası

### 4.1 Addımın xarici ardıcıllığı (`engine.run`)

```
1.  RATE paylanması      λ_bağlantı = λw + λo (yığılmış vəziyyətdən)
                         → rate_share hər perforasiyaya
2.  Səth debiti proqnozu surface_rate.predict(state)
3.  NYUTON               time_stepper.advance() → §4.2 dövrəsi
4.  THP dövrəsi          addımın ÖZ debitləri ilə IPR ∩ VLP → yeni BHP
                         → dəyişiklik > 1 bar-dırsa addım YENİDƏN həll olunur
5.  Səth debiti dövrəsi  sapma > 1e-3 → hədəf miqyaslanır → yenidən həll
6.  BHP limit dövrəsi    limit pozulubsa RATE → BHP keçidi → yenidən həll
7.  QƏBUL                state ← new_state,  t ← t + Δt
8.  SIRALAR              qo, qw, qwi, Np, Wp
                         WCT = qw/(qo+qw)·100,  RF = Np/OOIP·100,  p̄ = mean(p_aktiv)
9.  ANLAR                snapshot_interval-a çatdıqda 3D an yazılır
10. ADAPTİV Δt           az iterasiya → ×1.5;  çox → ×0.5;  yığılmadı → kəs və TƏKRARLA
```

4–6-cı dövrələr **qalığı və Jakobianı dəyişmir** — yalnız bağlantının
`mode`/`target`-ini dəyişib addımı yenidən həll edirlər.

### 4.2 Bir Nyuton iterasiyasının içi — ƏSAS ARDICILLIQ

`newton.solve` → `residual.residual`. Hər iterasiyada tam olaraq bu sıra:

```
1. FLÜİD VƏZİYYƏTİ          cari (p, Sw) → cədvəldən interpolyasiya
      μw(p), μo(p), Bw(p), Bo(p)            ← PVT (np.interp)
      krw(Sw), kro(Sw)                      ← Corey / cədvəl / Stone II
      Pc(Sw)                                ← Brooks-Corey / cədvəl
      λw = krw/μw,   λo = kro/μo            ← MOBİLLİK

2. AKKUMULYASİYA (yeni və köhnə qat)
      PV(p) = PV_ref·[1 + c_r·(p − p_ref)]  ← süxur sıxılması
      N_w = PV·Sw/Bw,     N_o = PV·(1−Sw)/Bo

3. ÜZ POTENSİALLARI
      ΔΦ_o = Δp − ½(ρo_a + ρo_b)·g·ΔD·1e−5
      ΔΦ_w = ΔΦ_o(ρw ilə) − (Pc_a − Pc_b)   ← kapilyar YALNIZ su qolunda

4. UPSTREAM SEÇİMİ          up = (ΔΦ ≥ 0) ? cell_a : cell_b     ← hər faza AYRI

5. ÜZ AXINI                 q_p = compute_flux(ΔΦ_p) · (λ_p/B_p)[up]
      (TPFA-da compute_flux = T·ΔΦ;  MPFA-O-da Σ_c T_cell[F,c]·p_c)

6. XALİS AXIN               hüceyrəyə: −q üzün A tərəfinə, +q B tərəfinə
                            (np.add.at → LOKAL saxlanma avtomatik)

7. QUYU DEBİTLƏRİ           BHP: q_p = WI·(λ_p/B_p)·(p_hüceyrə − BHP)
                            RATE: hədəf × rate_share (payı addımın əvvəlindən)

8. QALIQ                    R_p = (N_p^new − N_p^old)/Δt − influx_p − q_p

9. KONVERGENSİYA            CNV_p = max|R_p|·Δt/PV      ← lokal
                            MB_p  = |ΣR_p|·Δt/Σ(PV·S/B) ← qlobal
                            ikisi də tolerantlıqdadırsa → QAYIT (yığıldı)

10. JAKOBİAN                J = ∂R/∂x — ANALİTİK (sonlu fərq YOX)
                            törəmələr: ∂B/∂p, ∂μ/∂p, ∂kr/∂S, ∂Pc/∂S, ∂WI·λ/∂x

11. AKTİV DARALMA           ACTNUM = 0 olan hüceyrələr sistemdən çıxarılır

12. XƏTTİ HƏLL              J·δ = −R
                            BiCGStab + CPR ön-şərtçi (təzyiq ayrılır → ILU)
                            kiçik sistemdə birbaşa LU

13. YENİLƏMƏ (chop)         δp və δSw iterasiyaya görə HƏDDƏ kəsilir
                            x ← x + chop(δ)

14. XƏTTİ AXTARIŞ           addım qalığı artırırsa δ kiçildilir (line search)

   → 1-ə qayıt (maks. 20 iterasiya)
```

**Diqqət:** 1-ci addımdakı `previous_fluid` bir dəfə, dövrədən ƏVVƏL
hesablanır — köhnə qat iterasiya boyu dəyişmir.

### 4.3 Üç fazalı halda fərq

Yuxarıdakı 14 addım eyni qalır, üç yerdə genişlənir:

```
1-də  ƏLAVƏ:  Bg(p), μg(p), Rs(p) və ya Rs (3-cü dəyişəndən),  krg, Pcog
2-də  ƏLAVƏ:  N_g = PV·( Sg/Bg + So·Rs/Bo )         ← qazın İKİ mənbəyi
8-də  ƏLAVƏ:  üçüncü tənlik R_g
```

**Dəyişən keçid (Sg ↔ Rs) harada baş verir:** iterasiyanın İÇİNDƏ YOX —
yalnız Nyuton **yığılandan SONRA**, bir dəfə. Səbəb ölçülüb: hər
iterasiyada çağırılanda hüceyrə doymuş↔doymamış arasında dövr edir və
Nyuton heç vaxt yığılmır.

---

## 5 · IMPES: bir addımın ardıcıllığı

Tam fərqli sıra — burada Nyuton dövrəsi YOXDUR, hər kəmiyyət **köhnə**
vəziyyətdən oxunur:

```
1. FLÜİD (KÖHNƏ vəziyyətdən)   μw, μo, Bw, Bo, ct
2. MOBİLLİK                     λw, λo,  λt = λw + λo
3. POTENSİAL + UPSTREAM         ΔΦ_o, ΔΦ_w → up_o, up_w
4. ÜZ ƏMSALLARI                 T_o = T·λo[up_o],  T_w = T·λw[up_w],  TT = T_o + T_w
5. CAZİBƏ/KAPİLYAR              SAĞ TƏRƏFƏ keçir (implicit DEYİL):
                                q_grav+cap = T_o(ΔΦ_o − Δp) + T_w(ΔΦ_w − Δp)
6. AKKUMULYASİYA                acc = PV·ct/Δt        → matrisin diaqonalı
7. QUYU ÜZVLƏRİ                 BHP: a = WI·λ → diaqonala,  a·p_hədəf → RHS
                                RATE: birbaşa RHS-ə (payla)
8. XƏTTİ HƏLL                   (acc + ΣTT)·p = RHS    → CG + ILU (simmetrik!)
9. AXINLAR                      yeni təzyiqlə üz axınları
10. DOYUMLULUQ (EXPLICIT)       Sw ← Sw + Δt/PV · (net su axını + quyu)
11. CFL YOXLAMASI               Δt_CFL = min PV/(throughput·max|dfw/dSw|)
                                Δt > cfl_factor·Δt_CFL olsa → addım ATILIR,
                                Δt kiçildilir və TƏKRARLANIR
12. QƏBUL                       p, Sw yazılır, sıralar hesablanır
```

Fərqin mahiyyəti: fully implicit-də mobillik **yeni** qatdan (ona görə
Nyuton lazımdır), IMPES-də **köhnə** qatdan (ona görə CFL lazımdır).

---

## 6 · Post-proses — mühərrik bitəndən sonra

```
1. LÜLƏ HİDRAVLİKASI (lülə həndəsəsi verilmiş quyularda)
   hər addım və hər quyu üçün, quyu dibindən quyu başına, seqment-seqment:
       ρ_ns, v_m  → Re = ρvd/μ  → f = Chen(Re, ε/d)
       dp/dz = [ρ_qrav·g + f·ρ_ns·v_m²/(2d)]/1e5
       THP = BHP − Σ Δp_seqment           ← sabit-nöqtə iterasiyası ilə
   (THP İDARƏSİNDƏ isə bu, addım daxilində TƏRSİNƏ işləyir — §4.1 addım 4)

2. NƏTİCƏ SIRALARI     debit, WCT, GOR, RF, p̄, BHP/THP, 3D anlar
3. GÜNLÜK CƏDVƏL       addım sıraları gün-gün interpolyasiya olunur
4. ÇƏKMƏ / İXRAC       matplotlib · VTK · CSV/JSON · PDF · Eclipse deck · GIF
```

---

## 7 · Doğrulama blokunun ardıcıllığı

```
BUCKLEY-LEVERETT
  1. fw(Sw) = λw/(λw + λo)                 ← SCAL-dan
  2. dfw/dSw ədədi olaraq
  3. Welge: Swi-dən çəkilən toxunan → chord maksimumu → Sw_shock
  4. x(Sw) = (q·t)/(φ·A) · dfw/dSw          ← rarefaksiya yelpiyi
  5. x_front = velocity · slope_shock · t
  6. profil: [süpürülmüş | yelpik | sıçrayış | toxunulmamış]
  7. simulyasiya profili ilə müqayisə: RMS, cəbhə mövqeyi, həcm balansı

SPE1CASE2
  1. deck oxunur (PVTO/PVDG/PVTW/SWOF/SGOF) → model
  2. üç fazalı mühərrik işlədilir
  3. OPM Flow-un summary faylı oxunur (SMSPEC + UNSMRY)
  4. eyni vaxt nöqtələrinə interpolyasiya → fərq cədvəli
```

---

## 8 · History matching — bütün dövrənin təkrarı

Ən xarici dövrə. Hər bir qiymətləndirmə **§2–§6-nın tamamıdır**:

```
təkrarla (optimizator dayanana qədər):
   1. u ∈ [0,1]ⁿ          ← Nelder-Mead / Powell / DE növbəti nöqtəni verir
   2. u → fiziki parametr (hədlər arasında miqyaslama)
   3. MODEL DƏYİŞDİRİLİR   (k çarpanı, φ, Swc, nw, …)
   4. TAM SİMULYASİYA      ← §3 → §4/§5 → §6
   5. UYĞUNSUZLUQ:
         model nəticəsi müşahidə VAXTLARINA interpolyasiya olunur
         SSE = Σ w(hesab − ölçü)²,  RMSE = √(SSE/n),  NRMSE = RMSE/diapazon
         yekun = NRMSE-lərin çəkili ortası
   6. yığılmadısa → cərimə 1e6 (istisna ATILMIR, axtarış davam edir)
   7. nəticə keşlənir (eyni nöqtə təkrar sınana bilər)
```

Həssaslıq analizi eyni dövrəni işlədir, sadəcə nöqtələri optimizator yox,
**plan** seçir: Tornado — hər parametr öz hədləri arasında bir-bir (OAT);
elastiklik — baza nöqtəsi ətrafında ±kiçik addım, mərkəzi fərq.

---

## 9 · Bir cümlə ilə: nə nədən asılıdır

```
variogram        ← QC-dən sonrakı ÇEVRİLMİŞ data        (əvvəl çevir, sonra fit et)
kriging          ← variogram + anizotropluq + qonşuluq  (üçü də ondan ƏVVƏL)
həqiqi həndəsə   ← struktur səthlər                     (qalan xassələrdən ƏVVƏL)
T_ij             ← həndəsə + K                          (bir dəfə, zamandan asılı deyil)
WI               ← həndəsə + K + rw + skin              (bir dəfə)
mobillik λ       ← kr(S) + μ(p)                         (HƏR iterasiyada)
axın             ← T + λ_upstream + ΔΦ                  (HƏR iterasiyada)
ΔΦ               ← p + ρ(B)·g·D + Pc(S)                 (HƏR iterasiyada)
Jakobian         ← qalığın BÜTÜN yuxarıdakı üzvlərinin analitik törəmələri
Δt               ← Nyuton iterasiyalarının sayı (implicit) / CFL (IMPES)
BHP (THP rejimi) ← ƏVVƏLKİ addımın debitləri                (açıq birləşmə)
dəyişən keçid    ← YIĞILMIŞ vəziyyət                        (iterasiya daxilində YOX)
```

---

## 10 · Açıq qalanlar (⏳)

- ⏳ MPFA-O seçiləndə §4.2-nin 5-ci addımı çoxnöqtəli stensilə keçir,
  lakin bu sənəddə lokal sistemin öz daxili ardıcıllığı (interaction
  region → lokal həll → qlobal yığım) açılmır — `docs/mpfa_o_phase5a.md`.
- ⏳ `coupled_newton.py` / `standard_well.py` (OPM tipli birləşmiş quyu
  naməlumları) bu sənədin ardıcıllığına DAXİL DEYİL — onlar ayrı bir
  yoldur və defolt axında işləmir.
