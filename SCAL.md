# Nisbi keçiricilik cədvəlləri (B4 + B7)

## Niyə lazımdır

Corey düsturu analitik və hamardır:

    krw = krw_end · Sn^nw

Real kern məlumatı belə davranmır: əyrilər asimmetrik olur, son
nöqtələr ayrıca ölçülür, bəzən orta hissədə əyilmə görünür. Bir
düsturla ifadə edilmir.

Üstəlik bir yataqda litologiya dəyişir — qumdaşı, əhəngdaşı, gilli
zona — və hər birinin öz əyrisi var.

## Texniki borcun bağlanması

`CoreyRelativePermeabilityAdapter` refaktorinqdə **müvəqqəti körpü**
kimi yazılmışdı (öz sənədləşməsində belə deyilir). Bu modul onun əsl
əvəzidir.

Adapter **silinmədi** — sadə modellər üçün hələ də faydalıdır və
bütün mövcud testlər ondan asılıdır. İkisi eyni interfeysi paylaşır;
model hansını işlədəcəyini özü seçir:

    model.scal_tables varsa   ->  TableRelativePermeabilityProvider
    yoxdursa                  ->  CoreyRelativePermeabilityAdapter

Doğrulama: eyni əyri iki yolla verildikdə RF **53.9045 %** və
**53.9049 %** — fərq yalnız cədvəlin diskretləşməsindən gəlir.

## Regionlar

`SATNUM` massivi hər hüceyrəni bir SCAL regionuna bağlayır. GRDECL
faylından oxunur (bax `ECLIPSE_IO.md`).

Ölçmə — eyni `Sw = 0.55`-də:

| Region | krw |
|---|---|
| Qum (Swc 0.15, krw_end 0.5) | 0.190 |
| Gilli (Swc 0.30, krw_end 0.2) | 0.062 |

İki regionlu modeldə RF 44.27 % — tək regionlu 53.90 %-dən aşağı,
çünki gilli zonada su daha çətin hərəkət edir.

### Doyumluluq hədləri — ən dar interval

Məhərrik doyumluluğu **bir** hədd cütü ilə kəsir. Regionlarda `Swc`
və `Sor` fərqli olduqda ən məhdudlaşdırıcı interval götürülür:

    Region 1: 0.15 … 0.80
    Region 2: 0.30 … 0.65
    Nəticə:   0.30 … 0.65

Əks halda bəzi hüceyrələrdə `Sw` cədvəldən kənara çıxardı.

### CFL limiti — ən sərt region

`max |dfw/dSw|` bütün regionlar üzrə **ən böyük** dəyəri qaytarır.
Zaman addımı modelin ən sərt hüceyrəsi ilə məhdudlaşmalıdır, orta ilə
yox — əks halda gilli zonada həll qeyri-stabil olardı.

## Fayl formatları

### CSV

    region,sw,krw,kro,pc
    1,0.20,0.000,0.800,0.50
    1,0.50,0.100,0.300,0.10
    2,0.30,0.000,0.750,0.80

`region` sütunu olmasa, hamısı 1-ci regiona düşür. `pc` istəyə görədir.

### Eclipse SWOF

    SWOF
    -- Sw    krw    kro    Pc
      0.20  0.000  0.800  0.5
      0.80  0.350  0.000  0.0
    /
      0.30  0.000  0.750  0.8
      0.65  0.200  0.000  0.0
    /

Hər `/` yeni regionu bitirir — Eclipse-in öz qaydası. Sətirlər
avtomatik `Sw` üzrə sıralanır.

Proqramın öz ixrac etdiyi deck geri oxuna bilir — test bunu qoruyur.

## Monotonluq yoxlaması

`krw` artan, `kro` azalan olmalıdır. Pozulsa diskretizasiya
qeyri-stabil olur, ona görə yoxlama **sərtdir** və fayl oxunarkən
dərhal xəta verilir:

    Region 1: krw azalır — monoton artan olmalıdır.

Həmçinin yoxlanılır: `Sw` artan sıralı, `[0, 1]` intervalında, `kr`
mənfi deyil və 1-dən böyük deyil, hərəkətli interval boş deyil.

## Törəmələr

Cədvəlin **dəqiq** parçalı-xətti meyli işlədilir, `np.gradient`-in
verdiyi hamar qiymət yox. Səbəb `A6_PLAN.md`-də ölçülüb: hamar
törəmə Nyutonun iterasiya sayını dəyişmir, lakin Jakobianı sonlu
fərqdən 10⁶ dəfə uzaqlaşdırır.

## Kapilyar təzyiq

Cədvəldə `Pc` sütunu varsa, o, Brooks-Corey analitik modelinin yerinə
işlədilir. Sütun yoxdursa və ya hamısı sıfırdırsa, kapilyar təzyiq
söndürülür.

## İnterfeys

`4 · NİSBİ KEÇİRİCİLİK` bölməsində mənbə seçimi:

| Rejim | Davranış |
|---|---|
| Corey düsturu | Aşağıdakı parametrlər işlədilir |
| Laboratoriya cədvəli | `CSV yüklə…` və ya `Eclipse SWOF yüklə…` |

Yüklənən cədvəlin xülasəsi göstərilir: region sayı, hər regionun
`Swc`, `Sor`, son nöqtələri və `Pc` sütununun olub-olmaması.
