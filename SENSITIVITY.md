# Həssaslıq analizi (C6)

Uyğunlaşdırmadan (C5) əvvəl soruşulmalı sual: **hansı parametr
nəticəyə nə qədər təsir edir?** Cavab olmadan optimallaşdırıcı bütün
parametrləri eyni ciddiyyətlə axtarır, halbuki bəziləri praktik
olaraq əhəmiyyətsiz ola bilər — bu, axtarış fəzasını lüzumsuz
genişləndirir və yığılmanı çətinləşdirir.

## İki üsul

| Üsul | Sual | Necə |
|---|---|---|
| **Tornado** | Hansı parametr öz TAM diapazonunda ən çox təsir edir? | Hər parametr öz hədləri arasında, digərləri baza dəyərində (one-at-a-time) |
| **Yerli elastiklik** | Hazırkı modeldə kiçik dəyişiklik nəyə təsir edir? | Baza nöqtəsi ətrafında kiçik ±addım |

Fərq vacibdir: geniş hədli, lakin lokal olaraq az təsirli parametr
Tornado-da yuxarıda, elastiklikdə aşağıda görünə bilər.

## Çıxış ölçüləri

RF, kumulyativ neft, son sulaşma, su gəlişi vaxtı, son orta təzyiq —
istənilənini seçmək olar.

`Su gəlişi vaxtı` üçün xüsusi qayda: su gəlişi baş verməyibsə,
simulyasiya müddəti qaytarılır. `None` qaytarmaq hesablamanı
sındırardı — "hələ gəlməyib" "sonsuz gecikmə" kimi işlədilir.

## Ölçülmüş nəticələr — fizika ilə uyğunluq

Beş parametrlə (`PERM_MULT`, `PORO_MULT`, `SOR`, `SWC`, `KRW_END`),
`RATE` idarə olunan vurucu, 1200 gün:

### Tornado — RF

| Parametr | Yayılma |
|---|---|
| `SOR` | **25.5** |
| `KRW_END` | 18.1 |
| `PORO_MULT` | 9.8 |
| `PERM_MULT` | 5.8 |
| `SWC` | 4.1 |

`SOR` ən böyük hərəkətvericidir — məntiqlidir, çünki o, RF-in
nəzəri yuxarı həddini birbaşa təyin edir.

### Tornado — su gəlişi vaxtı: gözlənilməyən nəticə

| Parametr | Yayılma (gün) |
|---|---|
| `SOR` | 580.7 |
| `PORO_MULT` | 560.2 |
| `KRW_END` | 550.4 |
| `SWC` | 509.5 |
| **`PERM_MULT`** | **10.7** |

`PERM_MULT`-un demək olar heç təsiri yoxdur (10 gün, digərləri
500+ gün). Bu, səhv deyil — **Buckley-Leverett nəzəriyyəsinin klassik
nəticəsidir**: `RATE` ilə idarə olunan vurucuda cəbhə sürəti vurulan
həcmdən (məsamə həcminə nisbətdə) asılıdır, mütləq keçiricilikdən
demək olar asılı deyil. Eyni sürətlə su vurulursa, keçiricilik 10 dəfə
dəyişsə də cəbhə oxşar zamanda çatır — yalnız təzyiq profili dəyişir.

Test bunu qoruyur: `PERM_MULT`-un yayılması `SOR`-un yayılmasının
1/10-dan azdır.

## Uğursuz hədd

Optimallaşdırıcı kimi, həssaslıq analizi də hədlərin kənarını sınayır
və oradakı model yığılmaya bilər (məsələn `SOR` maksimuma yaxın olanda
hərəkətli interval demək olar boşalır).

Belə hal **istisna atmır** — baza dəyəri ilə əvəzlənir və cədvəldə
"UĞURSUZ" işarələnir, tornado qrafikində isə həmin zolaq ştrixlənir.
Susmaq təhlükəlidir: istifadəçi bütün parametrlərin sınandığını
düşünərdi.

## Tornado qrafiki

Zolaqlar yayılmaya görə sıralanır (ən böyüyü yuxarıda — klassik
tərtibat). Rəng istiqaməti göstərir: parametr artanda çıxış artırsa
bir rəng, azaltsa digəri. Şaquli xətt baza dəyərini göstərir.

## İnterfeys

**"Həssaslıq"** tabı: üsul (Tornado / Yerli elastiklik), çıxış ölçüsü,
addım ölçüsü (yalnız yerli üsul üçün), başlat/dayandır, nəticə
cədvəli, tornado qrafiki. Axtarış fon axınında gedir.

## Vaxt xərci

Hər parametr üçün **iki** əlavə simulyasiya (aşağı və yuxarı hədd).
5 parametrlə 11 simulyasiya (1 baza + 10). Böyük modellərdə fully
implicit mühərrik burada da həlledicidir.
