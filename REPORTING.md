# PDF hesabat (B6)

## Yeni asılılıq yoxdur

matplotlib onsuz da tələb olunur; `PdfPages` onun öz bir hissəsidir.
Məvcud renderer-lər (`rendering/renderers.py`) birbaşa işlədilir —
hesabatdakı qrafiklər interfeysdəki qrafiklərlə **eyni koddan** çıxır,
ikinci "hesabat üçün" versiyası yoxdur.

## Səhifələr

Hər bölmə **istəyə görədir** — model tək başına da, tam nəticə
dəstiylə də hesabat yaradıla bilər:

| Bölmə | Şərt |
|---|---|
| Başlıq | həmişə |
| Model xülasəsi | grid, süxur, flüid, SCAL, PVT, quyular, **faylar** |
| Diaqnostika | yalnız xəta/xəbərdarlıq varsa |
| Areal xəritələr | məsaməlilik + keçiricilik |
| SCAL əyriləri | Corey əsaslıdırsa (cədvəl əsaslı SCAL bu bölmədə deyil) |
| PVT əyriləri | cədvəl varsa |
| Nəticələr | simulyasiya nəticəsi verilibsə |
| Tarixçə uyğunluğu | müşahidə uyğunsuzluğu hesablanıbsa |

`ReportSections` ilə hər bölmə ayrıca söndürülə bilər.

## Tapılan və düzəldilən problem — tünd tema ağ səhifədə

Renderer-lər interfeysin **tünd** teması üçün yazılıb (`PALETTE.text`
kimi açıq rənglər tünd fon üçün nəzərdə tutulub). Onları birbaşa ağ
PDF səhifəsinə köçürəndə başlıqlar demək olar görünməz olurdu.

Həll: hər renderer-i palitra qəbul edəcək şəkildə dəyişmək əvəzinə
(unudulma riski yaradardı), çəkildikdən **sonra** bütün oxlar məcburi
işıqlandırılır (`_print_friendly()`) — fon ağ, mətn tünd, grid xətləri
açıq boz. Bu, hər renderer-in daxili detallarından asılı olmur.

## Tapılan ikinci problem — konfiqurasiya işləmirdi

İlk versiyada `ReportGenerator(sections=...)` konstruktoru
konfiqurasiyanı saxlayırdı, amma `write()` yalnız
`context.sections`-a baxırdı — ikisi arasında əlaqə yox idi.
`ReportSections(...)` konstruktora versəniz, heç bir təsiri olmurdu.

Düzəliş: `ReportGenerator` tam **vəziyyətsiz** (stateless) edildi —
konstruktor parametri yoxdur, bölmə seçimi yalnız
`ReportContext.sections` üzərindədir. Test bunu qoruyur: bütün
bölmələr söndürüləndə hesabat yalnız 1 səhifə (başlıq) olmalıdır.

## Metadata

PDF-in özündə başlıq, müəllif, mövzu və yaradılma tarixi yazılır —
`pypdf` ilə oxunub yoxlanılır (test).

## İnterfeys

`Layihə → PDF hesabat yaz…`. Cari model, son simulyasiya nəticəsi və
(varsa) tarixçə uyğunluğu avtomatik daxil edilir.
