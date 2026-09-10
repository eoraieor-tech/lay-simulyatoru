# Lay Simulyatoru

> Sıfırdan qurulan lay (rezervuar) simulyasiya layihəsi.
> Bu repo layihənin **bütün detallarını** — kodu, qərarları, mərhələləri və iş
> jurnalını — bir yerdə saxlayır.

**Status:** 🟡 Başlanğıc mərhələsi — struktur quruldu, texniki tələblər toplanır.
**Başlama tarixi:** 10 sentyabr 2026

---

## 1. Layihə nədir

Lay Simulyatoru — neft-qaz yataqlarında lay (rezervuar) proseslərinin
riyazi modelləşdirilməsi və simulyasiyası üçün proqram təminatıdır.

**Məqsəd:** kəşfiyyat quyularının 3D koordinat və geoloji məlumatlarından
(məsaməlik, keçiricilik tenzoru) başlayaraq 3D heterogen lay modelini quran,
5-spot quyu şəbəkəsində su ilə sıxışdırma və qazın neftdən ayrılmasını (Rs)
günbəgün simulyasiya edən, RF / THP / BHP / lay təzyiqini hesablayan və
nəticələri interaktiv 3D-də göstərən platforma.

**Əsas texniki fərq:** klassik TPFA əvəzinə **MPFA-O** (Multi-Point Flux
Approximation, O-sxemi) — anizotrop keçiricilik tenzorlarında və
qeyri-ortoqonal / corner-point gridlərdə axın dəqiqliyini qorumaq üçün.

**Əhatə dairəsi (scope):**

| Var | Yoxdur |
|---|---|
| 3D corner-point grid | Kompozisiya (EOS) modeli |
| 3 faza: neft, su, qaz (black-oil) | Termal (buxar) proseslər |
| MPFA-O diskretizasiyası | Qeyri-struktur (PEBI) grid |
| Geostatistik interpolyasiya (Kriging, SGS) | Paralel / GPU hesablama |
| IMPES + tam implicit (FIM) həlledici | |
| Peaceman quyuları, BHP/THP, VFP | |
| 3D vizualizasiya və dashboard | |

> ⚠️ **10 sentyabr 2026:** audit göstərdi ki, bu işin böyük hissəsi
> sahibkarın əvvəlki `lay-simulyatoru` layihəsində (65 175 sətir,
> 1 993 keçən test) **artıq mövcuddur**. Layihə sıfırdan yazılmır.
> Bax: [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md) və
> [QARARLAR.md](QARARLAR.md) → Q-04.

---

## 2. Sənədlərin xəritəsi

Layihənin hər aspekti ayrıca sənəddə qeyd olunur. Bir şey axtarırsansa,
başlanğıc nöqtəsi budur:

| Sənəd | Nə var içində |
|---|---|
| [README.md](README.md) | Bu fayl — ümumi baxış və naviqasiya |
| [ARCHITECTURE.md](ARCHITECTURE.md) | Sistemin quruluşu, qatlar, modullar, məlumat axını |
| [ROADMAP.md](ROADMAP.md) | Mərhələlər (fazalar), nə bitib / nə qalıb |
| [ISH_HESABATI.md](ISH_HESABATI.md) | İş jurnalı — hər addım, hər dəyişiklik, tarixlə |
| [QARARLAR.md](QARARLAR.md) | Texniki qərarlar və **niyə** məhz belə seçildi |
| [AUDIT_2026-09-10.md](AUDIT_2026-09-10.md) | Mövcud `lay-simulyatoru` kod bazasının tam auditi |
| [berpa/](berpa/) | Git tarixçəsindən bərpa edilmiş kod (A7 qaz fazası) |
| [docs/](docs/) | Dərin mövzular: fizika, riyaziyyat, formatlar, alqoritmlər |

---

## 3. Quraşdırma

⏳ *Texnologiya seçimindən sonra doldurulacaq.*

---

## 4. İstifadə

⏳ *İlk işlək versiyadan sonra doldurulacaq.*

---

## 5. Bu repo necə işləyir

Qayda sadədir: **layihədə baş verən hər şey burada qeyd olunur.**

- Yeni tələb gəldi → [ROADMAP.md](ROADMAP.md)-ə mərhələ kimi yazılır
- Texniki seçim edildi → [QARARLAR.md](QARARLAR.md)-a səbəbi ilə yazılır
- İş görüldü → [ISH_HESABATI.md](ISH_HESABATI.md)-na tarixlə yazılır
- Struktur dəyişdi → [ARCHITECTURE.md](ARCHITECTURE.md) yenilənir
- Hər mərhələdən sonra → `git commit` + `git push`

Detallı qaydalar: [CLAUDE.md](CLAUDE.md)

---

## 6. Əlaqə

**Müəllif:** eoraieor-tech
**Repo:** https://github.com/eoraieor-tech/LAY-SIMULYATIR-MODELI-
