# Bu layihədə iş qaydaları

Bu fayl Claude Code (və digər AI köməkçiləri) üçün təlimatdır.

## Əsas prinsip

**Layihədə baş verən HƏR ŞEY bu repoda qeyd olunur.** Sahibkarın tələbi budur.
Heç bir qərar, dəyişiklik və ya tapıntı yalnız söhbətdə qalmamalıdır.

## Hər iş seansında

1. **Əvvəl oxu:** `ISH_HESABATI.md` (son bölmə) və `ROADMAP.md` — harada
   qaldığımızı bil.
2. **İş gör.**
3. **Sonra yaz:**
   - `ISH_HESABATI.md` → yeni bölmə: tarix, nə edildi, hansı qərar verildi,
     nə buraxıldı, hansı sual açıq qaldı
   - `QARARLAR.md` → texniki seçim edilibsə, səbəbi ilə
   - `ROADMAP.md` → mərhələ statusu yenilənir
   - `ARCHITECTURE.md` → struktur dəyişibsə
4. **Commit + push:** hər mənalı addımdan sonra, işi yarımçıq saxlama.

## Sənəd qaydaları

- Dil: **Azərbaycan dili**. Kod və identifikatorlar ingiliscə.
- Doldurulmamış yerlər ⏳ işarəsi ilə açıq göstərilir — uydurma məzmun YAZILMIR.
- Köhnə hesabat qeydləri **silinmir və dəyişdirilmir**, yalnız yenisi əlavə olunur.
- Səhv qərar verilibsə: `QARARLAR.md`-da köhnəsi "Əvəzlənib" işarələnir,
  yenisi altına yazılır — tarixçə itmir.

## Commit mesajları

Qısa, məzmunlu, Azərbaycanca və ya ingiliscə:
```
faza 1.2: grid strukturu əlavə edildi
docs: arxitektura sənədi yeniləndi
fix: interpolyasiya sərhəd xətası
```

## Nə etmə

- Tələb aydın deyilsə, uydurma — soruş.
- Sənədləri yeniləmədən push etmə.
- Böyük dəyişikliyi bir commit-ə yığma.

---

## graphify (bilik qrafiki)

Bu layihədə `graphify-out/` altında bilik qrafiki var — god node-lar,
icma strukturu, fayllararası əlaqələr.

Qaydalar:

- Kod bazası ilə bağlı sual olanda əvvəlcə `graphify query "<sual>"`
  işlət (`graphify-out/graph.json` varsa). Əlaqələr üçün
  `graphify path "<A>" "<B>"`, anlayışlar üçün `graphify explain "<anlayış>"`.
  Bunlar dar altqraf qaytarır — `GRAPH_REPORT.md` və ya çılpaq grep-dən
  xeyli kiçikdir.
- `graphify-out/wiki/index.md` varsa, geniş naviqasiya üçün onu işlət.
- `graphify-out/GRAPH_REPORT.md`-i yalnız geniş arxitektura baxışı üçün oxu.
- Kod dəyişdikdən sonra `graphify update .` işlət (AST-only, API xərci yoxdur).
