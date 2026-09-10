"""engine — Simulyasiya həlledicisi.

3D 3-fazalı IMPES (Implicit Pressure, Explicit Saturation) sxemi,
sonrakı mərhələdə FIM (Fully Implicit).

Addım: MPFA-O matrisi ilə A·P = B təzyiq həlli → upwind sxemi ilə
Sw, So, Sg yenilənməsi → P < Psat olduqda qazın ayrılması (Rs nəzarəti)
→ CFL şərtinə görə dinamik zaman addımı (dt).

Status: ⏳ Mərhələ 4-də yazılacaq
"""
