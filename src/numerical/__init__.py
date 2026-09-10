"""numerical — MPFA-O axın diskretizasiyası.

Multi-Point Flux Approximation (O-scheme) transmissibilite mühərriki.
TPFA-dan fərqli olaraq üz axınını yalnız iki hüceyrə mərkəzindən deyil,
təpə nöqtəsi ətrafındakı interaction volume-un bütün iştirakçılarından
hesablayır — anizotrop tenzor və qeyri-ortoqonal gridlərdə dəqiqlik
buna görə qorunur.

Nəticə: scipy.sparse formatında transmissibilite matrisi.

Status: ⏳ Mərhələ 3-də yazılacaq
"""
