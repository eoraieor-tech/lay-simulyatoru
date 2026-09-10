"""
Lay Simulyatoru — 3D 3-Fazalı Hidrodinamik Simulyasiya Platforması.

Modul xəritəsi:
    geostats   — kəşfiyyat quyularından 3D interpolyasiya (RBF / Kriging)
    grid       — 3D structured və corner-point grid, həndəsə, ACTNUM
    numerical  — MPFA-O transmissibilite mühərriki
    pvt        — PVT cədvəlləri və interpolyasiya (Bo, Bg, Rs, µ)
    relperm    — nisbi keçiricilik (2-faza + Stone II 3-faza)
    wells      — Peaceman quyu indeksi, BHP/THP, VFP, rejim nəzarəti
    engine     — IMPES / FIM həlledicisi
    analytics  — STOIIP, RF, GOR, WaterCut, kumulyativ hasilat
    viz        — PyVista 3D render, slice plane, dashboard

Ətraflı: ARCHITECTURE.md
"""

__version__ = "0.1.0.dev0"
