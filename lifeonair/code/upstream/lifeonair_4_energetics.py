#!/usr/bin/env python3
"""lifeonair step 4 -- ENERGETIC FEASIBILITY + GLOBAL-SINK ESTIMATE.

Reproduces every number in Supplementary Note 1 of the manuscript:
  (1) Gibbs free energy of H2 and CO oxidation at ambient hyperarid mixing ratios (Nernst);
  (2) diffusive supply into the soil (Fick's first law, Millington-Quirk effective diffusivity);
  (3) community power supply vs per-cell basal maintenance, across the cell-density envelope;
  (4) biome-scale contribution to the global soil H2/CO sinks, and the CO2-equivalent.

Method adapted from Ortiz et al. 2021 PNAS (per-cell power balance); global denominators from
Ehhalt & Rohrer 2009 (H2) and Cordero et al. 2019 (CO); H2 GWP100 from Ocko & Hamburg 2022.

EVERY parameter below is flagged [measured] or [lit] (assumed from literature). The dominant
uncertainty is the absolute cell density, which is [lit] and should be replaced by a measured
Empty-Quarter value (16S rRNA qPCR copies/g) when one exists. This is an energetic FEASIBILITY
model, not a measured flux.

  python3 scripts/lifeonair_4_energetics.py
"""
import math

# ---- parameters -------------------------------------------------------------------------
T        = 303.15          # K, active soil temperature                                [lit]
R        = 8.3145e-3       # kJ/mol/K                                              [constant]
P_H2     = 0.53e-6         # atm, ambient H2                                            [lit]
P_CO     = 0.10e-6         # atm, ambient CO                                            [lit]
P_O2     = 0.2095          # atm                                                        [lit]
P_CO2    = 420e-6          # atm                                                        [lit]
DG0_H2   = -237.2          # kJ/mol, H2 + 0.5 O2 -> H2O                                 [lit]
DG0_CO   = -257.2          # kJ/mol, CO + 0.5 O2 -> CO2                                 [lit]
D0_H2    = 6.1e-5          # m2/s, free-air diffusion of H2 (light gas; NOT 2.4e-5)     [lit]
D0_CO    = 2.0e-5          # m2/s                                                       [lit]
POROSITY = 0.38            # dry desert sand                                            [lit]
DEPTH    = 0.10            # m, active horizon                                      [assumed]
BULK     = 1.6e6           # g/m3 dry bulk density (1.6 g/cm3)                          [lit]
GUILD_FR = 0.184           # guild share of community (step 3 v2, MAG-direct)      [measured]
CELLS_G  = (1e6, 1e7, 1e8) # cells/g dry soil envelope                                  [lit]
AREA     = 6.5e11          # m2, Rub' al-Khali (650,000 km2)                       [measured]
H2_RATE  = (5, 25)         # mg/m2/yr hyperarid H2 uptake (downscaled from Buzzard 2022)[lit]
CO_RATE  = (20, 100)       # mg/m2/yr                                                   [lit]
SINK_H2  = (60, 80)        # Tg/yr global soil H2 sink (Ehhalt & Rohrer 2009)           [lit]
SINK_CO  = 250             # Tg/yr global soil CO sink (Cordero et al. 2019)            [lit]
GWP100   = 11              # indirect GWP100 of H2 (Ocko & Hamburg 2022)                [lit]

# ---- 1. free energy at ambient mixing ratios (Nernst) -----------------------------------
RT = R * T
dG_H2 = DG0_H2 + RT * math.log(1.0 / (P_H2 * P_O2 ** 0.5))
dG_CO = DG0_CO + RT * math.log(P_CO2 / (P_CO * P_O2 ** 0.5))
print("== 1. free energy at ambient mixing ratios ==")
print(f"  dG(H2) = {dG_H2:7.1f} kJ/mol      dG(CO) = {dG_CO:7.1f} kJ/mol   (both strongly exergonic)")

# ---- 2. diffusive supply (Fick + Millington-Quirk) --------------------------------------
mol_m3 = 101325.0 / (8.314 * T)           # mol/m3 per atm at T
f_eff  = POROSITY ** (4.0 / 3.0)          # De/D0
J_H2 = D0_H2 * f_eff * (P_H2 * mol_m3) / DEPTH     # mol/m2/s
J_CO = D0_CO * f_eff * (P_CO * mol_m3) / DEPTH
print("== 2. diffusive supply (max potential, C=0 at depth) ==")
print(f"  De/D0 = {f_eff:.3f}   J(H2) = {J_H2*1e9:.2f} nmol/m2/s   J(CO) = {J_CO*1e9:.2f} nmol/m2/s")

# ---- 3. power balance vs maintenance ----------------------------------------------------
P_area = J_H2 * abs(dG_H2) * 1e3 + J_CO * abs(dG_CO) * 1e3     # W/m2
soil_m2 = BULK * DEPTH                                          # g soil per m2 column
print("== 3. community power balance ==")
print(f"  supply = {P_area:.2e} W/m2   ({soil_m2:.3g} g soil per m2 column)")
for d in CELLS_G:
    cells_m2 = soil_m2 * d * GUILD_FR
    print(f"  {d:.0e} cells/g -> guild {cells_m2:.2e} cells/m2 -> {P_area/cells_m2:.2e} W/cell")
print("  thresholds: basal survival ~1e-18..1e-19 W/cell; active maintenance ~1e-13..1e-15 W/cell")
print("  => supply exceeds the survival requirement by 2-4 orders across the envelope")
print("     (a maintenance/persistence budget, NOT growth)")

# ---- 4. contribution to global sinks ----------------------------------------------------
print("== 4. biome contribution to global sinks ==")
for lo, hi, sink, gas in ((H2_RATE[0], H2_RATE[1], SINK_H2, "H2"), (CO_RATE[0], CO_RATE[1], (SINK_CO, SINK_CO), "CO")):
    tg_lo, tg_hi = lo * AREA / 1e15, hi * AREA / 1e15           # mg/m2/yr * m2 -> Tg/yr
    pc_lo, pc_hi = 100 * tg_lo / sink[1], 100 * tg_hi / sink[0]
    print(f"  {gas}: {tg_lo:.4f}-{tg_hi:.4f} Tg/yr = {pc_lo:.3f}-{pc_hi:.3f}% of the global soil sink")
    if gas == "H2":
        print(f"      CO2-eq (GWP100={GWP100}): {tg_lo*GWP100*1e3:.0f}-{tg_hi*GWP100*1e3:.0f} kt CO2-eq/yr")
print("  => the sink contribution is minor (<1%); report it as one honest sentence, not the headline.")
print("     The result is FEASIBILITY (section 3); the novelty is the dataset's scale and biome.")
