"""Physical, atomic, mineralogical, and body constants for Mnemosyne.

All quantities are CGS.  Each value carries an inline reference comment.

The mineral vapor-pressure database uses the Clausius-Clapeyron relation in
the exact form consumed by dishoom's ``functions4.Pvap_rocks``:
    P_vap = exp(-A / T + B)   [dyne cm⁻²],  T in K
where ``A`` [K], ``B`` (dimensionless), and ``mu`` [amu] are stored alongside
their primary literature source in :data:`MINERALS`.  This is the natural-exp,
cgs-pressure convention of van Lieshout+ 2014 (A&A 572, A76), NOT a base-10
log or a bar-referenced pressure — feeding base-10 / bar coefficients here
would corrupt every sublimation rate.
"""

from __future__ import annotations

from typing import NamedTuple

# ── Physical constants ───────────────────────────────────────────────────────
G     = 6.67384e-8         # gravitational constant [cm³ g⁻¹ s⁻²]  NIST CODATA 2014
k_B   = 1.3806488e-16      # Boltzmann constant [erg K⁻¹]           NIST CODATA 2014
amu   = 1.66053892e-24     # atomic mass unit [g]                    NIST CODATA 2014
au2cm = 1.49597870700e13   # 1 AU [cm]                               IAU 2012

# ── Atomic masses ────────────────────────────────────────────────────────────
m_Na = 22.9897 * amu   # sodium [g]     IUPAC atomic weights 2021
m_K  = 39.0983 * amu   # potassium [g]  IUPAC atomic weights 2021
m_Fe = 55.845  * amu   # iron [g]       IUPAC atomic weights 2021
m_Mg = 24.305  * amu   # magnesium [g]  IUPAC atomic weights 2021

# ── Photoionization rates at 1 AU ────────────────────────────────────────────
k_hv_Na_Gstar = 5.92e-6   # Na, G2-star [s⁻¹]  Huebner & Mukherjee 2015
k_hv_Na_Kstar = 9.5e-7    # Na, K2-star [s⁻¹]  Huang+ 2017 (MUSCLES spectrum)
k_hv_Na_Fstar = 1.3e-5    # Na, F/K0-star [s⁻¹] Oza (dishoom)
k_hv_Na_Astar = 0.0298    # Na, A-star [s⁻¹]   Oza (KELT-9, T_eff=10170 K)
k_hv_Na       = k_hv_Na_Gstar   # Na default (solar G-type)  Huebner & Mukherjee 2015

# ── Solar-system body parameters ─────────────────────────────────────────────
r_Io = 1.821e8   # Io radius [cm]  NASA fact sheet / Anderson+ 2001
m_Io = 8.93e25   # Io mass [g]     Anderson+ 2001

# ── Volcanic / compositional mass fractions ──────────────────────────────────
X_Na_volcanic = 0.015  # Io volcanic Na mass fraction  Lellouch+ 2003, 2015
#   Notes: ~0.3% of Io's atmosphere is NaCl (Lellouch 2003);
#   1/5 of atmosphere is volcanic, 100% of sublimation is SO2 (not NaCl).
#   Range in dishoom: Xna_0_total_low=0.003, Xna_0_total_high=0.013,
#   Xna_0_total_chondrite=0.13 (Fegley & Zolotov 2000 upper limit).


# ── Mineral vapor-pressure database ─────────────────────────────────────────
class MineralEntry(NamedTuple):
    """Single entry in the mineral vapor-pressure database.

    Coefficients satisfy: P_vap = exp(-A / T + B) [dyne cm⁻²], with T in K.

    Attributes:
        A: Clausius-Clapeyron slope [K]  (A = ΔH_sub / R for a 2-term fit).
        B: Dimensionless offset (sets the pressure normalisation in dyne cm⁻²).
        mu: Molar mass [amu].
        reference: Primary literature source for these coefficients.
    """
    A: float
    B: float
    mu: float
    reference: str


MINERALS: dict[str, MineralEntry] = {
    # ── Silicate / refractory minerals ────────────────────────────────────────
    # All from van Lieshout+ 2014, A&A 572, A76 (Table 3) — verbatim.
    # The paper's own relation IS P_v = exp(-A/T + B) in dyne cm⁻², so these
    # A, B drop straight into dishoom's Pvap_rocks with no unit conversion.
    # (Cross-check: their Fe row, A=48354/B=29.2, gives a 1-atm boiling point
    #  of ~3146 K, matching iron's measured 3134 K — convention confirmed.)
    "MgSiO3":  MineralEntry(68908, 38.1, 100.389, "van Lieshout+ 2014"),            # enstatite
    "Al2O3":   MineralEntry(77365, 39.3, 101.961, "van Lieshout+ 2014; Schaefer & Fegley 2004"),  # corundum
    "SiO":     MineralEntry(49520, 32.5,  44.085, "van Lieshout+ 2014"),
    "SiC":     MineralEntry(78462, 37.8,  40.10,  "van Lieshout+ 2014"),            # silicon carbide
    "SiO2":    MineralEntry(69444, 33.1,  60.084, "van Lieshout+ 2014"),            # quartz
    "C":       MineralEntry(93646, 36.7,  12.011, "van Lieshout+ 2014"),            # graphite
    "Mg2SiO4": MineralEntry(65308, 34.1, 140.694, "van Lieshout+ 2014"),            # forsterite (Io mantle 76–85%, Sohl 2002)
    "Fe2SiO4": MineralEntry(60377, 37.7, 203.774, "van Lieshout+ 2014"),            # fayalite

    # ── Pure metallic elements ────────────────────────────────────────────────
    # Fit to Alcock+ 1984 (Can. Metall. Q. 23:309) full 4-term equation
    #   log₁₀(P/atm) = a + b/T + c·log₁₀T + d/T³
    # over the high-T (liquid / near-boiling) regime, then collapsed to the
    # 2-term  exp(-A/T + B)  [dyne cm⁻²]  form dishoom expects.  For the alkali
    # metals and the liquid-phase transition metals Alcock gives c=d=0, so the
    # collapse is exact (≤0.1%); for Cr/Mn (only solid-phase coeffs published)
    # the 2-term fit holds to ~1–4%, well inside Alcock's stated ±5%.
    #
    # NOTE: the previous values here were wrong — they used a bogus 2-term
    # reduction that dropped the c·log₁₀T curvature term AND mixed up atm vs.
    # bar in the offset constant (used +13.82 = ln(10⁶) for bar while claiming
    # an atm-based source). Under dishoom's exp(-A/T+B), those shifted boiling
    # points by 60–270 K and skewed every metal sublimation rate.
    #
    # Fe is taken verbatim from van Lieshout+ 2014's own Table 3 row (iron is
    # tabulated there), keeping it consistent with the silicate block above.
    "Fe":      MineralEntry(48354, 29.2,   55.845, "van Lieshout+ 2014 (Table 3)"),  # T_boil ≈ 3146 K
    "Ni":      MineralEntry(47813, 29.178, 58.693, "fit to Alcock+ 1984 (liquid)"),  # T_boil ≈ 3115 K
    "Cr":      MineralEntry(48748, 33.328, 51.996, "fit to Alcock+ 1984 (solid)"),   # solid-phase fit; bp~2944 K extrapolated
    "Mn":      MineralEntry(31535, 28.087, 54.938, "fit to Alcock+ 1984 (solid)"),   # T_boil ≈ 2334 K
    "Co":      MineralEntry(47383, 28.768, 58.933, "fit to Alcock+ 1984 (liquid)"),  # T_boil ≈ 3172 K
    "Na":      MineralEntry(12381, 24.660, 22.990, "fit to Alcock+ 1984 (liquid)"),  # T_boil ≈ 1143 K
    "K":       MineralEntry(10253, 23.965, 39.098, "fit to Alcock+ 1984 (liquid)"),  # T_boil ≈ 1012 K

    # ── Alkali halides (Io-relevant: NaCl/KCl are the main surface Na/K reservoirs)
    # 2-term Clausius-Clapeyron from the sublimation enthalpy and normal
    # boiling point:  A = ΔH_sub / R,  B = ln(1 atm in dyne cm⁻²) + A / T_boil
    # so the curve passes through 1 atm = 1.01325×10⁶ dyne cm⁻² at T_boil.
    # ΔH_sub and T_boil from NIST-JANAF / CRC, consistent with the standard
    # alkali-halide reference Lamoreaux & Hildenbrand 1984, JPCRD 13(4).
    # (Earlier values used the wrong offset constant and a too-low T_boil.)
    "NaCl":    MineralEntry(27627, 29.724, 58.443, "NIST-JANAF / CRC; Lamoreaux & Hildenbrand 1984"),  # T_boil ≈ 1738 K
    "KCl":     MineralEntry(25979, 29.201, 74.551, "NIST-JANAF / CRC; Lamoreaux & Hildenbrand 1984"),  # T_boil ≈ 1690 K

    # ── Metal oxides ─────────────────────────────────────────────────────────
    # 2-term Clausius-Clapeyron fits anchored on the effective vaporization
    # enthalpy and vaporization/boiling temperature (A = ΔH/R, B set to 1 atm
    # in dyne cm⁻² at that T).  These oxides vaporize incongruently (partial
    # decomposition to suboxides + O₂), so treat them as order-of-magnitude
    # surface-sublimation estimates rather than precision curves.  ΔH/T from
    # NIST-JANAF (Chase 1998); CaSiO₃ slope from Sossi+ 2019 evaporation data.
    "NiO":     MineralEntry(55325, 38.638, 74.692, "NIST-JANAF (Chase 1998); decomposes ≈2230 K"),  # bunsenite
    "CaO":     MineralEntry(69758, 36.187, 56.077, "NIST-JANAF (Chase 1998)"),     # lime, T_vap ≈ 3120 K
    "CaSiO3":  MineralEntry(74600, 34.902, 116.162, "Sossi+ 2019 (incongruent)"),  # wollastonite
    "TiO2":    MineralEntry(67954, 34.770,  79.866, "NIST-JANAF (Chase 1998)"),    # rutile, T_vap ≈ 3245 K

    # ── Metal sulfides (chondritic surface phases; dominant Fe/Ni UV source candidates)
    # 2-term Clausius-Clapeyron estimates (A = ΔH_vap / R, B anchored to a
    # measured (T, P) point on the congruent-vaporization curve).  Both sulfides
    # vaporize INCONGRUENTLY near/above their melting points (FeS → Fe + ½S₂,
    # NiS → Ni + ½S₂), so these are effective surface-sublimation curves valid in
    # the sub-decomposition regime, not true molecular-FeS/NiS sublimation lines.
    #   FeS: ΔH_vap ≈ 315 kJ/mol, anchored at the 1463 K melting point where the
    #        congruent vapour pressure is ≈ 10⁻³ atm (JANAF; Ferrante+ 1971).
    #   NiS: ΔH_vap ≈ 330 kJ/mol, anchored at ≈ 1070 K / ≈ 10⁻⁴ atm near its
    #        decomposition point (JANAF; Rosenqvist 1954).
    "FeS":     MineralEntry(37886, 32.817,  87.910, "JANAF; Ferrante+ 1971 (congruent, incongruent above ~1460 K)"),  # troilite
    "NiS":     MineralEntry(39690, 41.712,  90.758, "JANAF; Rosenqvist 1954 (decomposes near ~1070 K)"),               # millerite
}
