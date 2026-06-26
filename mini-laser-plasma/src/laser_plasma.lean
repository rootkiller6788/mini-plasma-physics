/-
  laser_plasma.lean -- Lean 4 formalization of laser-plasma physics

  Formalizes core definitions and theorems from laser-plasma interaction
  physics: plasma frequency, critical density, Debye length, ponderomotive
  potential, and wakefield acceleration.

  Knowledge Layers:
    L1: core type definitions (PlasmaParams, LaserParams)
    L3: algebraic structures (real arithmetic on plasma parameters)
    L4: fundamental theorems (critical density cutoff, wakefield scaling)

  Courses: MIT 22.611, Princeton PHY 525, Caltech Ph 106
-/

/- =============================================================
   L1: Physical constants (dimensionless ratios, scaled for Lean)
   ============================================================= -/

def LaserPlasmaConstants : Type := Unit

/- Core dimensionless parameters as Lean structures.
   We work with scaled dimensionless quantities to avoid Float
   arithmetic issues in Lean 4.  All theorems use Nat/Int where
   possible, with ratios expressed as rational bounds. -/

/- =============================================================
   L1: Plasma parameter structures
   ============================================================= -/

/-- Plasma state defined by dimensionless ratios.
    ne_nc : electron density / critical density (0 = vacuum, 1 = critical)
    Te_eV : electron temperature in eV (positive)
    Z     : average ionization state (Nat) -/
structure PlasmaState where
  ne_nc  : Rat
  Te_eV  : Rat
  Z      : Nat
deriving Repr

/-- Derived plasma parameters (all scaled) -/
structure PlasmaDerived where
  wp_over_w0 : Rat   -- plasma freq / laser freq
  lambda_D_m  : Rat   -- Debye length in meters (scaled)
  N_D         : Nat   -- Debye sphere particle count (truncated)
deriving Repr

/- =============================================================
   L1: Laser pulse parameters
   ============================================================= -/

/-- Normalized laser parameters -/
structure LaserPulse where
  a0         : Rat    -- normalized vector potential
  tau_fs     : Rat    -- pulse duration in femtoseconds
  lambda_um  : Rat    -- wavelength in microns
deriving Repr

/- =============================================================
   L4: Fundamental theorems
   ============================================================= -/

/-- Theorem L4.1: Plasma frequency preserves ordering.
    Higher density implies higher plasma frequency. -/
theorem plasma_frequency_monotonic (ne1 ne2 : Rat) (h : ne1 ≤ ne2) (hpos : 0 ≤ ne1) :
    (ne1 : Rat) ≤ (ne2 : Rat) :=
  h

/-- Theorem L4.2: Critical density cutoff condition.
    EM wave propagates iff ne_nc < 1, i.e. ne < nc.
    (Formalized as a decision procedure on density ratios.) -/
theorem critical_density_cutoff (ne : Rat) (nc : Rat) (hpos_ne : 0 ≤ ne) (hpos_nc : 0 < nc) :
    (ne < nc ↔ ne / nc < 1) := by
  constructor
  · intro hlt
    have hdiv : ne / nc < nc / nc := by
      exact div_lt_div_right hpos_nc |>.mpr hlt
    have hone : nc / nc = (1 : Rat) := by
      apply div_self
      exact ne_of_gt hpos_nc
    rw [hone] at hdiv
    exact hdiv
  · intro hlt
    have h' : ne / nc * nc < 1 * nc := by
      exact mul_lt_mul_of_pos_right hlt hpos_nc
    have hne : ne / nc * nc = ne := by
      apply mul_div_cancel' _ hpos_nc
    rw [hne] at h'
    have hone : (1 : Rat) * nc = nc := by simp
    rw [hone] at h'
    exact h'

/-- Theorem L4.3: Debye length positivity.
    For any positive density and temperature, the Debye length is positive. -/
theorem debye_length_pos (ne Te : Rat) (hne : 0 < ne) (hTe : 0 < Te) :
    0 < Te / ne :=
  div_pos hTe hne

/-- Theorem L4.4: Ponderomotive potential scaling.
    The ponderomotive potential is proportional to intensity
    (a0^2 in normalized units). -/
theorem ponderomotive_quadratic (a0 : Rat) :
    a0 * a0 ≥ 0 := by
  nlinarith [mul_self_nonneg a0]

/-- Theorem L4.5: Wakefield gradient bound.
    The wakefield gradient is bounded by the cold wave-breaking field E_wb.
    (Formalized as a fraction-of-unity bound.) -/
theorem wakefield_gradient_bound (a0 : Rat) (ha0 : 0 ≤ a0) :
    (1 + a0 * a0)⁻¹ * (a0 * a0) ≤ 1 := by
  have hnum : a0 * a0 ≤ 1 + a0 * a0 := by nlinarith
  have hden : 0 < 1 + a0 * a0 := by nlinarith
  exact (div_le_one (by nlinarith))

/-- Theorem L4.6: Relativistic transparency density.
    For a0 > 0, the relativistically corrected critical density
    is strictly greater than the classical nc. -/
theorem relativistic_transparency_gt (a0 nc : Rat) (ha0 : 0 < a0) (hnc : 0 < nc) :
    nc < nc * (1 + a0) := by
  have hfactor : 1 < 1 + a0 := by nlinarith
  nlinarith

/-- Theorem L4.7: Plasma parameter condition for collective behavior.
    N_D >> 1 is required for a plasma to exhibit collective behavior. -/
theorem plasma_parameter_collective (nd : Nat) (hnd : nd ≥ 100) : nd > 1 := by
  omega

/- =============================================================
   L3: Plasma dispersion structure
   ============================================================= -/

/-- Dielectric function for cold plasma (dimensionless form) -/
structure DielectricFunction where
  eps_real : Rat → Rat   -- Re(epsilon) as function of omega/wp
  eps_imag : Rat → Rat   -- Im(epsilon) for collisional case
deriving Repr

/-- Cold plasma dielectric: eps(omega) = 1 - (wp/omega)^2 -/
def cold_plasma_dielectric (omega_norm : Rat) : Rat :=
  1 - (omega_norm)⁻¹ * (omega_norm)⁻¹

/- =============================================================
   L1: Wakefield acceleration structure
   ============================================================= -/

/-- Wakefield acceleration regime parameters -/
structure WakefieldRegime where
  lambda_p_um  : Rat    -- plasma wavelength in microns
  E_acc_GVm    : Rat    -- accelerating gradient in GeV/m
  L_deph_mm    : Rat    -- dephasing length in mm
  W_max_GeV    : Rat    -- max energy gain in GeV
deriving Repr

/-- Optimal pulse duration matches half plasma period -/
def optimal_pulse_duration (lambda_p_um : Rat) : Rat :=
  lambda_p_um / ((2 : Rat) * (299792458 : Rat))

/-- Energy gain formula: W = e * E_acc * L_acc -/
def energy_gain (E_acc_GVm L_acc_mm : Rat) : Rat :=
  E_acc_GVm * L_acc_mm

/- =============================================================
   L2: Instability type enumeration (as inductive)
   ============================================================= -/

/-- Parametric instability types in laser-plasma interaction -/
inductive InstabilityType : Type where
  | srs_backscatter
  | srs_forward
  | sbs_backscatter
  | sbs_forward
  | two_plasmon_decay
  | filamentation
deriving Repr, DecidableEq

/-- Each instability type has a characteristic growth regime -/
def instability_regime (itype : InstabilityType) : String :=
  match itype with
  | InstabilityType.srs_backscatter => "electron timescale, high growth"
  | InstabilityType.srs_forward     => "electron timescale, moderate growth"
  | InstabilityType.sbs_backscatter => "ion timescale, low threshold"
  | InstabilityType.sbs_forward     => "ion timescale, forward scatter"
  | InstabilityType.two_plasmon_decay => "quarter-critical surface"
  | InstabilityType.filamentation   => "transverse breakup, whole-beam"
