/**
 * magnetosphere.h — Planetary Magnetospheres
 *
 * References:
 *   Kivelson & Russell §7-§10 (Earth's Magnetosphere)
 *   Chapman & Ferraro (1931) "Terrestrial Magnetic Storms"
 *   Dungey (1961) Phys. Rev. Lett. 6:47 "Interplanetary Magnetic Field"
 *   MIT 22.611 Lectures 9-14
 *
 * Knowledge: L1 (definitions), L2 (magnetosphere concepts),
 *   L4 (Chapman-Ferraro, dipole), L6 (Earth's magnetosphere),
 *   L7 (space weather)
 */

#ifndef MAGNETOSPHERE_H
#define MAGNETOSPHERE_H

#include "space_plasma.h"

/*──────────────────────────────────────────────────────────────────────
 * L1: Dipole Field
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Earth's magnetic dipole field
 *
 * In spherical coordinates (r, θ, φ) with dipole moment M aligned with -z:
 *   B_r = -(2 μ₀ M cos θ) / (4π r³)
 *   B_θ = -(μ₀ M sin θ) / (4π r³)
 *   B_φ = 0
 *
 * Or in Cartesian (dipole along z), for position vector x:
 *   B = (μ₀/(4πr⁵)) * [3(m·x̂)x̂ - m̂] * M
 *
 * @param x     Position vector in Earth radii [x, y, z] (units of R_E)
 * @param M     Dipole moment [A·m²]
 * @param B_out Output: B field [T] in Cartesian
 */
void earth_dipole_field(const double x[3], double M, double B_out[3]);

/**
 * @brief Dipole field line equation
 *
 * r = L R_E sin²(θ)
 *
 * Where L is the McIlwain L-parameter (equatorial crossing distance in R_E).
 *
 * @param L      L-shell parameter
 * @param theta  Colatitude [rad]
 * @return       Radial distance [m]
 */
double dipole_fieldline_r(double L, double theta);

/**
 * @brief Dipole field magnitude along field line
 * B = (B_E / L³) * sqrt(1 + 3 cos²(θ)) / sin⁶(θ)
 *
 * @param L       L-shell
 * @param theta   Colatitude [rad]
 * @return        |B| [T]
 */
double dipole_field_magnitude(double L, double theta);

/**
 * @brief Compute McIlwain L-parameter from (r, θ)
 * L = r / (R_E sin²(θ))
 */
double mcilwain_L(double r, double theta);

/**
 * @brief Compute invariant latitude from L-shell
 * Λ = acos(1/√L)
 */
double invariant_latitude(double L);

/*──────────────────────────────────────────────────────────────────────
 * L4: Chapman-Ferraro Magnetopause
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Chapman-Ferraro magnetopause standoff distance
 *
 * Balance between solar wind dynamic pressure and magnetic pressure:
 *   p_sw = ½ ρ_sw v_sw²  (ram pressure)
 *   p_B  = B²/(2μ₀) = (B_E²/(2μ₀)) * (R_E/R_mp)⁶  (dipole at nose)
 *
 * Setting p_sw = f * p_B (f ≈ 2 for spherical, ≈ 0.88 for realistic):
 *   R_mp/R_E = (B_E² / (2 μ₀ f ρ_sw v_sw²))^(1/6)
 *
 * @param n_sw   Solar wind density [m⁻³]
 * @param v_sw   Solar wind speed [m/s]
 * @param f      Compression factor (≈ 2.44 for gasdynamic, ~1.0 empirical)
 * @return       Magnetopause standoff distance [R_E]
 */
double chapman_ferraro_standoff(double n_sw, double v_sw, double f);

/**
 * @brief Chapman-Ferraro magnetopause current
 *
 * The magnetopause current (Chapman-Ferraro current) flows on the
 * dayside boundary and cancels the geomagnetic field outside.
 *
 * Current density: J_CF ≈ 2 B_dipole / (μ₀ δ) where δ ≈ d_i (ion skin depth)
 *
 * @param B_dipole  Dipole field at magnetopause [T]
 * @param delta     Boundary thickness [m]
 * @return          Surface current density [A/m]
 */
double chapman_ferraro_current(double B_dipole, double delta);

/*──────────────────────────────────────────────────────────────────────
 * L2: Magnetospheric Convection
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Magnetospheric E×B convection velocity
 *
 * v_E = (E × B) / B²
 *
 * For uniform electric field E (dawn-dusk) and dipole B:
 * produces sunward convection in inner magnetosphere.
 *
 * @param E      Electric field [V/m]
 * @param B      Magnetic field [T]
 * @param v_out  Output: E×B drift velocity [m/s]
 */
void magnetosphere_exb_drift(const double E[3], const double B[3], double v_out[3]);

/**
 * @brief Cross-polar cap potential
 *
 * Φ_PC = E_sw * L_PC where E_sw ≈ v_sw B_z (southward IMF)
 * and L_PC ≈ 30 R_E (polar cap width).
 *
 * Empirical: Φ_PC [kV] ≈ 30 + 15 * v_sw[km/s] * B_s[nT] (Boyle et al. 1997)
 *
 * @param v_sw   Solar wind speed [km/s]
 * @param Bz_imf IMF Bz (negative = southward) [nT]
 * @return       Cross-polar cap potential [kV]
 */
double cross_polar_cap_potential(double v_sw, double Bz_imf);

/**
 * @brief Dungey cycle reconnection voltage
 *
 * Φ_Dungey ≈ v_sw * B_s * L_MP / 2
 *
 * @param v_sw   Solar wind speed [m/s]
 * @param B_s    Southward IMF component [T]
 * @param L_MP   Magnetopause width [m]
 * @return       Reconnection voltage [V]
 */
double dungey_cycle_voltage(double v_sw, double B_s, double L_MP);

/*──────────────────────────────────────────────────────────────────────
 * L2/L6: Magnetic Field Models
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Parabolic magnetopause shape
 *
 * R(θ) = R_mp * (2 / (1 + cos θ))       [dayside]
 * Cylindrical tail: √(y²+z²) = R_tail   [nightside]
 *
 * @param theta  Angle from Sun-Earth line [rad] (0 = subsolar)
 * @param R_mp   Standoff distance [R_E]
 * @param R_tail Tail radius [R_E]
 * @param x      Output: x-coordinate [R_E]
 * @param r      Output: radial distance from axis [R_E]
 */
void magnetopause_shape(double theta, double R_mp, double R_tail,
                        double *x, double *r);

/**
 * @brief Tsyganenko-like stretched field model
 *
 * Simple 2D model incorporating dipole + tail current sheet.
 *
 * @param x      GSM x [R_E] (sunward positive)
 * @param z      GSM z [R_E] (northward positive)
 * @param Dst    Disturbance storm-time index [nT]
 * @param B_out  Output: B field [T]
 */
void tsyganenko_simple(double x, double z, double Dst, double B_out[3]);

/*──────────────────────────────────────────────────────────────────────
 * L2: Ring Current
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Dessler-Parker-Sckopke relation
 *
 * Dst [nT] ≈ -3.98 × 10⁻³⁰ * U_R [J]
 *
 * Where U_R is the total ring current energy.
 *
 * @param total_energy  Ring current particle energy [J]
 * @return              Dst contribution [nT]
 */
double dessler_parker_sckopke(double total_energy);

/**
 * @brief Ring current energy density
 *
 * For a bi-Maxwellian distribution with density n, temperature T,
 * and pressure anisotropy A = T_perp/T_par - 1:
 *
 * U_R = ³⁄₂ n k_B T * (1 + ⅔ A)
 *
 * @param n       Density [m⁻³]
 * @param T_perp  Perpendicular temperature [K]
 * @param T_par   Parallel temperature [K]
 * @return        Energy density [J/m³]
 */
double ring_current_energy_density(double n, double T_perp, double T_par);

/*──────────────────────────────────────────────────────────────────────
 * L2: Plasmapause & Plasmasphere
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Plasmapause location (Carpenter & Anderson 1992)
 *
 * L_pp = 5.6 - 0.46 * Kp_max
 *
 * Where Kp_max is the maximum Kp in preceding 24 hours.
 * For average conditions (Kp ≈ 2): L_pp ≈ 4.7 R_E
 *
 * @param Kp_max  Maximum Kp index (0-9)
 * @return        Plasmapause L-shell
 */
double plasmapause_L(double Kp_max);

/**
 * @brief Plasmasphere density profile
 *
 * n_e(L) = n_0 * (L_pp / L)^α  for L < L_pp
 * n_e(L) = n_trough  for L > L_pp
 *
 * Typical: n_0 ≈ 1000 cm⁻³, α ≈ 3-4, n_trough ≈ 1 cm⁻³
 *
 * @param L       L-shell
 * @param L_pp    Plasmapause L-shell
 * @param n0      Plasmasphere base density [cm⁻³]
 * @param n_trough Trough density [cm⁻³]
 * @param alpha   Power-law exponent
 * @return        Density [cm⁻³]
 */
double plasmasphere_density(double L, double L_pp, double n0,
                            double n_trough, double alpha);

/*──────────────────────────────────────────────────────────────────────
 * L6: Bow Shock
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Bow shock standoff distance
 *
 * R_bs ≈ 1.3 * R_mp (empirical, gasdynamic)
 * More precisely (Farris & Russell 1994):
 *   R_bs/R_mp = 1 + 1.1 * (ρ_up/ρ_down - 1)⁻¹
 *
 * @param R_mp    Magnetopause standoff [R_E]
 * @param M_ms    Upstream magnetosonic Mach number
 * @param gamma   Polytropic index
 * @return        Bow shock standoff [R_E]
 */
double bow_shock_standoff(double R_mp, double M_ms, double gamma);

/**
 * @brief Rankine-Hugoniot jump conditions for perpendicular MHD shock
 *
 * For a shock with normal along x:
 *   ρ₂/ρ₁ = (γ+1)M₁² / ((γ-1)M₁² + 2)           [density jump]
 *   p₂/p₁ = (2γM₁² - (γ-1)) / (γ+1)              [pressure jump]
 *
 * @param M1      Upstream fast Mach number
 * @param gamma   Polytropic index
 * @param rho_ratio Output: ρ₂/ρ₁
 * @param p_ratio   Output: p₂/p₁
 * @param B_ratio   Output: B₂/B₁ (= ρ₂/ρ₁ for perp shock)
 */
void mhd_shock_jump(double M1, double gamma, double *rho_ratio,
                    double *p_ratio, double *B_ratio);

/**
 * @brief Magnetosheath flow deflection
 *
 * Flow is deflected around magnetopause obstacle.
 * Deflection angle ≈ atan(vy/vx) depends on position.
 *
 * @param x      Distance from nose along Sun-Earth [R_E]
 * @param y      Distance from axis [R_E]
 * @param R_mp   Magnetopause standoff [R_E]
 * @param vx_out Output: x-velocity (normalized to upstream)
 * @param vy_out Output: y-velocity (normalized to upstream)
 */
void magnetosheath_flow(double x, double y, double R_mp,
                        double *vx_out, double *vy_out);

/*──────────────────────────────────────────────────────────────────────
 * L6: Particle Drifts in Magnetosphere
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Gradient-B drift velocity
 *
 * v_∇B = (m v_perp² / (2qB³)) B × ∇B
 *
 * @param v_perp  Perpendicular velocity [m/s]
 * @param B       Magnetic field vector [T]
 * @param gradB   Gradient of |B| [T/m], 3 components
 * @param q       Particle charge [C] (signed)
 * @param m       Particle mass [kg]
 * @param v_drift Output: drift velocity [m/s]
 */
void gradB_drift(double v_perp, const double B[3], const double gradB[3],
                 double q, double m, double v_drift[3]);

/**
 * @brief Curvature drift velocity
 *
 * v_c = (m v_par² / (q B²)) B × (b̂·∇)b̂
 *
 * where b̂ = B/|B| and (b̂·∇)b̂ = curvature vector.
 *
 * @param v_par    Parallel velocity [m/s]
 * @param B        Magnetic field [T]
 * @param curv     Curvature vector [1/m], 3 components
 * @param q        Particle charge [C] (signed)
 * @param m        Particle mass [kg]
 * @param v_drift  Output: drift velocity [m/s]
 */
void curvature_drift(double v_par, const double B[3], const double curv[3],
                     double q, double m, double v_drift[3]);

/**
 * @brief Total guiding center drift (grad-B + curvature)
 *
 * For dipole field with isotropic pressure: v_d ∝ L²
 * Electrons drift eastward, ions westward.
 *
 * @param E_kin   Particle kinetic energy [J]
 * @param L       L-shell
 * @param q       Charge [C] (signed)
 * @param v_drift Output: azimuthal drift speed [m/s]
 */
void dipole_drift_speed(double E_kin, double L, double q, double v_drift[3]);

/**
 * @brief Bounce period for particles in dipole field
 *
 * τ_b ≈ (4 L R_E / v) * (1.30 - 0.56 sin(α_eq))
 *
 * Where α_eq is the equatorial pitch angle.
 *
 * @param E_kin     Particle kinetic energy [J]
 * @param m         Particle mass [kg]
 * @param L         L-shell
 * @param alpha_eq  Equatorial pitch angle [rad]
 * @return          Bounce period [s]
 */
double dipole_bounce_period(double E_kin, double m, double L, double alpha_eq);

#endif /* MAGNETOSPHERE_H */
