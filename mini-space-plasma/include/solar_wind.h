/**
 * solar_wind.h — Solar Wind and Heliospheric Physics
 *
 * References:
 *   Parker (1958) ApJ 128:664 "Dynamics of the Interplanetary Gas"
 *   Kivelson & Russell §6 (Solar Wind)
 *   Hundhausen "Coronal Expansion and Solar Wind" (1972)
 *   MIT 22.611 Lectures 5-8
 *
 * Knowledge: L1 (definitions), L2 (solar wind concepts),
 *   L4 (Parker equation), L6 (Parker spiral, solar wind at 1 AU)
 */

#ifndef SOLAR_WIND_H
#define SOLAR_WIND_H

#include "space_plasma.h"

/*──────────────────────────────────────────────────────────────────────
 * L1: Parker Solar Wind Model Parameters
 *──────────────────────────────────────────────────────────────────────*/

/** @brief Parker wind solution type */
typedef enum {
    PARKER_BREEZE_I,      /* Subsonic everywhere (class I)       */
    PARKER_BREEZE_II,     /* Subsonic→supersonic→subsonic (II)  */
    PARKER_WIND_III,      /* Subsonic→supersonic (solar wind!)  */
    PARKER_ACCRETION_IV,  /* Supersonic→subsonic (accretion)     */
    PARKER_SUPERSONIC_V    /* Supersonic everywhere (class V)     */
} parker_solution_type_t;

/** @brief Parker solar wind solution at a radial point */
typedef struct {
    double r;              /* Radial distance [m]                 */
    double v;              /* Radial velocity [m/s]              */
    double n;              /* Number density [m⁻³]               */
    double T;              /* Temperature [K]                    */
    double Mach;           /* Mach number v/c_s                  */
    double c_s;            /* Sound speed [m/s]                  */
    parker_solution_type_t type;
} parker_profile_t;

/*──────────────────────────────────────────────────────────────────────
 * L1/L4: Parker Solar Wind Equation
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Parker's isothermal solar wind critical radius
 *
 * r_c = G M_sun / (2 c_s²)
 *
 * At r = r_c, the transonic solution has v = c_s.
 *
 * @param c_s   Sound speed [m/s]
 * @return      Critical radius [m]
 */
double parker_critical_radius(double c_s);

/**
 * @brief Solve Parker's isothermal wind equation for v(r)
 *
 * (v² - c_s²)(1/v) dv/dr = 2c_s²/r - GM_sun/r²
 *
 * Integrated form:
 *   (v/c_s)² - ln(v/c_s)² = 4 ln(r/r_c) + 4 r_c/r + C
 *
 * Uses the Lambert W function for the transonic solution.
 *
 * @param r       Radial distance [m] (must be > 0)
 * @param r_c     Critical radius [m]
 * @param v_guess Initial guess for velocity [m/s]
 * @param c_s     Isothermal sound speed [m/s]
 * @param tol     Convergence tolerance
 * @param max_iter Maximum iterations
 * @return        Solar wind velocity [m/s], or -1 on failure
 */
double parker_isothermal_velocity(double r, double r_c, double v_guess,
                                  double c_s, double tol, int max_iter);

/**
 * @brief Classify Parker solution at radius r
 *
 * Determines which of the 5 solution classes the point belongs to.
 *
 * @param r    Radial distance [m]
 * @param v    Velocity [m/s]
 * @param r_c  Critical radius [m]
 * @param c_s  Sound speed [m/s]
 * @return     Solution type enum
 */
parker_solution_type_t parker_classify(double r, double v, double r_c, double c_s);

/**
 * @brief Full Parker transonic wind profile
 *
 * Computes v(r), n(r), T(r) along the transonic (class III) solution.
 * Uses continuity: n(r) v(r) r² = n_0 v_0 r_0²
 *
 * @param r       Array of radial distances [N]
 * @param N       Number of points
 * @param r_c     Critical radius [m]
 * @param c_s     Sound speed [m/s]
 * @param r0      Reference radius [m]
 * @param n0      Density at reference radius [m⁻³]
 * @param profile Output: profile array [N]
 * @return        0 on success, -1 if solution fails
 */
int parker_transonic_profile(const double *r, size_t N, double r_c, double c_s,
                             double r0, double n0, parker_profile_t *profile);

/*──────────────────────────────────────────────────────────────────────
 * L2: Solar Wind Properties at 1 AU
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Typical solar wind parameters at 1 AU
 *
 * Fills a solar_wind_t with average quiet-sun values:
 *   v_sw ≈ 400 km/s, n_sw ≈ 5 cm⁻³, T_sw ≈ 1e5 K
 *   B_imf ≈ 5 nT
 *
 * @param sw    Output: filled solar_wind_t
 */
void solar_wind_quiet_sun(solar_wind_t *sw);

/**
 * @brief Fast solar wind (coronal hole) parameters at 1 AU
 *   v_sw ≈ 750 km/s, n_sw ≈ 3 cm⁻³, T_sw ≈ 8e5 K
 */
void solar_wind_fast(solar_wind_t *sw);

/**
 * @brief CME-like solar wind parameters at 1 AU
 *   v_sw up to 2000 km/s, n_sw elevated, strong B
 */
void solar_wind_cme(solar_wind_t *sw);

/*──────────────────────────────────────────────────────────────────────
 * L6: Parker Spiral
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Parker spiral magnetic field
 *
 * In spherical coordinates (r, θ, φ):
 *   B_r(r, θ, φ) = B_0 (r_0/r)²
 *   B_φ(r, θ, φ) = -B_r * (ω_sun * r * sin θ) / v_sw
 *   B_θ = 0
 *
 * φ is heliographic latitude (equator: θ = π/2).
 * ω_sun = 2π / (25.38 * 86400) ≈ 2.865e-6 rad/s (sidereal).
 *
 * @param r        Heliocentric distance [m]
 * @param theta    Colatitude [rad] (equator = π/2)
 * @param B0       Reference field at r0 [T]
 * @param r0       Reference radius [m]
 * @param omega    Solar rotation rate [rad/s]
 * @param v_sw     Solar wind speed [m/s]
 * @param B_out    Output: (B_r, B_theta, B_phi) in spherical [T]
 */
void parker_spiral_field(double r, double theta, double B0, double r0,
                         double omega, double v_sw, double B_out[3]);

/**
 * @brief Parker spiral angle
 *
 * The angle between the radial direction and the IMF:
 *   tan(ψ) = B_φ/B_r = -ω_sun r sin(θ) / v_sw
 *
 * @param r        Heliocentric distance [m]
 * @param theta    Colatitude [rad]
 * @param omega    Solar rotation rate [rad/s]
 * @param v_sw     Solar wind speed [m/s]
 * @return         Spiral angle [rad], in [-π/2, π/2]
 */
double parker_spiral_angle(double r, double theta, double omega, double v_sw);

/**
 * @brief IMF polarity from Parker spiral
 *
 * Returns +1 for away sector (positive B_r), -1 for toward sector.
 *
 * @param B_imf   IMF vector in RTN coordinates
 * @return        +1 (away) or -1 (toward)
 */
int imf_sector_polarity(const double B_imf[3]);

/*──────────────────────────────────────────────────────────────────────
 * L4: Solar Wind Number Density from Continuity
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Solar wind density from mass continuity
 *
 * n(r) = n_0 (r_0/r)² * (v_0 / v(r))
 *
 * For constant speed: n(r) ∝ r⁻² (spherical expansion).
 *
 * @param r     Radial distance [m]
 * @param r0    Reference radial distance [m]
 * @param n0    Reference density [m⁻³]
 * @param v0    Reference velocity [m/s]
 * @param v     Velocity at r [m/s]
 * @return      Density at r [m⁻³]
 */
double solar_wind_density(double r, double r0, double n0, double v0, double v);

/*──────────────────────────────────────────────────────────────────────
 * L2: Solar Wind Mach Numbers
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Alfven Mach number
 * M_A = v_sw / v_A
 */
double sp_alfven_mach(double v_sw, double B, double rho);

/**
 * @brief Sonic Mach number
 * M_S = v_sw / c_s
 */
double sp_sonic_mach(double v_sw, double c_s);

/**
 * @brief Magnetosonic Mach number
 * M_ms = v_sw / v_ms
 */
double sp_magnetosonic_mach(double v_sw, double v_A, double c_s);

/**
 * @brief Fast magnetosonic Mach number for arbitrary angle
 *
 * M_f = v_sw / v_f where v_f depends on propagation angle.
 *
 * @param v_sw   Solar wind speed [m/s]
 * @param v_A    Alfven speed [m/s]
 * @param c_s    Sound speed [m/s]
 * @param theta  Angle between k and B [rad]
 * @return       Fast Mach number
 */
double sp_fast_mach(double v_sw, double v_A, double c_s, double theta);

/*──────────────────────────────────────────────────────────────────────
 * L6: Co-rotating Interaction Regions (CIR)
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Check if fast wind catches slow wind (CIR condition)
 *
 * A CIR forms when fast wind from a coronal hole overtakes
 * slow wind emitted earlier from the same longitude.
 *
 * @param v_fast  Fast wind speed [m/s]
 * @param v_slow  Slow wind speed [m/s]
 * @param r_fast  Release radius of fast wind [m]
 * @param r_slow  Release radius of slow wind [m]
 * @return        r at which CIR forms [m], or INFINITY if no interaction
 */
double cir_formation_radius(double v_fast, double v_slow, double r_fast, double r_slow);

/**
 * @brief Compute solar wind travel time from Sun to distance r
 *
 * t = ∫_{r0}^{r} dr / v(r)
 * For constant v: t = (r - r0) / v
 *
 * @param r    Target distance [m]
 * @param r0   Starting distance [m]
 * @param v    Constant solar wind speed [m/s]
 * @return     Travel time [s]
 */
double solar_wind_travel_time(double r, double r0, double v);

/**
 * @brief Solar rotation (Carrington) period
 * @return 25.38 days in seconds
 */
double carrington_period(void);

/**
 * @brief Heliospheric current sheet tilt angle at distance r
 *
 * Approximate model: tilt increases with solar activity.
 *
 * @param r        Distance [m]
 * @param tilt_0   Tilt at source surface [rad]
 * @return         Tilt angle [rad]
 */
double hcs_tilt_angle(double r, double tilt_0);

#endif /* SOLAR_WIND_H */
