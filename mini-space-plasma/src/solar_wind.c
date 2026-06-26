/**
 * solar_wind.c -- Solar Wind and Heliospheric Physics Implementation
 *
 * Every function implements one independent solar wind concept.
 * References:
 *   Parker (1958) ApJ 128:664
 *   Kivelson & Russell S6
 *   Hundhausen "Coronal Expansion and Solar Wind" (1972)
 *   MIT 22.611 Lectures 5-8
 *
 * Knowledge: L1-L6 solar wind physics
 */
#include "../include/space_plasma.h"
#include "../include/solar_wind.h"
#include <math.h>
#include <string.h>
#include <float.h>

/*========================================================================
 * L1/L4: Parker Solar Wind Equation
 *========================================================================*/

double parker_critical_radius(double c_s) {
    /* Parker's isothermal solar wind critical radius.
     * r_c = G*M_sun / (2*c_s^2)
     * At r=r_c the transonic solution has v=c_s.
     * For corona T~1e6 K, c_s~165 km/s, r_c~5.5 R_sun.
     * Reference: Parker (1958) ApJ 128:664 Eq. (12) */
    if (c_s <= 0.0) return INFINITY;
    return SP_G * SP_MSUN / (2.0 * c_s * c_s);
}

double parker_isothermal_velocity(double r, double r_c, double v_guess,
                                  double c_s, double tol, int max_iter) {
    /* Solve Parker's isothermal wind equation for v(r):
     * (v^2 - c_s^2) * dv/dr = 2*c_s^2*v/r - G*M_sun*v/r^2
     *
     * Integrated transcendental form:
     * (v/c_s)^2 - ln(v/c_s)^2 = 4*ln(r/r_c) + 4*r_c/r + C
     *
     * We solve by Newton-Raphson on the integrated equation.
     * Define F(M) = M^2 - ln(M^2) - 4*ln(r/r_c) - 4*r_c/r + 3
     * where M = v/c_s.
     *
     * The constant C=3 corresponds to the transonic solution (class III)
     * that passes through the critical point (r=r_c, M=1).
     *
     * Reference: Parker (1958) ApJ 128:664 Eq. (15)-(19)
     *            Kivelson & Russell S6.2 */
    if (r <= 0.0 || r_c <= 0.0 || c_s <= 0.0 || max_iter <= 0) return -1.0;
    if (tol <= 0.0) tol = 1e-8;

    double M = (v_guess > 0.0) ? v_guess / c_s : 0.5;
    /* For r > r_c, expect supersonic (M>1). For r < r_c, subsonic (M<1). */
    if (r > r_c && M < 1.0) M = 2.0;
    if (r < r_c && M > 1.0) M = 0.5;

    double target = 4.0 * log(r / r_c) + 4.0 * r_c / r - 3.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double M_sq = M * M;
        if (M_sq < 1e-30) { M = 0.01; continue; }

        double F = M_sq - log(M_sq) - target;
        double dF_dM = 2.0 * M - 2.0/M;  /* derivative: 2M - 2/M */

        if (fabs(dF_dM) < 1e-30) return -1.0;

        double dM = -F / dF_dM;

        /* Under-relax near the critical point */
        if (fabs(M - 1.0) < 0.01) dM *= 0.3;

        M += dM;
        if (M <= 0.0) M = 0.01;

        if (fabs(dM) < tol) break;
    }

    return M * c_s;
}

parker_solution_type_t parker_classify(double r, double v, double r_c, double c_s) {
    /* Classify Parker solution at radius r into 5 classes:
     * Class I:   subsonic everywhere (solar breeze type I)
     * Class II:  subsonic->supersonic->subsonic (unphysical, double-valued)
     * Class III: subsonic->supersonic (transonic, the observed solar wind)
     * Class IV:  supersonic->subsonic (accretion/bondi solution)
     * Class V:   supersonic everywhere
     *
     * Classification based on (v/c_s) compared to critical solution at r.
     * Reference: Parker (1958) Fig. 2, Kivelson & Russell Fig. 6.5 */
    if (r <= 0.0 || c_s <= 0.0) return PARKER_BREEZE_I;

    double M = v / c_s;
    double M_crit;

    if (r >= r_c) {
        /* Supersonic branch: r > r_c, the transonic solution has M > 1 */
        M_crit = 2.0;  /* approximate supersonic branch value */
    } else {
        /* Subsonic branch: r < r_c, M < 1 */
        M_crit = 0.5;
    }

    /* Solve for the exact transonic Mach number at this r */
    double v_trans = parker_isothermal_velocity(r, r_c, M_crit*c_s, c_s, 1e-8, 50);
    if (v_trans < 0.0) v_trans = M_crit * c_s;
    double M_trans = v_trans / c_s;

    if (r >= r_c) {
        /* r >= r_c: transonic solution is supersonic (M > 1) */
        if (M < 1.0) return PARKER_BREEZE_I;
        if (fabs(M - M_trans) < 0.1) return PARKER_WIND_III;
        if (M > M_trans) return PARKER_SUPERSONIC_V;
        return PARKER_BREEZE_II;
    } else {
        /* r < r_c: transonic solution is subsonic (M < 1) */
        if (M > 1.0) return PARKER_SUPERSONIC_V;
        if (fabs(M - M_trans) < 0.1) return PARKER_WIND_III;
        if (M > M_trans) return PARKER_BREEZE_II;
        return PARKER_BREEZE_I;
    }
}

int parker_transonic_profile(const double *r, size_t N, double r_c, double c_s,
                             double r0, double n0, parker_profile_t *profile) {
    /* Compute full Parker transonic wind profile v(r), n(r), T(r).
     * Mass continuity: n(r)*v(r)*r^2 = n0*v0*r0^2
     * Isothermal: T = constant = c_s^2*m_p/(gamma*kB)
     * Reference: Parker (1958), Kivelson & Russell S6.2 */
    if (!r || !profile || N == 0 || r_c <= 0.0 || c_s <= 0.0) return -1;

    /* Get velocity at reference point r0 */
    double v0_guess = (r0 > r_c) ? 2.0*c_s : 0.5*c_s;
    double v0 = parker_isothermal_velocity(r0, r_c, v0_guess, c_s, 1e-8, 100);
    if (v0 < 0.0) return -1;

    for (size_t i = 0; i < N; i++) {
        double ri = r[i];
        if (ri <= 0.0) return -1;

        double v_guess = (ri > r_c) ? 2.0*c_s : 0.5*c_s;
        double vi = parker_isothermal_velocity(ri, r_c, v_guess, c_s, 1e-8, 100);
        if (vi < 0.0) vi = v0 * r0 / ri;  /* fallback */

        profile[i].r    = ri;
        profile[i].v    = vi;
        profile[i].Mach = vi / c_s;
        profile[i].c_s  = c_s;

        /* Density from continuity (steady, spherical): n*v*r^2 = const */
        double mass_flux_const = n0 * v0 * r0 * r0;
        profile[i].n = mass_flux_const / (vi * ri * ri);

        /* Temperature from isothermal assumption */
        double m_proton = SP_MP;
        double gamma_eff = 5.0/3.0;
        profile[i].T = c_s * c_s * m_proton / (gamma_eff * SP_KB);

        profile[i].type = parker_classify(ri, vi, r_c, c_s);
    }

    return 0;
}

/*========================================================================
 * L2: Solar Wind Properties at 1 AU
 *========================================================================*/

void solar_wind_quiet_sun(solar_wind_t *sw) {
    /* Average quiet-Sun solar wind at 1 AU.
     * v ~ 400 km/s, n ~ 5 cm^-3, T ~ 1e5 K, B ~ 5 nT
     * Reference: Hundhausen (1972), Kivelson & Russell S6.3 */
    if (!sw) return;
    sw->v_sw = 400.0e3;
    sw->n_sw = 5.0e6;
    sw->T_sw = 1.0e5;
    sw->B_imf[0] = -3.0e-9;
    sw->B_imf[1] = 2.0e-9;
    sw->B_imf[2] = -1.0e-9;
    double rho = sw->n_sw * SP_MP;
    double vA = SP_MAG3(sw->B_imf) / sqrt(SP_MU0 * rho);
    double cs = sqrt(1.6666 * SP_KB * sw->T_sw / SP_MP);
    sw->Ma = sw->v_sw / vA;
    sw->Ms = sw->v_sw / cs;
    sw->M_ms = sw->v_sw / sqrt(vA*vA + cs*cs);
    sw->beta_sw = 2.0 * SP_MU0 * sw->n_sw * SP_KB * sw->T_sw
                  / (SP_MAG3(sw->B_imf) * SP_MAG3(sw->B_imf));
}

void solar_wind_fast(solar_wind_t *sw) {
    /* Fast solar wind from coronal holes at 1 AU.
     * v ~ 750 km/s, n ~ 3 cm^-3, T ~ 8e5 K
     * Reference: McComas et al. (1998) Geophys. Res. Lett. 25:1 */
    if (!sw) return;
    sw->v_sw = 750.0e3;
    sw->n_sw = 3.0e6;
    sw->T_sw = 8.0e5;
    sw->B_imf[0] = -2.0e-9;
    sw->B_imf[1] = 3.0e-9;
    sw->B_imf[2] = 0.5e-9;
    double rho = sw->n_sw * SP_MP;
    double vA = SP_MAG3(sw->B_imf) / sqrt(SP_MU0 * rho);
    double cs = sqrt(1.6666 * SP_KB * sw->T_sw / SP_MP);
    sw->Ma = sw->v_sw / vA;
    sw->Ms = sw->v_sw / cs;
    sw->M_ms = sw->v_sw / sqrt(vA*vA + cs*cs);
    sw->beta_sw = 2.0 * SP_MU0 * sw->n_sw * SP_KB * sw->T_sw
                  / (SP_MAG3(sw->B_imf) * SP_MAG3(sw->B_imf));
}

void solar_wind_cme(solar_wind_t *sw) {
    /* CME-like (Coronal Mass Ejection) parameters at 1 AU.
     * v up to 2000 km/s, n elevated, strong B, low beta.
     * Magnetic cloud structure with rotating B field.
     * Reference: Burlaga et al. (1981) J. Geophys. Res. 86:6673 */
    if (!sw) return;
    sw->v_sw = 1200.0e3;
    sw->n_sw = 20.0e6;
    sw->T_sw = 2.0e5;
    sw->B_imf[0] = 0.0;
    sw->B_imf[1] = 15.0e-9;
    sw->B_imf[2] = 15.0e-9;
    double rho = sw->n_sw * SP_MP;
    double Bmag = SP_MAG3(sw->B_imf);
    double vA = Bmag / sqrt(SP_MU0 * rho);
    double cs = sqrt(1.6666 * SP_KB * sw->T_sw / SP_MP);
    sw->Ma = sw->v_sw / vA;
    sw->Ms = sw->v_sw / cs;
    sw->M_ms = sw->v_sw / sqrt(vA*vA + cs*cs);
    sw->beta_sw = 2.0 * SP_MU0 * sw->n_sw * SP_KB * sw->T_sw / (Bmag * Bmag);
}

/*========================================================================
 * L6: Parker Spiral
 *========================================================================*/

void parker_spiral_field(double r, double theta, double B0, double r0,
                         double omega, double v_sw, double B_out[3]) {
    /* Parker spiral magnetic field in heliospheric coordinates.
     * B_r(r)   = B0 * (r0/r)^2
     * B_phi(r) = -B_r * omega*r*sin(theta)/v_sw
     * B_theta = 0
     *
     * The spiral angle psi satisfies tan(psi) = B_phi/B_r = -omega*r*sin(theta)/v_sw.
     * At 1 AU with v_sw=400 km/s: psi ~ 45 degrees.
     *
     * Reference: Parker (1958) ApJ 128:664 Eq. (50)-(52)
     *            Kivelson & Russell S6.4 */
    if (!B_out) return;
    if (r <= 0.0 || r0 <= 0.0 || v_sw <= 0.0) {
        B_out[0] = B_out[1] = B_out[2] = 0.0;
        return;
    }

    double Br = B0 * (r0 / r) * (r0 / r);
    double Bphi = -Br * omega * r * sin(theta) / v_sw;

    B_out[0] = Br;     /* radial */
    B_out[1] = 0.0;    /* theta */
    B_out[2] = Bphi;   /* phi */
}

double parker_spiral_angle(double r, double theta, double omega, double v_sw) {
    /* Parker spiral angle psi = atan(B_phi/B_r)
     * tan(psi) = -omega * r * sin(theta) / v_sw
     * At Earth (1 AU, equator): psi ~ 45 deg for slow wind, ~30 deg for fast.
     * Reference: Parker (1958), Kivelson & Russell S6.4 */
    if (r <= 0.0 || v_sw <= 0.0) return 0.0;
    double tan_psi = -omega * r * sin(theta) / v_sw;
    return atan(tan_psi);
}

int imf_sector_polarity(const double B_imf[3]) {
    /* IMF sector polarity: +1 = away from Sun (positive B_r),
     * -1 = toward Sun (negative B_r).
     * Sectors map to heliospheric current sheet crossings.
     * Reference: Wilcox & Ness (1965) J. Geophys. Res. 70:5793 */
    if (!B_imf) return 0;
    return (B_imf[0] >= 0.0) ? 1 : -1;
}

/*========================================================================
 * L4: Solar Wind Continuity
 *========================================================================*/

double solar_wind_density(double r, double r0, double n0, double v0, double v) {
    /* Mass continuity: n(r)*v(r)*r^2 = n0*v0*r0^2
     * For steady, spherically symmetric expansion.
     * Reference: Parker (1958), Hundhausen (1972) SII */
    if (r <= 0.0 || r0 <= 0.0 || v <= 0.0 || v0 <= 0.0) return 0.0;
    return n0 * (r0/r)*(r0/r) * (v0/v);
}

/*========================================================================
 * L2: Solar Wind Mach Numbers
 *========================================================================*/

double sp_alfven_mach(double v_sw, double B, double rho) {
    /* Alfven Mach number: M_A = v_sw / v_A
     * M_A > 1: super-Alfvenic (solar wind at 1 AU typically M_A ~ 8-10)
     * Reference: Kivelson & Russell S4.5 */
    if (rho <= 0.0) return INFINITY;
    double v_A = B / sqrt(SP_MU0 * rho);
    if (v_A <= 0.0) return INFINITY;
    return v_sw / v_A;
}

double sp_sonic_mach(double v_sw, double c_s) {
    /* Sonic Mach number: M_S = v_sw / c_s
     * M_S > 1: supersonic. Solar wind at 1 AU has M_S ~ 5-8.
     * Reference: Parker (1958) */
    if (c_s <= 0.0) return INFINITY;
    return v_sw / c_s;
}

double sp_magnetosonic_mach(double v_sw, double v_A, double c_s) {
    /* Magnetosonic Mach number: M_ms = v_sw / sqrt(v_A^2 + c_s^2)
     * Used for bow shock characterization.
     * Reference: Spreiter et al. (1966) Planet. Space Sci. 14:223 */
    double v_ms = sqrt(v_A*v_A + c_s*c_s);
    if (v_ms <= 0.0) return INFINITY;
    return v_sw / v_ms;
}

double sp_fast_mach(double v_sw, double v_A, double c_s, double theta) {
    /* Fast magnetosonic Mach number for arbitrary propagation angle.
     * v_f^2 solves v^4 - (vA^2+cs^2)v^2 + vA^2*cs^2*cos^2(theta) = 0
     * Reference: Kivelson & Russell S4.5 */
    if (v_sw <= 0.0) return 0.0;
    double vA2 = v_A * v_A;
    double cs2 = c_s * c_s;
    double sum_sq = vA2 + cs2;
    double disc = sum_sq*sum_sq - 4.0 * vA2 * cs2 * cos(theta)*cos(theta);
    if (disc < 0.0) disc = 0.0;
    double vf2 = 0.5 * (sum_sq + sqrt(disc));
    double v_f = sqrt(vf2);
    if (v_f <= 0.0) return INFINITY;
    return v_sw / v_f;
}

/*========================================================================
 * L6: Co-rotating Interaction Regions (CIR)
 *========================================================================*/

double cir_formation_radius(double v_fast, double v_slow, double r_fast, double r_slow) {
    /* CIR formation radius: where fast wind overtakes slow wind.
     * Fast wind emitted at r_fast from longitude phi at time t-r_fast/v_fast
     * catches slow wind emitted earlier at r_slow from same phi.
     * Formation when: r/v_fast + (r-r_fast)/v_fast = r/v_slow + (r-r_slow)/v_slow
     * -> r_cir = (v_fast*v_slow*(r_slow/v_slow - r_fast/v_fast)) / (v_fast - v_slow)
     *
     * Reference: Gosling & Pizzo (1999) Space Sci. Rev. 89:21
     *            Belcher & Davis (1971) J. Geophys. Res. 76:3534 */
    if (v_fast <= v_slow || v_slow <= 0.0) return INFINITY;

    double dt = (r_slow/v_slow) - (r_fast/v_fast);
    if (dt <= 0.0) return r_fast;

    double r_cir = v_fast * v_slow * dt / (v_fast - v_slow);
    return (r_cir > 0.0) ? r_cir : INFINITY;
}

double solar_wind_travel_time(double r, double r0, double v) {
    /* Solar wind travel time from r0 to r at constant speed v.
     * For typical speeds, Sun-to-Earth travel time is 2-5 days.
     * CME travel times: 1-3 days (fast), 3-5 days (slow).
     * Reference: Kivelson & Russell S6.5 */
    if (v <= 0.0 || r <= r0) return INFINITY;
    return (r - r0) / v;
}

double carrington_period(void) {
    /* Carrington rotation period (sidereal): 25.38 days.
     * Synodic period (as seen from Earth): 27.2753 days.
     * Reference: Carrington (1863) "Observations of the Spots on the Sun" */
    return 25.38 * 24.0 * 3600.0;  /* seconds */
}

double hcs_tilt_angle(double r, double tilt_0) {
    /* Heliospheric current sheet tilt angle.
     * Dipole tilt warps the current sheet into a "ballerina skirt" shape.
     * At solar minimum: tilt ~ 5-10 deg; at solar max: tilt ~ 50-70 deg.
     * Ballistic propagation: phi = phi_0 + omega*(r-r0)/v_sw creates spiral.
     *
     * Simplified model: tilt increases slowly with distance due to
     * decreasing magnetic pressure relative to solar wind ram pressure.
     * Reference: Smith (2001) J. Geophys. Res. 106:15819
     *            Hoeksema "The Heliospheric Current Sheet" (1995) */
    if (r <= 0.0) return tilt_0;

    /* Simple model: tilt grows with ln(r/r0) due to field-line curvature */
    double r0 = SP_RSUN * 2.5;  /* source surface at ~2.5 Rsun */
    if (r <= r0) return tilt_0;

    double growth = 0.1 * log(r / r0);
    return tilt_0 * (1.0 + growth);
}
