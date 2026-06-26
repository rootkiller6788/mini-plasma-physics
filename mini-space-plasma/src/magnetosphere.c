/**
 * magnetosphere.c -- Planetary Magnetosphere Implementation
 *
 * Every function implements one independent magnetospheric concept.
 * References:
 *   Kivelson & Russell S7-S10
 *   Chapman & Ferraro (1931) Terr. Mag. Atmos. Electr. 36:77,171
 *   Dungey (1961) Phys. Rev. Lett. 6:47
 *   MIT 22.611 Lectures 9-14
 *
 * Knowledge: L1-L6 magnetospheric physics
 */
#include "../include/space_plasma.h"
#include "../include/magnetosphere.h"
#include <math.h>

/*========================================================================
 * L1: Earth Dipole Field
 *========================================================================*/

void earth_dipole_field(const double x[3], double M, double B_out[3]) {
    /* Earth's magnetic dipole field in Cartesian coordinates.
     * Dipole moment aligned with -z (geographic axis, south magnetic pole).
     *
     * B = (mu0/(4*pi*r^5)) * [3(m.r)r - m*r^2]
     *
     * where m = (0, 0, -M) is the dipole moment vector.
     * For Earth: M_E = 7.78e22 A.m^2, B_E(equator) ~ 3.12e-5 T.
     *
     * In spherical coords (r, theta measured from +z):
     *   B_r = -(mu0*M*cos(theta)) / (2*pi*r^3)
     *   B_theta = -(mu0*M*sin(theta)) / (4*pi*r^3)
     *   B_phi = 0
     *
     * Reference: Chapman & Bartels "Geomagnetism" (1940)
     *            Kivelson & Russell S7.2 */
    if (!x || !B_out) return;

    /* Convert from R_E units to meters */
    double xm[3];
    xm[0] = x[0] * SP_RE;
    xm[1] = x[1] * SP_RE;
    xm[2] = x[2] * SP_RE;

    double r2 = xm[0]*xm[0] + xm[1]*xm[1] + xm[2]*xm[2];
    double r = sqrt(r2);
    if (r < 1.0) r = 1.0;
    double r5 = r2 * r2 * r;

    /* Dipole moment: m = (0, 0, -M) */
    double m_dot_r = -M * xm[2];
    double coeff = SP_MU0 / (4.0 * M_PI * r5);

    B_out[0] = coeff * (3.0 * m_dot_r * xm[0]);
    B_out[1] = coeff * (3.0 * m_dot_r * xm[1]);
    B_out[2] = coeff * (3.0 * m_dot_r * xm[2] + M * r2);  /* +M*r2 from -(-M)*r2 */
}

double dipole_fieldline_r(double L, double theta) {
    /* Dipole field line equation: r = L * R_E * sin^2(theta)
     * L is the McIlwain L-parameter (equatorial crossing in R_E).
     * theta is colatitude (0 at north pole, pi/2 at equator).
     * Reference: McIlwain (1961) J. Geophys. Res. 66:3681 */
    if (L <= 0.0) return 0.0;
    return L * SP_RE * sin(theta) * sin(theta);
}

double dipole_field_magnitude(double L, double theta) {
    /* Dipole field magnitude along field line:
     * B = (B_E/L^3) * sqrt(1+3*cos^2(theta)) / sin^6(theta)
     * where B_E = mu0*M/(4*pi*R_E^3) ~ 3.12e-5 T.
     * Reference: Schulz & Lanzerotti "Particle Diffusion in Radiation Belts" (1974) */
    if (L <= 0.0) return 0.0;
    double B_E = SP_B0_EARTH;
    double st = sin(theta);
    if (st < 1e-10) return INFINITY;  /* avoid division by zero at poles */
    double ct = cos(theta);
    double numerator = sqrt(1.0 + 3.0 * ct * ct);
    double denominator = st * st * st * st * st * st;
    return (B_E / (L*L*L)) * numerator / denominator;
}

double mcilwain_L(double r, double theta) {
    /* McIlwain L-parameter from (r, theta):
     * L = r / (R_E * sin^2(theta))
     * For a dipole field line in a centered dipole.
     * Reference: McIlwain (1961) J. Geophys. Res. 66:3681 */
    if (theta <= 0.0 || sin(theta) < 1e-10) return INFINITY;
    return r / (SP_RE * sin(theta) * sin(theta));
}

double invariant_latitude(double L) {
    /* Invariant latitude Lambda from L-shell:
     * Lambda = acos(1/sqrt(L))
     * This is the footpoint latitude of the field line at Earth's surface.
     * Reference: O'Brien (1963) J. Geophys. Res. 68:989 */
    if (L < 1.0) return 0.0;
    return acos(1.0 / sqrt(L));
}

/*========================================================================
 * L4: Chapman-Ferraro Magnetopause
 *========================================================================*/

double chapman_ferraro_standoff(double n_sw, double v_sw, double f) {
    /* Chapman-Ferraro magnetopause standoff distance.
     *
     * Balance: solar wind ram pressure = magnetic pressure at nose
     *   p_sw = n_sw * m_p * v_sw^2 (dynamic pressure)
     *   p_B  = B_mp^2 / (2*mu0) = B_E^2*(R_E/R_mp)^6 / (2*mu0) (dipole compression)
     *
     * With compression factor f (typically 2.0-2.44 for spherical obstacle):
     *   R_mp/R_E = (B_E^2 / (2*mu0*f*p_sw))^(1/6)
     *
     * For typical solar wind at 1 AU: R_mp ~ 10 R_E.
     * During geomagnetic storms: R_mp can shrink to 6-7 R_E.
     *
     * Reference: Chapman & Ferraro (1931)
     *            Spreiter et al. (1966) Planet. Space Sci. 14:223 */
    if (n_sw <= 0.0 || v_sw <= 0.0 || f <= 0.0) return INFINITY;

    double p_sw = n_sw * SP_MP * v_sw * v_sw;  /* ram pressure */
    double B_E = SP_B0_EARTH;
    double p_B0 = B_E * B_E / (2.0 * SP_MU0);  /* magnetic pressure at 1 R_E */

    return pow(p_B0 / (f * p_sw), 1.0/6.0);
}

double chapman_ferraro_current(double B_dipole, double delta) {
    /* Chapman-Ferraro magnetopause current density.
     * The CF current flows on the magnetopause boundary and cancels
     * the geomagnetic field outside the magnetosphere.
     *
     * Ampere's law across boundary: J_s = delta_B / (mu0*delta)
     * where delta_B = 2*B_dipole (field doubles inside, zeroes outside).
     * J_s ~ 2*B_dipole/(mu0*delta), where delta ~ d_i (ion inertial length).
     *
     * For B_dipole ~ 5e-8 T at magnetopause, delta~100 km:
     *   J_s ~ 0.1 A/m (total current ~ 1 MA)
     *
     * Reference: Chapman & Ferraro (1931), Kivelson & Russell S7.4 */
    if (delta <= 0.0) return INFINITY;
    return 2.0 * B_dipole / (SP_MU0 * delta);
}

/*========================================================================
 * L2: Magnetospheric Convection
 *========================================================================*/

void magnetosphere_exb_drift(const double E[3], const double B[3], double v_out[3]) {
    /* E x B drift velocity: v_E = (E x B) / B^2
     *
     * For uniform dawn-dusk electric field E_y and dipole B:
     * produces sunward convection in inner magnetosphere
     * and anti-sunward return flow on flanks.
     *
     * This drives the Dungey cycle: dayside reconnection ->
     * anti-sunward lobe convection -> nightside reconnection ->
     * sunward inner magnetosphere return flow.
     *
     * Reference: Dungey (1961) Phys. Rev. Lett. 6:47
     *            Kivelson & Russell S9.3 */
    if (!E || !B || !v_out) return;

    /* E x B */
    SP_CROSS3(E, B, v_out);

    double B2 = SP_DOT3(B, B);
    if (B2 < 1e-30) {
        v_out[0] = v_out[1] = v_out[2] = 0.0;
        return;
    }

    v_out[0] /= B2;
    v_out[1] /= B2;
    v_out[2] /= B2;
}

double cross_polar_cap_potential(double v_sw, double Bz_imf) {
    /* Cross-polar cap potential (CPCP).
     *
     * Empirical formula (Boyle et al. 1997):
     *   Phi_PC [kV] = 30 + 15 * v_sw[km/s] * B_s[nT]
     *
     * where B_s = -Bz_imf for southward IMF (Bz < 0).
     * For northward IMF (Bz > 0), B_s = 0.
     *
     * CPCP is a key measure of solar wind-magnetosphere coupling.
     * Typical values: 30-50 kV (quiet), 100-150 kV (storm).
     *
     * Reference: Boyle et al. (1997) J. Geophys. Res. 102:111
     *            Reiff & Luhmann "Solar Wind-Magnetosphere Coupling" (1986) */
    if (v_sw < 0.0) v_sw = 0.0;
    double v_sw_kms = v_sw / 1000.0;
    double B_s = (Bz_imf < 0.0) ? -Bz_imf : 0.0;  /* Bz in nT */
    return 30.0 + 15.0 * v_sw_kms * B_s;
}

double dungey_cycle_voltage(double v_sw, double B_s, double L_MP) {
    /* Dungey cycle reconnection voltage.
     * Phi_Dungey = v_sw * B_s * L_MP_eff
     * where L_MP_eff is the width of the reconnection line.
     *
     * For southward IMF (B_s > 0): reconnection at dayside magnetopause.
     * Typical Phi_Dungey ~ 50-200 kV for average solar wind.
     *
     * Reference: Dungey (1961) Phys. Rev. Lett. 6:47
     *            Cowley "Magnetospheric Reconnection" (1985) */
    if (v_sw <= 0.0 || B_s <= 0.0 || L_MP <= 0.0) return 0.0;
    return v_sw * B_s * L_MP;
}

/*========================================================================
 * L2/L6: Magnetic Field Models
 *========================================================================*/

void magnetopause_shape(double theta, double R_mp, double R_tail,
                        double *x, double *r) {
    /* Parabolic magnetopause shape (dayside, x > 0):
     * x = R_mp * 2*cos(theta)/(1+sin^2(theta))
     * r = sqrt(y^2+z^2) = R * sin(theta)
     *
     * Simplified Shue et al. (1997) shape:
     *   r = R_mp * (2/(1+cos(theta)))^alpha,  x = r*cos(theta)
     *
     * For theta=0 (subsolar): x = R_mp, r = 0
     * For theta=pi/2 (terminator): x = 0, r = k*R_mp
     *
     * Nightside (x < 0): cylindrical approximation
     *   sqrt(y^2+z^2) = R_tail (constant)
     *
     * Reference: Shue et al. (1997) J. Geophys. Res. 102:9497
     *            Roelof & Sibeck (1993) J. Geophys. Res. 98:21421 */
    if (!x || !r) return;
    if (R_mp <= 0.0) { *x = 0.0; *r = 0.0; return; }

    if (theta <= M_PI/2.0) {
        /* Dayside: r = R_mp * (2/(1+cos(theta)))^0.5 */
        double shape = pow(2.0/(1.0+cos(theta)), 0.5);
        *r = R_mp * shape * sin(theta);
        *x = R_mp * shape * cos(theta);
    } else {
        /* Nightside: cylindrical tail */
        *x = -R_tail * cos(theta - M_PI/2.0);
        *r = R_tail;
    }
}

void tsyganenko_simple(double x, double z, double Dst, double B_out[3]) {
    /* Simplified Tsyganenko-like magnetospheric magnetic field model (2D).
     * Includes dipole + tail current sheet + ring current depression.
     *
     * B_x = B_dipole_x + B_tail_x
     * B_z = B_dipole_z + B_tail_z + B_ring_z
     *
     * Tail field: B_tail_x = B_t0 * tanh(z/H) * exp(x/L_tail)  for x < 0
     * Ring current: Dst depression in B_z (uniform dipole-like)
     *
     * Reference: Tsyganenko (1987) Planet. Space Sci. 35:1347
     *            Tsyganenko & Sitnov (2005) J. Geophys. Res. 110:A03208 */
    if (!B_out) return;

    double xm = x * SP_RE;  /* convert to meters */
    double zm = z * SP_RE;
    double r2 = xm*xm + zm*zm;
    double r = sqrt(r2);
    if (r < SP_RE) r = SP_RE;

    /* Dipole field at (x,z) in 2D (magnetospheric coordinates)
     * Using dipole with moment along -z */
    double r5 = r2 * r2 * r;
    double M = SP_MU_EARTH;
    double coeff_dip = SP_MU0 * M / (4.0 * M_PI * r5);

    B_out[0] = coeff_dip * 3.0 * xm * (-zm);  /* Bx from B = (mu0/(4pi r^5))[3(m.r)r - m r^2] */
    B_out[1] = 0.0;  /* By = 0 in 2D noon-midnight meridian */
    B_out[2] = coeff_dip * (3.0 * (-zm) * zm - M*r2/M); /* Bz simplified */

    /* Tail current sheet (for x < 0, nightside): B_x enhancement */
    if (x < 0.0) {
        double H_tail = 5.0 * SP_RE;   /* current sheet half-thickness */
        double L_tail = 30.0 * SP_RE;  /* tail length scale */
        double B_t0 = 20.0e-9;         /* tail lobe field at x~-20 R_E */
        double tail_factor = B_t0 * tanh(zm / H_tail) * exp(xm / L_tail);
        B_out[0] += tail_factor;
    }

    /* Ring current depression (uniform Dst in inner magnetosphere) */
    if (r < 10.0 * SP_RE) {
        double Dst_T = Dst * 1e-9;  /* convert nT to T */
        double rc_factor = exp(-r2/(2.0*4.0*SP_RE*4.0*SP_RE));
        B_out[2] += Dst_T * rc_factor;  /* Dst depresses Bz at Earth */
    }
}

/*========================================================================
 * L2: Ring Current
 *========================================================================*/

double dessler_parker_sckopke(double total_energy) {
    /* Dessler-Parker-Sckopke relation:
     * Delta_B [nT] = -3.98e-30 * U_R [J]
     * where U_R is total ring current particle energy.
     *
     * For a major storm with Dst = -200 nT:
     * U_R ~ 5e22 J (equivalent to ~12 megatons TNT)
     *
     * This is the fundamental energy-magnetic field relation for
     * magnetically trapped particle populations.
     *
     * Reference: Dessler & Parker (1959) J. Geophys. Res. 64:2239
     *            Sckopke (1966) J. Geophys. Res. 71:3125 */
    return -3.98e-30 * total_energy;
}

double ring_current_energy_density(double n, double T_perp, double T_par) {
    /* Ring current energy density for bi-Maxwellian plasma.
     *
     * U_R = (3/2)*n*kB*T_perp + (1/2)*n*kB*T_par
     *     = n*kB*T_perp*(3/2 + T_par/(2*T_perp))
     *
     * Including pressure anisotropy: A = T_perp/T_par - 1
     * U_R = (3/2)*n*kB*T*(1 + (2/3)*A) where T = (2*T_perp+T_par)/3
     *
     * Anisotropic distributions (A>0, mirror unstable) have enhanced energy.
     *
     * Reference: Daglis et al. (1999) Rev. Geophys. 37:407
     *            Williams "Ring Current and Radiation Belts" (1987) */
    if (n <= 0.0) return 0.0;
    return n * SP_KB * (1.5*T_perp + 0.5*T_par);
}

/*========================================================================
 * L2: Plasmapause & Plasmasphere
 *========================================================================*/

double plasmapause_L(double Kp_max) {
    /* Plasmapause location from Carpenter & Anderson (1992):
     * L_pp = 5.6 - 0.46 * Kp_max
     *
     * Kp = 0-1 (quiet):  L_pp ~ 5.1-5.6 R_E
     * Kp = 4-5 (active): L_pp ~ 3.3-3.8 R_E
     * Kp = 8-9 (storm):  L_pp ~ 1.5-1.9 R_E
     *
     * The plasmapause forms where co-rotation and convection
     * electric fields produce a stagnation point.
     *
     * Reference: Carpenter & Anderson (1992) J. Geophys. Res. 97:1097 */
    if (Kp_max < 0.0) Kp_max = 0.0;
    if (Kp_max > 9.0) Kp_max = 9.0;
    double Lpp = 5.6 - 0.46 * Kp_max;
    return (Lpp > 1.0) ? Lpp : 1.0;
}

double plasmasphere_density(double L, double L_pp, double n0,
                            double n_trough, double alpha) {
    /* Plasmasphere density profile.
     *
     * Inside plasmapause (L < L_pp): n_e(L) = n0 * (L_pp/L)^alpha
     * Outside plasmapause (L > L_pp): n_e(L) = n_trough
     *
     * Typical: n0=1000 cm^-3, alpha=3-4, n_trough=1 cm^-3.
     *
     * This sharp density gradient at the plasmapause supports
     * various wave-particle interactions (whistlers, EMIC waves).
     *
     * Reference: Carpenter & Anderson (1992)
     *            Gallagher et al. (2000) Space Sci. Rev. 91:235 */
    if (L <= 0.0 || L_pp <= 0.0) return 0.0;
    if (L < L_pp) {
        return n0 * pow(L_pp / L, alpha);
    } else {
        return n_trough;
    }
}

/*========================================================================
 * L6: Bow Shock
 *========================================================================*/

double bow_shock_standoff(double R_mp, double M_ms, double gamma) {
    /* Bow shock standoff distance using gasdynamic theory.
     *
     * Farris & Russell (1994):
     *   Delta/R_mp = 0.8 * (rho_up/rho_down - 1)^(-1)
     *
     * For strong shocks (M >> 1): rho_down/rho_up -> (gamma+1)/(gamma-1)
     * For gamma=5/3: rho_down/rho_up ~ 4, Delta/R_mp ~ 0.27
     * So R_bs ~ 1.27 * R_mp
     *
     * More generally: R_bs/R_mp = 1 + 1.1*(gamma-1)/((gamma+1)(M_ms^2-1))
     *
     * Reference: Farris & Russell (1994) J. Geophys. Res. 99:17681
     *            Spreiter et al. (1966) Planet. Space Sci. 14:223 */
    if (R_mp <= 0.0) return 0.0;
    if (M_ms <= 1.0) return R_mp * 1.5;  /* no shock for sub-fast flow */
    double M2 = M_ms * M_ms;
    double rho_ratio = (gamma + 1.0) * M2 / ((gamma - 1.0) * M2 + 2.0);
    double delta_ratio = 1.1 / (rho_ratio - 1.0);
    return R_mp * (1.0 + delta_ratio);
}

void mhd_shock_jump(double M1, double gamma, double *rho_ratio,
                    double *p_ratio, double *B_ratio) {
    /* Rankine-Hugoniot jump conditions for perpendicular MHD shock.
     *
     * For shock normal perpendicular to B:
     *   rho_2/rho_1 = v_1/v_2 = (gamma+1)M1^2 / ((gamma-1)M1^2 + 2)
     *   p_2/p_1     = (2*gamma*M1^2 - (gamma-1)) / (gamma+1)
     *   B_2/B_1     = rho_2/rho_1  (flux freezing, perpendicular shock)
     *
     * For M1 -> infinity (strong shock):
     *   rho_ratio -> (gamma+1)/(gamma-1) = 4 for gamma=5/3
     *   p_ratio -> infinity (quadratic in M1)
     *
     * Reference: Tidman & Krall "Shock Waves in Collisionless Plasmas" (1971)
     *            Kivelson & Russell S7.3 */
    if (!rho_ratio || !p_ratio || !B_ratio) return;
    if (M1 <= 1.0) {
        *rho_ratio = 1.0;
        *p_ratio   = 1.0;
        *B_ratio   = 1.0;
        return;
    }

    double M12 = M1 * M1;
    double gp1 = gamma + 1.0;
    double gm1 = gamma - 1.0;

    *rho_ratio = gp1 * M12 / (gm1 * M12 + 2.0);
    *p_ratio   = (2.0 * gamma * M12 - gm1) / gp1;
    *B_ratio   = *rho_ratio;  /* perpendicular shock: B proportional to rho */
}

void magnetosheath_flow(double x, double y, double R_mp,
                        double *vx_out, double *vy_out) {
    /* Magnetosheath flow deflection around magnetopause obstacle.
     *
     * Using gasdynamic potential flow approximation:
     * Flow is deflected around the magnetopause with velocity:
     *   v_x ~ 1 - (R_mp/r)^2 * cos(2*phi) (to first order)
     * where r = sqrt(x^2+y^2), phi = atan2(y,x).
     *
     * At the nose (x=R_mp, y=0): v_x = 0 (stagnation point)
     * At the flank (x=0, y=R_mp): v_y increases
     *
     * Reference: Spreiter et al. (1966) Planet. Space Sci. 14:223
     *            Kivelson & Russell S7.4 */
    if (!vx_out || !vy_out) return;

    double r2 = x*x + y*y;
    double r = sqrt(r2);
    if (r < 1e-10) { *vx_out = 0.0; *vy_out = 0.0; return; }

    double cos_phi = x / r;
    double sin_phi = y / r;
    double R_mp_r = (r > 0.0) ? R_mp / r : 0.0;
    double R_mp_r3 = R_mp_r * R_mp_r * R_mp_r;

    /* Potential flow around sphere: v_r = (1-R^3/r^3)*cos(phi) */
    double v_r = (1.0 - R_mp_r3) * cos_phi;

    /* Approximate 2D deflection (cylindrical obstacle) */
    double v_theta = -(1.0 + 0.5*R_mp_r3) * sin_phi;

    *vx_out = v_r * cos_phi - v_theta * sin_phi;
    *vy_out = v_r * sin_phi + v_theta * cos_phi;
}

/*========================================================================
 * L6: Particle Drifts in Magnetosphere
 *========================================================================*/

void gradB_drift(double v_perp, const double B[3], const double gradB[3],
                 double q, double m, double v_drift[3]) {
    /* Gradient-B drift velocity:
     * v_gradB = (m*v_perp^2 / (2*q*B^3)) * B x grad(B)
     *
     * Direction: ions and electrons drift in opposite directions.
     * In Earth's dipole: ions drift westward, electrons eastward.
     * Drift speed increases with L-shell (~L^2).
     *
     * Reference: Northrop "The Adiabatic Motion of Charged Particles" (1963)
     *            Goldston & Rutherford S2.5 */
    if (!B || !gradB || !v_drift) return;
    if (q == 0.0 || m <= 0.0) {
        v_drift[0] = v_drift[1] = v_drift[2] = 0.0;
        return;
    }

    double Bmag = SP_MAG3(B);
    if (Bmag < 1e-30) {
        v_drift[0] = v_drift[1] = v_drift[2] = 0.0;
        return;
    }
    double B3 = Bmag * Bmag * Bmag;

    double coeff = m * v_perp * v_perp / (2.0 * q * B3);

    /* B x gradB */
    SP_CROSS3(B, gradB, v_drift);
    v_drift[0] *= coeff;
    v_drift[1] *= coeff;
    v_drift[2] *= coeff;
}

void curvature_drift(double v_par, const double B[3], const double curv[3],
                     double q, double m, double v_drift[3]) {
    /* Curvature drift velocity:
     * v_curv = (m*v_par^2 / (q*B^2)) * B x (b.grad)b
     * where (b.grad)b = curvature vector.
     *
     * Combined with grad-B drift, the total guiding center drift for
     * isotropic pressure in a dipole is westward for ions, eastward for electrons.
     *
     * Reference: Northrop (1963), Roederer "Dynamics of Geomagnetically
     *            Trapped Radiation" (1970) */
    if (!B || !curv || !v_drift) return;
    if (q == 0.0 || m <= 0.0) {
        v_drift[0] = v_drift[1] = v_drift[2] = 0.0;
        return;
    }

    double Bmag = SP_MAG3(B);
    if (Bmag < 1e-30) {
        v_drift[0] = v_drift[1] = v_drift[2] = 0.0;
        return;
    }
    double B2 = Bmag * Bmag;

    double coeff = m * v_par * v_par / (q * B2);

    /* B x curv */
    SP_CROSS3(B, curv, v_drift);
    v_drift[0] *= coeff;
    v_drift[1] *= coeff;
    v_drift[2] *= coeff;
}

void dipole_drift_speed(double E_kin, double L, double q, double v_drift[3]) {
    /* Total guiding center drift speed in dipole field.
     *
     * For equatorially mirroring particles in a dipole:
     * v_d ~ (3*E_kin*L^2) / (q*B_E*R_E)
     *
     * Drift period: tau_d ~ 2*pi*L*R_E / v_d
     * For 1 MeV electrons at L=4: tau_d ~ hours
     * For 100 keV protons at L=6: tau_d ~ hours
     *
     * Reference: Hamlin et al. (1961) J. Geophys. Res. 66:1
     *            Schulz & Lanzerotti (1974) */
    if (!v_drift) return;
    if (L <= 0.0 || q == 0.0) {
        v_drift[0] = v_drift[1] = v_drift[2] = 0.0;
        return;
    }

    /* Azimuthal drift speed (phi direction) */
    double B_E = SP_B0_EARTH;
    double v_d = 3.0 * E_kin * L * L / (fabs(q) * B_E * SP_RE);

    /* Direction: eastward for electrons (q<0), westward for ions (q>0) */
    double sign = (q < 0.0) ? 1.0 : -1.0;

    v_drift[0] = 0.0;
    v_drift[1] = sign * v_d;  /* azimuthal (phi) direction */
    v_drift[2] = 0.0;
}

double dipole_bounce_period(double E_kin, double m, double L, double alpha_eq) {
    /* Bounce period for particles in dipole magnetic field.
     *
     * tau_b ~ (4*L*R_E/v) * (1.30 - 0.56*sin(alpha_eq))
     *
     * where v = sqrt(2*E_kin/m) is particle speed,
     * alpha_eq is equatorial pitch angle.
     *
     * Bounce period is independent of particle charge (only E, m).
     * For 1 MeV electrons at L=4: tau_b ~ 0.5 s
     * For 100 keV protons at L=4: tau_b ~ 10 s
     *
     * Reference: Hamlin et al. (1961)
     *            Schulz & Lanzerotti (1974)
     *            Roederer (1970) */
    if (E_kin <= 0.0 || m <= 0.0 || L <= 0.0) return 0.0;

    double v = sqrt(2.0 * E_kin / m);
    if (v <= 0.0) return INFINITY;

    double sin_alpha = sin(alpha_eq);
    double geom_factor = 1.30 - 0.56 * sin_alpha;

    return (4.0 * L * SP_RE / v) * geom_factor;
}
