/*
 * instabilities.c -- Parametric instabilities in laser-plasma interaction
 *
 * Implements growth rates, thresholds, and matching conditions
 * for SRS, SBS, TPD, and filamentation instabilities.
 *
 * References:
 *   - Kruer (1988), Ch. 7-8
 *   - Liu, Rosenbluth & White (1974), Phys. Fluids 17, 1211
 *   - Drake et al. (1974), Phys. Fluids 17, 778
 *   - Forslund, Kindel & Lindman (1975), Phys. Fluids 18, 1002
 *
 * Knowledge Layers: L1, L2, L4, L5
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "laser_plasma.h"
#include "laser_propagation.h"
#include "instabilities.h"

static const double c   = PLASMA_C;

/* ============================================================
 *  L2: srs_matching -- Stimulated Raman Scattering matching
 *
 * Three-wave parametric coupling:
 *   omega0 = omega_s + omega_epw   (energy)
 *   k0     = k_s     + k_epw       (momentum)
 *
 * For backscatter SRS:
 *   omega_s ~ omega0 - omega_p (Stokes shift)
 *   k_epw ~ 2 k0              (short-wavelength EPW)
 *
 * The EPW carries the momentum difference between incident
 * and scattered photons.
 *
 * Source: Kruer Sec 7.2; Liu et al. 1974
 *
 * Complexity: O(1)
 * Theorem: Manley-Rowe relations for parametric three-wave
 *          coupling: action conservation N0 + N_s = const
 * ============================================================ */
void srs_matching(double omega0, double k0, double wp,
                  double *omega_s, double *k_epw)
{
    if (!omega_s || !k_epw || omega0 <= 0.0 || wp >= omega0) {
        if (omega_s) *omega_s = 0.0;
        if (k_epw)  *k_epw = 0.0;
        return;
    }

    /* For backscatter: omega_s = omega0 - omega_p (approx) */
    double omega_s_val = omega0 - wp;

    /* And from EM dispersion: omega_s^2 = wp^2 + c^2 k_s^2 */
    double k_s = sqrt(omega_s_val * omega_s_val - wp * wp) / c;
    /* k_epw = k0 + k_s (backscatter: k_s negative, k_epw > k0) */
    double k_epw_val = k0 + k_s;

    *omega_s = omega_s_val;
    *k_epw  = k_epw_val;
}

/* ============================================================
 *  L2: sbs_matching -- Stimulated Brillouin Scattering matching
 *
 *   omega0 = omega_s + omega_iaw
 *   k0     = k_s     + k_iaw
 *
 * For backscatter SBS:
 *   omega_iaw ~ 2 omega0 cs/c (ion acoustic frequency)
 *   k_iaw     ~ 2 k0          (this is the sound wave)
 *
 * SBS involves a low-frequency ion acoustic wave (IAW) rather
 * than an electron plasma wave.  The frequency shift is very
 * small (~10^-3 omega0) but the growth can be significant.
 *
 * Source: Kruer Sec 7.3; Forslund et al. 1975
 *
 * Complexity: O(1)
 * ============================================================ */
void sbs_matching(double omega0, double k0, double cs,
                  double *omega_s, double *omega_iaw, double *k_iaw)
{
    if (!omega_s || !omega_iaw || !k_iaw || omega0 <= 0.0) {
        if (omega_s)  *omega_s = 0.0;
        if (omega_iaw) *omega_iaw = 0.0;
        if (k_iaw)    *k_iaw = 0.0;
        return;
    }

    /* IAW wavenumber ~ 2 k0 (backscatter) */
    double k_iaw_val = 2.0 * k0;
    /* IAW frequency: omega_iaw = k_iaw * cs */
    double omega_iaw_val = k_iaw_val * cs;

    /* Scattered light frequency */
    double omega_s_val = omega0 - omega_iaw_val;

    *omega_s   = omega_s_val;
    *omega_iaw = omega_iaw_val;
    *k_iaw     = k_iaw_val;
}

/* ============================================================
 *  L2: tpd_matching -- Two-Plasmon Decay matching
 *
 *   omega0 = omega_epw1 + omega_epw2
 *   k0     = k_epw1     + k_epw2
 *
 * At the quarter-critical surface (ne = nc/4 => omega_p = omega0/2),
 * the laser photon can decay into two electron plasma waves,
 * each at half the laser frequency.
 *
 * Source: Kruer Sec 7.4; Drake et al. 1974
 *
 * Complexity: O(1)
 * ============================================================ */
void tpd_matching(double omega0, double k0, double wp,
                  double *omega1, double *omega2,
                  double *k1, double *k2)
{
    if (!omega1 || !omega2 || !k1 || !k2) return;

    /* TPD requires omega0 = 2 omega_p => wp = omega0/2 */
    *omega1 = omega0 / 2.0;
    *omega2 = omega0 / 2.0;

    /* Symmetric configuration: k1 = k0/2 + k_perp, k2 = k0/2 - k_perp */
    *k1 = k0 / 2.0;
    *k2 = k0 / 2.0;
    (void)wp;   /* wp used only for TPD existence condition */
}

/* ============================================================
 *  L1: srs_growth_rate_backscatter
 *
 *  gamma0 = (k_epw v_osc / 4) * sqrt(wp / omega0)
 *
 * Homogeneous (infinite, uniform plasma) temporal growth rate
 * for SRS backscatter.
 *
 * For typical ICF parameters at nc/10:
 *   k_epw ~ 2 k0, v_osc ~ 0.1c (for I~10^16 W/cm^2)
 *   gamma0 ~ 10^13 s^-1  (growth time ~ 100 fs)
 *
 * Source: Kruer Sec 7.2; Liu et al. 1974, Eq. 20
 *
 * Complexity: O(1)
 * ============================================================ */
double srs_growth_rate_backscatter(double k_epw, double v_osc,
                                    double wp, double omega0)
{
    if (omega0 <= 0.0 || wp > omega0) return 0.0;
    return (k_epw * v_osc / 4.0) * sqrt(wp / omega0);
}

/* ============================================================
 *  L1: srs_growth_rate_forward
 *
 *  gamma0 = (wp^2 / (4 omega0)) * (v_osc / v_the)
 *
 * Forward SRS has lower growth rate than backscatter by
 * factor ~ (wp/omega0) * (c/v_the).  However, it can be
 * dangerous in large, hot plasmas where the convective
 * loss is small.
 *
 * Source: Kruer Sec 7.2
 *
 * Complexity: O(1)
 * ============================================================ */
double srs_growth_rate_forward(double wp, double omega0,
                                double v_osc, double v_the)
{
    if (omega0 <= 0.0 || v_the <= 0.0 || wp > omega0) return 0.0;
    return (wp * wp) / (4.0 * omega0) * (v_osc / v_the);
}

/* ============================================================
 *  L1: sbs_growth_rate -- SBS growth rate
 *
 *  gamma0 = k_iaw v_osc / (4 sqrt(2)) * sqrt(omega_pi / omega0)
 *
 * SBS growth rate is typically lower than SRS by factor
 * sqrt(me/mi) ~ 1/43, but SBS has lower threshold because
 * the IAW damping is much weaker than EPW damping.
 *
 * Source: Kruer Sec 7.3; Forslund et al. 1975
 *
 * Complexity: O(1)
 * ============================================================ */
double sbs_growth_rate(double k_iaw, double v_osc, double wpi,
                        double omega_s, double omega_iaw)
{
    if (omega_s <= 0.0 || omega_iaw <= 0.0) return 0.0;
    return (k_iaw * v_osc / 4.0)
           * sqrt(wpi * wpi / (omega_s * omega_iaw));
}

/* ============================================================
 *  L1: tpd_growth_rate -- Two-plasmon decay growth rate
 *
 *  gamma0 = k_epw v_osc / 4
 *
 * Maximum at quarter-critical surface (ne = nc/4) where
 * omega0 = 2 omega_p.
 *
 * TPD is a particular concern for direct-drive ICF because
 * it generates hot electrons that preheat the fuel.
 *
 * Source: Kruer Sec 7.4; Drake et al. 1974
 *
 * Complexity: O(1)
 * ============================================================ */
double tpd_growth_rate(double k_epw, double v_osc)
{
    return k_epw * v_osc / 4.0;
}

/* ============================================================
 *  L2: srs_convective_gain -- Rosenbluth gain exponent
 *
 *  G = 2 pi gamma0^2 / (|v_g1 v_g2| kappa')
 *
 * where kappa' = d/dx(k0 - k_s - k_epw) is the wavenumber
 * mismatch gradient due to plasma inhomogeneity.
 *
 * Convective SRS is amplified as the scattered light and EPW
 * propagate through the interaction region.  G >> 1 indicates
 * significant amplification.
 *
 * Source: Rosenbluth (1972), Phys. Rev. Lett. 29, 565
 *
 * Complexity: O(1)
 * ============================================================ */
double srs_convective_gain(double gamma0, double vg_s, double vg_epw,
                            double kappa_prime)
{
    if (vg_s <= 0.0 || vg_epw <= 0.0 || kappa_prime <= 0.0) return 0.0;
    return 2.0 * M_PI * gamma0 * gamma0
           / (vg_s * vg_epw * kappa_prime);
}

/* ============================================================
 *  L2: sbs_convective_gain -- Rosenbluth gain for SBS
 *
 * Same formula as SRS but with the scattered light and IAW
 * group velocities.
 *
 * Complexity: O(1)
 * ============================================================ */
double sbs_convective_gain(double gamma0, double vg_s, double vg_iaw,
                            double kappa_prime)
{
    if (vg_s <= 0.0 || vg_iaw <= 0.0 || kappa_prime <= 0.0) return 0.0;
    return 2.0 * M_PI * gamma0 * gamma0
           / (vg_s * vg_iaw * kappa_prime);
}

/* ============================================================
 *  L2: srs_absolute_threshold_a0
 *
 *  (v_osc/c)^2 > (8 / (k0 L_n)) * (nu_epw / omega0)
 *
 * The transition from convective to absolute SRS occurs when
 * the amplification in the inhomogeneous plasma exceeds the
 * damping.  Above this threshold, SRS grows exponentially
 * in time at a fixed position.
 *
 * Absolute SRS is particularly dangerous for ICF because it
 * can reach high levels and produce energetic electrons.
 *
 * Source: Liu et al. 1974; Kruer Sec 7.2
 *
 * Complexity: O(1)
 * ============================================================ */
double srs_absolute_threshold_a0(double L_n, double lambda0,
                                  double nu_epw, double omega0)
{
    if (L_n <= 0.0 || lambda0 <= 0.0 || omega0 <= 0.0) return DBL_MAX;

    double k0 = 2.0 * M_PI / lambda0;
    return sqrt((8.0 / (k0 * L_n)) * (nu_epw / omega0));
}

/* ============================================================
 *  L2: sbs_threshold_intensity
 *
 *  I_th_SBS ~ 4e15 * (Te_keV) / (L_n_um * lambda_um)  [W/cm^2]
 *
 * Approximate threshold for observable SBS in ICF plasmas.
 * SBS becomes important when the intensity exceeds this
 * value.
 *
 * Source: empirical scaling from ICF experiments
 *
 * Complexity: O(1)
 * ============================================================ */
double sbs_threshold_intensity(double Te_keV, double L_n_um,
                                double lambda_um)
{
    if (Te_keV <= 0.0 || L_n_um <= 0.0 || lambda_um <= 0.0)
        return DBL_MAX;
    return 4.0e15 * Te_keV / (L_n_um * lambda_um);
}

/* ============================================================
 *  L2: tpd_threshold -- Two-plasmon decay threshold
 *
 *  (v_osc / v_the)^2 > 12 / (k0 L_n)
 *
 * At the quarter-critical surface, TPD becomes absolutely
 * unstable above this threshold.
 *
 * Source: Drake et al. 1974
 *
 * Complexity: O(1)
 * ============================================================ */
double tpd_threshold(double v_osc, double v_the, double k0, double L_n)
{
    (void)v_osc;  /* threshold expressed in terms of v_the */
    if (k0 <= 0.0 || L_n <= 0.0) return DBL_MAX;
    double rhs = 12.0 / (k0 * L_n);
    return sqrt(rhs) * v_the;  /* returns threshold v_osc */
}

/* ============================================================
 *  L5: analyze_instability -- Complete instability characterization
 *
 * Given plasma state and laser parameters, compute growth
 * rate, threshold, and matching conditions for a specified
 * instability type.
 *
 * Returns 0 on success, -1 if instability is not possible.
 *
 * Complexity: O(1)
 * ============================================================ */
int analyze_instability(InstabilityType type,
                         const PlasmaState *ps,
                         const PlasmaDerived *pd,
                         double v_osc,
                         double density_scale_length,
                         InstabilityResult *result)
{
    if (!ps || !pd || !result) return -1;

    memset(result, 0, sizeof(InstabilityResult));
    result->convective = 1;  /* most are convective by default */

    double omega0 = pd->omega_laser;
    double k0 = 2.0 * M_PI / ps->lam_laser;
    double wp  = pd->wp;
    double omega_s, k_epw, k_iaw, omega_iaw;

    switch (type) {
    case INSTABILITY_SRS_BACKSCATTER:
        if (wp >= omega0) return -1;
        srs_matching(omega0, k0, wp, &omega_s, &k_epw);
        result->k_epw = k_epw;
        result->omega_epw = omega0 - omega_s;
        result->gamma0 = srs_growth_rate_backscatter(k_epw, v_osc,
                                                      wp, omega0);
        if (density_scale_length > 0.0) {
            double kappa_p = 1.0 / density_scale_length;
            double vg_s = c * sqrt(1.0 - wp*wp/(omega_s*omega_s));
            double vg_epw = 3.0 * pd->v_the * pd->v_the * k_epw
                            / result->omega_epw;
            if (vg_s > 0.0 && vg_epw > 0.0 && kappa_p > 0.0) {
                result->rosenbluth_gain = srs_convective_gain(
                    result->gamma0, vg_s, vg_epw, kappa_p);
            }
        }
        result->threshold = srs_absolute_threshold_a0(
            density_scale_length, ps->lam_laser,
            pd->nu_ei, omega0);
        if (result->threshold < 1e-6) result->convective = 0;
        break;

    case INSTABILITY_SRS_FORWARD:
        if (wp >= omega0 || pd->v_the <= 0.0) return -1;
        result->gamma0 = srs_growth_rate_forward(wp, omega0,
                                                  v_osc, pd->v_the);
        result->k_epw = wp / c;
        result->omega_epw = wp;
        break;

    case INSTABILITY_SBS_BACKSCATTER:
        sbs_matching(omega0, k0, pd->cs, &omega_s, &omega_iaw, &k_iaw);
        result->k_iaw = k_iaw;
        result->omega_iaw = omega_iaw;
        result->gamma0 = sbs_growth_rate(k_iaw, v_osc,
                                          pd->wpi, omega_s, omega_iaw);
        result->threshold = sbs_threshold_intensity(
            ps->Te * 1e-3, density_scale_length * 1e6,
            ps->lam_laser * 1e6);
        break;

    case INSTABILITY_TPD:
        if (fabs(ps->ne - pd->nc/4.0) > 0.1 * pd->nc) return -1;
        tpd_matching(omega0, k0, wp, &result->omega_epw,
                     &omega_iaw, &result->k_epw, &k_iaw);
        result->gamma0 = tpd_growth_rate(result->k_epw, v_osc);
        result->threshold = tpd_threshold(v_osc, pd->v_the, k0,
                                          density_scale_length);
        break;

    case INSTABILITY_FILAMENTATION:
        result->gamma0 = filamentation_growth_rate(wp, omega0,
                                                    pd->a0);
        result->threshold = relativistic_self_focusing_threshold(
            ps->ne, pd->nc);
        break;

    default:
        return -1;
    }

    return 0;
}

/* ============================================================
 *  L1: instability_name -- Human-readable names
 * ============================================================ */
const char *instability_name(InstabilityType type)
{
    static const char *names[] = {
        "SRS Backscatter",
        "SRS Forward Scatter",
        "SRS Sidescatter",
        "SBS Backscatter",
        "SBS Forward Scatter",
        "Two-Plasmon Decay",
        "Relativistic Filamentation",
        "Modulational Instability"
    };
    if (type < 0 || type >= INSTABILITY_COUNT) return "Unknown";
    return names[type];
}
