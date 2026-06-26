#ifndef INSTABILITIES_H
#define INSTABILITIES_H

/*
 * instabilities.h -- Parametric instabilities in laser-plasma interaction
 *
 * Models for the principal laser-driven instabilities that occur in
 * underdense plasma: Stimulated Raman Scattering (SRS), Stimulated
 * Brillouin Scattering (SBS), Two-Plasmon Decay (TPD), and
 * filamentation.
 *
 * References:
 *   - Kruer (1988), Ch. 7-8
 *   - Liu, Rosenbluth & White (1974), Phys. Fluids 17, 1211
 *   - Drake et al. (1974), Phys. Fluids 17, 778
 *   - Afeyan & Williams (1985), Phys. Fluids 28, 3397
 *
 * Knowledge Layers:
 *   L1: growth rates, thresholds, matching conditions
 *   L2: SRS, SBS, TPD instability mechanisms
 *   L4: Three-wave coupling equations, Manley-Rowe relations
 *   L5: Convective/absolute instability analysis
 *
 * Courses: MIT 22.611, Princeton PHY 535, Caltech Ph 230
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L2: Instability Type Definitions
 * ============================================================ */

/** InstabilityType enumerates the parametric processes */
typedef enum {
    INSTABILITY_SRS_BACKSCATTER = 0,  /* Stimulated Raman backscatter  */
    INSTABILITY_SRS_FORWARD     = 1,  /* Stimulated Raman forward     */
    INSTABILITY_SRS_SIDESCATTER = 2,  /* Stimulated Raman sidescatter */
    INSTABILITY_SBS_BACKSCATTER = 3,  /* Stimulated Brillouin back    */
    INSTABILITY_SBS_FORWARD     = 4,  /* Stimulated Brillouin forward */
    INSTABILITY_TPD             = 5,  /* Two-plasmon decay             */
    INSTABILITY_FILAMENTATION   = 6,  /* Relativistic filamentation    */
    INSTABILITY_MODULATIONAL    = 7,  /* Modulational instability      */
    INSTABILITY_COUNT           = 8
} InstabilityType;

/**
 * InstabilityResult -- output of instability analysis
 *
 * gamma0     : homogeneous growth rate [s^-1]
 * gamma_max  : maximum growth rate (with gradient) [s^-1]
 * threshold  : threshold intensity [W/m^2] or amplitude
 * convective : 1 = convective, 0 = absolute
 * k_epw      : wavenumber of the excited EPW [rad/m]
 * omega_epw  : frequency of excited EPW [rad/s]
 * k_iaw      : wavenumber of the IAW (SBS only) [rad/m]
 * omega_iaw  : frequency of IAW (SBS only) [rad/s]
 * rosenbluth_gain : convective gain exponent
 */
typedef struct {
    double gamma0;
    double gamma_max;
    double threshold;
    int    convective;
    double k_epw;
    double omega_epw;
    double k_iaw;
    double omega_iaw;
    double rosenbluth_gain;
} InstabilityResult;

/* ============================================================
 *  L1/L2: Phase-Matching Conditions
 * ============================================================ */

/**
 * srs_matching -- SRS frequency and wavenumber matching
 *
 *   omega0 = omega_s + omega_epw  (energy conservation)
 *   k0 = k_s + k_epw             (momentum conservation)
 *
 * Solves the matching conditions for SRS backscatter. Given incident
 * laser frequency omega0 and wavenumber k0, returns the scattered
 * light frequency omega_s and EPW wavenumber k_epw.
 *
 * Complexity: O(1)
 * Theorem: Three-wave parametric coupling from Maxwell-fluid equations
 */
void srs_matching(double omega0, double k0, double wp,
                  double *omega_s, double *k_epw);

/**
 * sbs_matching -- SBS frequency and wavenumber matching
 *
 *   omega0 = omega_s + omega_iaw
 *   k0 = k_s + k_iaw
 *
 * Frequency matching for stimulated Brillouin scattering.
 * Returns scattered light frequency omega_s and IAW frequency
 * omega_iaw, wavenumber k_iaw.
 */
void sbs_matching(double omega0, double k0, double cs,
                  double *omega_s, double *omega_iaw, double *k_iaw);

/**
 * tpd_matching -- Two-plasmon decay matching
 *
 *   omega0 = omega_epw1 + omega_epw2
 *   k0 = k_epw1 + k_epw2
 *
 * Laser photon decays into two electron plasma waves at the
 * quarter-critical surface. Returns EPW frequencies and wavenumbers.
 */
void tpd_matching(double omega0, double k0, double wp,
                  double *omega1, double *omega2,
                  double *k1, double *k2);

/* ============================================================
 *  L1/L4: Growth Rates
 * ============================================================ */

/**
 * srs_growth_rate_backscatter -- SRS backscatter growth rate
 *
 *   gamma0 = (k_epw v_osc / 4) sqrt(wp / (omega_s - wp^2/omega_s))
 *
 * Homogeneous temporal growth rate for SRS backscatter.
 * [Kruer Sec 7.2; Liu et al. 1974]
 */
double srs_growth_rate_backscatter(double k_epw, double v_osc,
                                    double wp, double omega0);

/**
 * srs_growth_rate_forward -- SRS forward scatter growth rate
 *
 *   gamma0 = (wp^2 / (4 omega0)) (v_osc / v_the)
 *
 * Forward SRS has typically lower growth rate than backscatter
 * but can be dangerous in large, underdense plasmas.
 */
double srs_growth_rate_forward(double wp, double omega0,
                                double v_osc, double v_the);

/**
 * sbs_growth_rate -- SBS growth rate
 *
 *   gamma0 = (k_iaw v_osc / 4) sqrt(omega_pi^2 / (omega_s omega_iaw))
 *
 * SBS is an ion-timescale instability. Usually lower growth rate
 * than SRS but lower threshold since it involves low-frequency IAW.
 * [Kruer Sec 7.3]
 */
double sbs_growth_rate(double k_iaw, double v_osc, double wpi,
                        double omega_s, double omega_iaw);

/**
 * tpd_growth_rate -- Two-plasmon decay growth rate
 *
 *   gamma0 = (k_epw v_osc / 4)
 *
 * Maximum at the quarter-critical surface (ne = nc/4) where
 * omega0 = 2 wp.
 * [Kruer Sec 7.4; Drake et al. 1974]
 */
double tpd_growth_rate(double k_epw, double v_osc);

/**
 * srs_convective_gain -- Rosenbluth gain for convective SRS
 *
 *   G = 2 pi gamma0^2 / (|v_g1 v_g2| kappa')
 *
 * where kappa' is the wavenumber mismatch gradient. G > 2pi
 * indicates significant convective amplification.
 * [Rosenbluth 1972]
 */
double srs_convective_gain(double gamma0, double vg_s, double vg_epw,
                            double kappa_prime);

/**
 * sbs_convective_gain -- Rosenbluth gain for convective SBS
 */
double sbs_convective_gain(double gamma0, double vg_s, double vg_iaw,
                            double kappa_prime);

/* ============================================================
 *  L2/L4: Threshold Conditions
 * ============================================================ */

/**
 * srs_absolute_threshold -- Threshold for absolute SRS
 *
 *   (v_osc^2 / c^2) > (8 / (k0 L_n)) * (nu_epw / omega0)
 *
 * Threshold for transition from convective to absolute SRS.
 * [Liu et al. 1974; Kruer Sec 7.2]
 */
double srs_absolute_threshold_a0(double L_n, double lambda0,
                                  double nu_epw, double omega0);

/**
 * sbs_threshold_intensity -- SBS threshold intensity
 *
 *   I_th_SBS ~ 4e15 * (Te_keV) / (L_n_mum * lambda_mum)  [W/cm^2]
 *
 * Approximate threshold for SBS in ICF-relevant plasmas.
 */
double sbs_threshold_intensity(double Te_keV, double L_n_mum,
                                double lambda_mum);

/**
 * tpd_threshold -- Two-plasmon decay threshold
 *
 *   (v_osc / v_the)^2 > 12 / (k0 L_n)
 *
 * Threshold at quarter-critical surface.
 * [Drake et al. 1974]
 */
double tpd_threshold(double v_osc, double v_the, double k0, double L_n);

/* ============================================================
 *  L5: Full Instability Analysis
 * ============================================================ */

/**
 * analyze_instability -- Compute instability parameters
 *
 * Given plasma state and laser parameters, compute the growth
 * rate, threshold, and matching conditions for a specified
 * instability type.
 *
 * Returns 0 on success, -1 if the instability is physically
 * impossible in the given conditions.
 *
 * Complexity: O(1)
 */
int analyze_instability(InstabilityType type,
                         const PlasmaState *ps,
                         const PlasmaDerived *pd,
                         double v_osc,
                         double density_scale_length,
                         InstabilityResult *result);

/**
 * instability_name -- Human-readable name for InstabilityType
 */
const char *instability_name(InstabilityType type);

#endif /* INSTABILITIES_H */
