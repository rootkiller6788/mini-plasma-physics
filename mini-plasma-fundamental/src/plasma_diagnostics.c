/**
 * plasma_diagnostics.c -- Plasma Diagnostics Methods
 *
 * Langmuir probe analysis, interferometry, Thomson scattering,
 * spectroscopy, and fusion/space plasma diagnostics.
 *
 * References:
 *   - Hutchinson, "Principles of Plasma Diagnostics" (2002)
 *   - Bosch & Hale (1992), Nucl. Fusion 32, 611
 *   - MIT 22.611 / Princeton PHY 521
 *
 * Knowledge Coverage:
 *   L5: Langmuir probe, interferometry, Thomson scattering
 *   L7: Space plasma (magnetopause), fusion diagnostics (Lawson)
 */
#include "plasma_constants.h"
#include "plasma_params.h"
#include <math.h>

typedef struct {
    double V_float;
    double V_plasma;
    double I_sat_ion;
    double I_sat_e;
    double Te_eV;
    double ne;
} LangmuirProbeResult;

double density_from_ion_saturation(double I_sat, double A_probe,
                                   double Te, double mi) {
    if (A_probe <= 0.0 || mi <= 0.0) return 0.0;
    double cs = sqrt(K_B * Te / mi);
    return fabs(I_sat) / (0.61 * E_CHARGE * A_probe * cs);
}

/**
 * Langmuir probe I-V analysis using exponential fit.
 *
 * Extracts Te, ne, V_float, V_plasma from probe I-V data.
 *
 * Reference: Mott-Smith & Langmuir (1926), Phys. Rev. 28, 727.
 */
int langmuir_probe_analyze(const double *V, const double *I, int n_points,
                           double A_probe, double mi,
                           LangmuirProbeResult *result) {
    if (!V || !I || !result || n_points < 10) return -1;

    double I_sat_ion = 0.0;
    int n_ion = 0;
    for (int i = 0; i < n_points; i++) {
        if (V[i] < V[n_points/2] - 10.0) {
            I_sat_ion += I[i]; n_ion++;
        }
    }
    if (n_ion > 0) I_sat_ion /= n_ion;

    double V_float = 0.0;
    for (int i = 1; i < n_points; i++) {
        if (I[i-1] * I[i] <= 0.0) {
            V_float = 0.5 * (V[i-1] + V[i]); break;
        }
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    int n_fit = 0;
    for (int i = 0; i < n_points && n_fit < 50; i++) {
        if (V[i] > V_float && I[i] < 0.5 * fabs(I_sat_ion)) {
            double I_corr = I[i] - I_sat_ion;
            if (I_corr > 1e-15) {
                sum_x += V[i]; sum_y += log(I_corr);
                sum_xy += V[i] * log(I_corr); sum_x2 += V[i] * V[i];
                n_fit++;
            }
        }
    }
    if (n_fit < 3) return -1;

    double slope = (n_fit * sum_xy - sum_x * sum_y)
                   / (n_fit * sum_x2 - sum_x * sum_x);
    double Te_eV = 1.0 / slope;
    double Te_K = Te_eV * E_CHARGE / K_B;

    double cs = sqrt(K_B * Te_K / mi);
    double ne = fabs(I_sat_ion) / (0.61 * E_CHARGE * A_probe * cs);

    result->V_float = V_float;
    result->V_plasma = V_float + Te_eV * log(1.0 / 0.61);
    result->I_sat_ion = I_sat_ion;
    result->I_sat_e = 0.25 * ne * E_CHARGE * A_probe
                      * sqrt(8.0 * K_B * Te_K / (M_PI * M_ELECTRON));
    result->Te_eV = Te_eV;
    result->ne = ne;
    return 0;
}

/* ================================================================
 * L5: Microwave Interferometry
 * ================================================================ */

double plasma_refractive_index(double omega, double omega_pe) {
    if (omega <= omega_pe) return 0.0;
    return sqrt(1.0 - (omega_pe * omega_pe) / (omega * omega));
}

double interferometer_phase_shift_per_m(double omega, double ne) {
    double omega_pe = electron_plasma_frequency(ne);
    if (omega <= omega_pe) return INFINITY;
    return omega * omega_pe * omega_pe / (2.0 * omega * omega * C_LIGHT);
}

double line_integrated_density(double frequency, double delta_phi) {
    double omega = 2.0 * M_PI * frequency;
    return (2.0 * omega * C_LIGHT * EPSILON_0 * M_ELECTRON
            / (E_CHARGE * E_CHARGE)) * delta_phi;
}

double cutoff_density(double frequency) {
    double omega = 2.0 * M_PI * frequency;
    return omega * omega * EPSILON_0 * M_ELECTRON / (E_CHARGE * E_CHARGE);
}

/* ================================================================
 * L5: Thomson Scattering
 * ================================================================ */

double thomson_cross_section(double theta) {
    double r_e = E_CHARGE * E_CHARGE
                 / (4.0 * M_PI * EPSILON_0 * M_ELECTRON * C_LIGHT * C_LIGHT);
    return r_e * r_e * sin(theta) * sin(theta);
}

double thomson_scattered_power(double P_incident, double ne,
                               double dsigma_dOmega,
                               double solid_angle, double L_scattering) {
    return P_incident * ne * dsigma_dOmega * solid_angle * L_scattering;
}

double thomson_te_from_width(double delta_lambda, double lambda_i,
                              double theta) {
    double sin_half = sin(theta / 2.0);
    if (sin_half <= 0.0) return 0.0;
    double v_th_sq = (delta_lambda * C_LIGHT / (2.0 * lambda_i * sin_half));
    v_th_sq *= v_th_sq;
    return M_ELECTRON * v_th_sq / (2.0 * K_B);
}

/* ================================================================
 * L5: Spectroscopy
 * ================================================================ */

double doppler_broadening_fwhm(double lambda_0, double T, double m) {
    double factor = 8.0 * K_B * T * log(2.0) / (m * C_LIGHT * C_LIGHT);
    return lambda_0 * sqrt(factor);
}

double stark_broadening_fwhm(double ne) {
    return 2.5e-14 * pow(ne, 2.0/3.0);
}

double zeeman_splitting(double lambda_0, double B) {
    return E_CHARGE * B * lambda_0 * lambda_0
           / (4.0 * M_PI * M_ELECTRON * C_LIGHT);
}

/* ================================================================
 * L7: Space and Fusion Plasma Diagnostics
 * ================================================================ */

void solar_wind_params(PlasmaRegime *regime) {
    double ne = 5.0e6;
    double Te = 1.0e5;
    classify_plasma(ne, Te, regime);
}

double magnetopause_distance(double B_0, double R_E,
                             double rho_sw, double v_sw) {
    double p_dyn = rho_sw * v_sw * v_sw;
    if (p_dyn <= 0.0) return INFINITY;
    double p_mag = B_0 * B_0 / (2.0 * MU_0);
    return R_E * pow(p_mag / p_dyn, 1.0/6.0);
}

double fusion_power_density(double n_D, double n_T, double sigmav) {
    double E_fusion = 17.6e6 * E_CHARGE;
    return n_D * n_T * sigmav * E_fusion;
}

double dt_fusion_reactivity(double Ti_keV) {
    if (Ti_keV <= 0.0) return 0.0;
    double theta = Ti_keV / (1.0 - Ti_keV * 1.72e-3
                             + Ti_keV * Ti_keV * 6.38e-6);
    double xi = 6.661e-3 / pow(theta, 5.0);
    xi = pow(xi, 1.0/6.0);
    return 1.1e-24 * pow(theta, 0.4) * exp(-3.0 * xi);
}

int lawson_criterion(double n, double T_keV, double tau_E,
                     double threshold) {
    return (n * T_keV * tau_E >= threshold) ? 1 : 0;
}
