/**
 * fusion_core.c — Fusion Plasma Core Physics Implementation
 *
 * Every function implements an independent physics knowledge point.
 * Refs: Goldston & Rutherford (1995), Wesson 4th ed (2011),
 *       Freidberg (2007), Bosch & Hale (1992),
 *       ITER Physics Basis (1999)
 *
 * L2: Debye shielding, plasma frequency, thermal velocities,
 *     gyroradii, cyclotron frequencies, Alfven/ion sound speeds,
 *     plasma beta, safety factor, aspect ratio, Greenwald limit,
 *     collision frequencies, Spitzer resistivity
 * L3: MHD wave speeds, Reynolds/Lundquist numbers,
 *     energy equilibration, trapped particle fraction
 * L4: Lawson criterion, triple product, ignition condition,
 *     bremsstrahlung/cyclotron/line radiation, breakeven,
 *     fusion power density (D-T, D-D, D-He3, p-B11),
 *     ohmic/NBI/alpha heating, fusion gain
 * L5: Bosch-Hale cross-section parameterizations (D-T, D-D, D-He3),
 *     Maxwellian reactivity Gauss-Laguerre integration
 * L6: Tritium breeding ratio, neutron wall loading,
 *     energy multiplication, alpha slowing-down
 */

#include "include/fusion_plasma.h"
#include <math.h>
#include <assert.h>

/* ================================================================
 * L2: Debye Shielding and Plasma Frequency — 6 functions
 * ================================================================ */

/**
 * debye_length — Debye shielding length [m]
 * lambda_D = sqrt(epsilon0 * kB * Te / (ne * e^2))
 * For Te=10keV, ne=1e20 m^-3: lambda_D ~ 7.4e-5 m
 */
double debye_length(double Te_eV, double ne) {
    assert(ne > 0.0);
    assert(Te_eV >= 0.0);
    double kT = Te_eV * E_CHARGE;  /* kT in Joules */
    return sqrt(EPSILON0 * kT / (ne * E_CHARGE * E_CHARGE));
}

/**
 * plasma_frequency — Electron plasma frequency [rad/s]
 * omega_pe = sqrt(ne * e^2 / (epsilon0 * me))
 * For ne=1e20 m^-3: f_pe ~ 90 GHz
 */
double plasma_frequency(double ne) {
    assert(ne >= 0.0);
    return sqrt(ne * E_CHARGE * E_CHARGE / (EPSILON0 * M_E));
}

/**
 * ion_plasma_frequency — Ion plasma frequency [rad/s]
 * omega_pi = sqrt(ni * Z^2 * e^2 / (epsilon0 * mi))
 * omega_pi/omega_pe = sqrt(Z^2 * me/mi)
 */
double ion_plasma_frequency(double ni, double mi_kg, double Z) {
    assert(ni >= 0.0);
    assert(mi_kg > 0.0);
    assert(Z >= 1.0);
    return sqrt(ni * Z * Z * E_CHARGE * E_CHARGE / (EPSILON0 * mi_kg));
}

/**
 * coulomb_logarithm — Coulomb logarithm ln(Lambda)
 * lnLambda = ln(4*pi * n * lambda_D^3)
 * Fusion plasmas: lnLambda ~ 15-20
 * Clamped to [2, 30] for physical validity.
 */
double coulomb_logarithm(double Te_eV, double ne) {
    double lambda_D = debye_length(Te_eV, ne);
    double N_D = (4.0 * M_PI / 3.0) * ne * lambda_D * lambda_D * lambda_D;
    if (N_D < 1.0) N_D = 1.0;
    double lnL = log(N_D);
    if (lnL < 2.0) lnL = 2.0;
    if (lnL > 30.0) lnL = 30.0;
    return lnL;
}

/**
 * plasma_parameter — Plasma coupling parameter Lambda
 * Lambda = 4*pi * n * lambda_D^3
 * Lambda >> 1: weakly coupled (kinetic theory valid)
 * Lambda <= 1: strongly coupled (liquid-like behavior)
 */
double plasma_parameter(double Te_eV, double ne) {
    double lambda_D = debye_length(Te_eV, ne);
    return 4.0 * M_PI * ne * lambda_D * lambda_D * lambda_D;
}

/* ================================================================
 * L2: Thermal Velocities and Gyroradii — 8 functions
 * ================================================================ */

/**
 * electron_thermal_velocity — Electron thermal speed [m/s]
 * v_the = sqrt(2 * kB * Te / me)
 * For Te=10keV: v_the ~ 5.9e7 m/s (~0.2c)
 */
double electron_thermal_velocity(double Te_eV) {
    double Te_J = Te_eV * E_CHARGE;
    return sqrt(2.0 * Te_J / M_E);
}

/**
 * ion_thermal_velocity — Ion thermal speed [m/s]
 * v_thi = sqrt(2 * kB * Ti / mi)
 * For D at 10keV: v_thi ~ 9.8e5 m/s
 */
double ion_thermal_velocity(double Ti_eV, double mi_kg) {
    double Ti_J = Ti_eV * E_CHARGE;
    return sqrt(2.0 * Ti_J / mi_kg);
}

/**
 * ion_gyroradius — Ion Larmor radius [m]
 * rho_i = sqrt(2*mi*kB*Ti) / (Z*e*B)
 * For D at 10keV, B=5T: rho_i ~ 4.1e-3 m
 */
double ion_gyroradius(double Ti_eV, double mi_kg, double B, double Z) {
    assert(B > 0.0);
    assert(Z >= 1.0);
    double v_thi = ion_thermal_velocity(Ti_eV, mi_kg);
    return mi_kg * v_thi / (Z * E_CHARGE * B);
}

/**
 * electron_gyroradius — Electron Larmor radius [m]
 * rho_e = sqrt(2*me*kB*Te) / (e*B)
 * For Te=10keV, B=5T: rho_e ~ 6.7e-5 m
 */
double electron_gyroradius(double Te_eV, double B) {
    assert(B > 0.0);
    double v_the = electron_thermal_velocity(Te_eV);
    return M_E * v_the / (E_CHARGE * B);
}

/**
 * ion_cyclotron_frequency — Ion cyclotron frequency [rad/s]
 * Omega_ci = Z * e * B / mi
 * For D at B=5T: f_ci ~ 38 MHz
 */
double ion_cyclotron_frequency(double B, double mi_kg, double Z) {
    assert(mi_kg > 0.0);
    return Z * E_CHARGE * B / mi_kg;
}

/**
 * electron_cyclotron_frequency — Electron cyclotron frequency [rad/s]
 * Omega_ce = e * B / me
 * For B=5T: f_ce ~ 140 GHz
 */
double electron_cyclotron_frequency(double B) {
    return E_CHARGE * B / M_E;
}

/**
 * ion_sound_speed — Ion sound speed [m/s]
 * c_s = sqrt((Z*Te + gamma_i*Ti) / mi)
 * For D plasma with Te=10keV: c_s ~ 9.8e5 m/s
 */
double ion_sound_speed(double Te_eV, double Ti_eV, double mi_kg, double Z) {
    double Te_J = Te_eV * E_CHARGE;
    double Ti_J = Ti_eV * E_CHARGE;
    double gamma_i = 3.0;
    return sqrt((Z * Te_J + gamma_i * Ti_J) / mi_kg);
}

/**
 * alfven_speed — Alfven speed [m/s]
 * v_A = B / sqrt(mu0 * ni * mi)
 * For D plasma, B=5T, ni=1e20: v_A ~ 8.7e6 m/s
 */
double alfven_speed(double B, double ni, double mi_kg) {
    assert(ni > 0.0);
    assert(mi_kg > 0.0);
    return B / sqrt(MU0 * ni * mi_kg);
}

/* ================================================================
 * L2: Plasma Pressure, Beta, Geometry — 10 functions
 * ================================================================ */

double plasma_beta(double n, double T_eV, double B) {
    assert(B > 0.0);
    double kT = T_eV * E_CHARGE;  /* kT in Joules */
    return 2.0 * MU0 * n * kT / (B * B);
}

double plasma_beta_poloidal(double n, double T_eV, double B_pol) {
    assert(B_pol > 0.0);
    double kT2 = T_eV * E_CHARGE;  /* kT in Joules */
    return 2.0 * MU0 * n * kT2 / (B_pol * B_pol);
}

double safety_factor_cylindrical(double a, double R, double B_phi, double B_theta) {
    assert(B_theta > 0.0);
    return (a / R) * (B_phi / B_theta);
}

double safety_factor_from_current(double a, double R, double B_phi, double Ip) {
    assert(Ip > 0.0);
    return (2.0 * M_PI * a * a * B_phi) / (MU0 * R * Ip);
}

double aspect_ratio(double R, double a) {
    assert(a > 0.0);
    return R / a;
}

double plasma_volume(double R, double a, double kappa) {
    return 2.0 * M_PI * M_PI * kappa * a * a * R;
}

double plasma_surface_area(double R, double a, double kappa) {
    double h = (kappa - 1.0) / (kappa + 1.0);
    double h2 = h * h;
    double C_elliptic = M_PI * (kappa + 1.0) * a *
                        (1.0 + 3.0*h2 / (10.0 + sqrt(4.0 - 3.0*h2)));
    return 2.0 * M_PI * R * C_elliptic;
}

double plasma_cross_section_area(double a, double kappa) {
    return M_PI * kappa * a * a;
}

double greenwald_density_limit(double Ip, double a) {
    assert(a > 0.0);
    /* n_G [10^20 m^-3] = Ip [MA] / (pi * a^2 [m^2]) */
    /* n_G [m^-3] = (Ip [A] / 1e6) / (pi * a^2) * 1e20 = Ip / (pi*a^2) * 1e14 */
    return Ip / (M_PI * a * a) * 1e14;
}

double greenwald_fraction(double ne, double n_GW) {
    assert(n_GW > 0.0);
    return ne / n_GW;
}

/* ================================================================
 * L3: Collision Frequencies and Transport Basics — 7 functions
 * ================================================================ */

double collisionality_nu_star(double R, double nu_ei, double v_the, double epsilon) {
    assert(v_the > 0.0);
    assert(epsilon > 0.0);
    double eps_pow = epsilon * sqrt(epsilon);
    return (R * nu_ei) / (v_the * eps_pow);
}

double trapped_particle_fraction(double epsilon) {
    assert(epsilon >= 0.0 && epsilon < 1.0);
    return sqrt(2.0 * epsilon);
}

double electron_electron_collision_freq(double ne, double Te_eV, double lnLambda) {
    double v_the = electron_thermal_velocity(Te_eV);
    double e4 = E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE;
    double eps02 = EPSILON0 * EPSILON0;
    return ne * e4 * lnLambda /
           (4.0 * M_PI * eps02 * M_E * M_E * v_the * v_the * v_the);
}

double electron_ion_collision_freq(double ni, double Z, double Te_eV, double lnLambda) {
    double v_the = electron_thermal_velocity(Te_eV);
    double e4 = E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE;
    double eps02 = EPSILON0 * EPSILON0;
    return ni * Z * Z * e4 * lnLambda /
           (4.0 * M_PI * eps02 * M_E * M_E * v_the * v_the * v_the);
}

double ion_ion_collision_freq(double ni, double Z, double Ti_eV, double mi_kg, double lnLambda) {
    double v_thi = ion_thermal_velocity(Ti_eV, mi_kg);
    double e4 = E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE;
    double eps02 = EPSILON0 * EPSILON0;
    double Z4 = Z * Z * Z * Z;
    return ni * Z4 * e4 * lnLambda /
           (4.0 * M_PI * eps02 * mi_kg * mi_kg * v_thi * v_thi * v_thi);
}

double energy_equilibration_time(double ne, double Z, double Te_eV, double Ti_eV,
                                  double mi_kg, double lnLambda) {
    double Te_J = Te_eV * E_CHARGE;
    double Ti_J = Ti_eV * E_CHARGE;
    double term = Te_J / M_E + Ti_J / mi_kg;  /* Te_J and Ti_J are already kT in J */
    double v_rel_eff = pow(term, 1.5);
    double e4 = E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE;
    double eps02 = EPSILON0 * EPSILON0;
    double num = 3.0 * pow(M_PI, 1.5) * eps02 * M_E * mi_kg;
    double den = ne * Z * Z * e4 * lnLambda;
    return num * v_rel_eff / den;
}

double spitzer_resistivity(double Z, double Te_eV, double lnLambda) {
    double Te_J = Te_eV * E_CHARGE;
    double num = M_PI * Z * E_CHARGE * E_CHARGE * sqrt(M_E) * lnLambda;
    double den = 16.0 * M_PI * M_PI * EPSILON0 * EPSILON0 *
                 pow(Te_J, 1.5);  /* Te_J is already kT in J */
    double gamma_E = 0.58 + 0.20 * (Z - 1.0) / (Z + 1.0);
    double correction = 0.295 * gamma_E;
    return correction * num / den;
}

/* ================================================================
 * L3: MHD Wave Speeds and Dimensionless Numbers — 6 functions
 * ================================================================ */

double magnetic_pressure(double B) {
    return B * B / (2.0 * MU0);
}

double thermal_pressure(double n, double T_eV) {
    return n * T_eV * E_CHARGE;  /* T_eV * e = kT in Joules */
}

double magnetic_reynolds_number(double L, double v, double eta) {
    assert(eta > 0.0);
    return MU0 * L * v / eta;
}

double lundquist_number(double L, double B, double ni, double mi_kg, double eta) {
    double v_A = alfven_speed(B, ni, mi_kg);
    return MU0 * L * v_A / eta;
}

double fast_magnetosonic_speed(double c_s, double v_A) {
    return sqrt(c_s * c_s + v_A * v_A);
}

double slow_magnetosonic_speed(double c_s, double v_A) {
    double den = sqrt(c_s * c_s + v_A * v_A);
    if (den < 1e-30) return 0.0;
    return c_s * v_A / den;
}

/* ================================================================
 * L4: Lawson Criterion and Ignition — 4 functions
 * ================================================================ */

double lawson_criterion(double n, double T, double tau_E) {
    assert(n >= 0.0);
    assert(T >= 0.0);
    assert(tau_E >= 0.0);
    return n * T * tau_E;
}

double triple_product(double n, double T, double tau_E) {
    return lawson_criterion(n, T, tau_E);
}

int ignition_condition(double n, double T, double tau_E) {
    double lp = lawson_criterion(n, T, tau_E);
    return (lp >= 3.0) ? 1 : 0;
}

int breakeven_condition(double P_fusion, double P_aux) {
    if (P_aux <= 0.0) return 0;
    return (P_fusion >= P_aux) ? 1 : 0;
}

/* ================================================================
 * L4: Fusion Power Density — 4 functions
 * ================================================================ */

double fusion_power_density_dt(double nD, double nT, double sigma_v, double Efus) {
    if (Efus <= 0.0) Efus = E_FUSION_DT_J;
    return nD * nT * sigma_v * Efus;
}

double fusion_power_density_dd(double nD, double sigma_v) {
    double E_DD_J = E_FUSION_DD * E_CHARGE;
    return 0.5 * nD * nD * sigma_v * E_DD_J;
}

double fusion_power_density_dhe3(double nD, double nHe3, double sigma_v) {
    double E_DHe3_J = E_FUSION_DHe3 * E_CHARGE;
    return nD * nHe3 * sigma_v * E_DHe3_J;
}

double fusion_power_density_pb11(double np, double nB11, double sigma_v) {
    double E_pB11 = 8.68e6 * E_CHARGE;
    return np * nB11 * sigma_v * E_pB11;
}

/* ================================================================
 * L4: Radiation Losses — 6 functions
 * ================================================================ */

double bremsstrahlung_power_density(double ne, double Zeff, double Te_eV) {
    const double C_B = 5.35e-37;
    return C_B * Zeff * ne * ne * sqrt(Te_eV);
}

double recombination_radiation_power(double ne, double Zeff, double Te_eV) {
    assert(Te_eV > 0.0);
    const double C_R = 5.0e-38;
    double Z3 = Zeff * Zeff * Zeff;
    return C_R * Z3 * ne * ne / sqrt(Te_eV);
}

double line_radiation_power(double ne, double nz, double L_z) {
    return L_z * ne * nz;
}

double cyclotron_power_density(double ne, double Te_eV, double B) {
    double omega_ce = electron_cyclotron_frequency(B);
    double Te_J = Te_eV * E_CHARGE;
    double e2 = E_CHARGE * E_CHARGE;
    double num = e2 * omega_ce * omega_ce;
    double den = 6.0 * M_PI * EPSILON0 * C_LIGHT * C_LIGHT * C_LIGHT;
    return (num / den) * ne * Te_J;  /* Te_J is already kT */
}

double bremsstrahlung_relativistic_correction(double Te_eV) {
    double Te_J = Te_eV * E_CHARGE;
    double me_c2 = M_E * C_LIGHT * C_LIGHT;
    double x = Te_J / me_c2;
    return 1.0 + 0.7936 * x + 1.874 * x * x;
}

double total_radiation_power(double ne, double Zeff, double Te_eV, double B,
                              double nz, double L_z) {
    double P_brem = bremsstrahlung_power_density(ne, Zeff, Te_eV);
    double P_cyc = cyclotron_power_density(ne, Te_eV, B);
    double P_line = line_radiation_power(ne, nz, L_z);
    return P_brem + P_cyc + P_line;
}

/* ================================================================
 * L4: Heating Power — 5 functions
 * ================================================================ */

double ohmic_heating_power(double eta, double Ip, double a, double kappa) {
    double A_cs = plasma_cross_section_area(a, kappa);
    if (A_cs <= 0.0) return 0.0;
    double j = Ip / A_cs;
    return eta * j * j;
}

double nbi_heating_power_deposited(double P_inj, double f_shine, double f_orbit) {
    assert(f_shine >= 0.0 && f_shine <= 1.0);
    assert(f_orbit >= 0.0 && f_orbit <= 1.0);
    assert(f_shine + f_orbit <= 1.0);
    return P_inj * (1.0 - f_shine - f_orbit);
}

double alpha_power_density(double nD, double nT, double sigma_v) {
    double E_alpha_J = E_ALPHA * E_CHARGE;
    return nD * nT * sigma_v * E_alpha_J;
}

double alpha_slowing_down_time(double ne, double Te_eV, double lnLambda) {
    double Te_J = Te_eV * E_CHARGE;
    double num = 3.0 * sqrt(M_PI) * EPSILON0 * EPSILON0 *
                 M_He4 * Te_J * sqrt(Te_J);  /* Te_J is already kT */
    double den = ne * E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE *
                 sqrt(M_E) * lnLambda;
    return num / den;
}

double alpha_birth_velocity(void) {
    double E_alpha_J = E_ALPHA * E_CHARGE;
    return sqrt(2.0 * E_alpha_J / M_He4);
}

double alpha_critical_energy(double Te_eV) {
    return 14.8 * 4.0 * Te_eV;
}

/* ================================================================
 * L6: Fusion Gain and Confinement — 6 functions
 * ================================================================ */

double fusion_gain_Q(double P_fusion, double P_aux) {
    if (P_aux <= 0.0) {
        return (P_fusion > 0.0) ? 1e10 : 0.0;
    }
    return P_fusion / P_aux;
}

double engineering_fusion_gain(double P_elec, double P_wall) {
    if (P_wall <= 0.0) return 0.0;
    return P_elec / P_wall;
}

double alpha_heating_fraction(double n, double T, double tau_E) {
    double sigma_v_dt = bosch_hale_sigma_v_dt(T);
    double n_SI = n * 1e20;
    double P_alpha_core = alpha_power_density(n_SI/2.0, n_SI/2.0, sigma_v_dt);
    double T_J = T * 1e3 * E_CHARGE;
    double P_loss = 3.0 * n_SI * K_BOLTZMANN * T_J / tau_E;
    if (P_loss <= 0.0) return 0.0;
    return P_alpha_core / P_loss;
}

double energy_confinement_time(double W_plasma, double P_loss) {
    if (P_loss <= 0.0) return 1e10;
    return W_plasma / P_loss;
}

double plasma_stored_energy(const PlasmaParameters *plasma) {
    double kTe = plasma->Te * E_CHARGE;  /* kT_e in Joules */
    double kTi = plasma->Ti * E_CHARGE;  /* kT_i in Joules */
    double energy_density = 1.5 * (plasma->ne * kTe +
                                    plasma->ni * kTi);
    return energy_density * plasma->V_p;
}

double energy_confinement_iter89p(double Ip_MA, double B, double n_20,
                                   double P_MW, double R, double a, double kappa) {
    assert(P_MW > 0.0);
    return 0.048 * pow(Ip_MA, 0.85) * pow(B, 0.2) * pow(n_20, 0.1) *
           pow(P_MW, -0.5) * pow(R, 1.2) * pow(a, 0.3) * pow(kappa, 0.5);
}

double energy_confinement_ipb98y2(double Ip_MA, double B, double n_20,
                                   double P_MW, double R, double a, double kappa) {
    assert(P_MW > 0.0);
    double eps = a / R;
    return 0.0562 * pow(Ip_MA, 0.93) * pow(B, 0.15) * pow(n_20, 0.41) *
           pow(P_MW, -0.69) * pow(R, 1.97) * pow(eps, 0.58) * pow(kappa, 0.78);
}

/* ================================================================
 * L5: Bosch-Hale Cross-Sections — 4 functions
 * ================================================================ */

double bosch_hale_sigma_v_dt(double Ti_keV) {
    if (Ti_keV <= 0.1) return 0.0;
    if (Ti_keV > 200.0) Ti_keV = 200.0;
    const double B_G = 34.3827;
    const double mRc2 = 1124656.0;
    const double C1 = 1.17302e-15;  /* cm^3/s -> m^3/s */
    const double C2 = 1.51361e-2;
    const double C3 = 7.51886e-2;
    const double C4 = 4.60643e-3;
    const double C5 = 1.35000e-2;
    const double C6 = -1.06750e-4;
    const double C7 = 1.36600e-5;
    double Ti3 = Ti_keV * Ti_keV * Ti_keV;
    double num = 1.0 - Ti_keV * (C2 + Ti_keV * (C4 + Ti_keV * C6));
    double den = 1.0 + Ti_keV * (C3 + Ti_keV * (C5 + Ti_keV * C7));
    if (den <= 0.0 || num <= 0.0) return 0.0;  /* formula breakdown guard */
    double theta = Ti_keV / (num / den);
    if (theta <= 0.0) return 0.0;
    double xi = pow(B_G * B_G / (4.0 * theta), 1.0/3.0);
    return C1 * theta * sqrt(xi / (mRc2 * Ti3)) * exp(-3.0 * xi);
}

double bosch_hale_sigma_v_dd(double Ti_keV) {
    if (Ti_keV <= 0.5) return 0.0;
    if (Ti_keV > 200.0) Ti_keV = 200.0;
    const double B_G = 31.3970;
    const double mRc2 = 937814.0;
    const double C1 = 5.43360e-18;  /* cm^3/s -> m^3/s */
    const double C2 = 5.85778e-3;
    const double C3 = 7.68222e-3;
    const double C5 = -2.96400e-6;
    double Ti3 = Ti_keV * Ti_keV * Ti_keV;
    double num = 1.0 - Ti_keV * C2;
    double den = 1.0 + Ti_keV * (C3 + Ti_keV * C5);
    if (den <= 0.0 || num <= 0.0) return 0.0;
    double theta = Ti_keV / (num / den);
    if (theta <= 0.0) return 0.0;
    double xi = pow(B_G * B_G / (4.0 * theta), 1.0/3.0);
    return C1 * theta * sqrt(xi / (mRc2 * Ti3)) * exp(-3.0 * xi);
}

double bosch_hale_sigma_v_dhe3(double Ti_keV) {
    if (Ti_keV <= 1.0) return 0.0;
    if (Ti_keV > 200.0) Ti_keV = 200.0;
    const double B_G = 68.7508;
    const double mRc2 = 1124572.0;
    const double C1 = 5.51036e-16;  /* cm^3/s -> m^3/s */
    const double C2 = 6.41918e-3;
    const double C3 = -2.02896e-3;
    const double C4 = -1.91080e-5;
    const double C5 = 1.35776e-4;
    double Ti3 = Ti_keV * Ti_keV * Ti_keV;
    double num = 1.0 - Ti_keV * (C2 + Ti_keV * C4);
    double den = 1.0 + Ti_keV * (C3 + Ti_keV * C5);
    if (den <= 0.0 || num <= 0.0) return 0.0;
    double theta = Ti_keV / (num / den);
    if (theta <= 0.0) return 0.0;
    double xi = pow(B_G * B_G / (4.0 * theta), 1.0/3.0);
    return C1 * theta * sqrt(xi / (mRc2 * Ti3)) * exp(-3.0 * xi);
}

double fusion_reactivity_maxwellian(double Ti_keV,
                                     double (*cross_section)(double E_J),
                                     double reduced_mass_kg) {
    if (Ti_keV <= 0.0 || cross_section == NULL || reduced_mass_kg <= 0.0) {
        return 0.0;
    }
    double T_J = Ti_keV * 1e3 * E_CHARGE;
    double kBT = K_BOLTZMANN * T_J;
    /* 4-point Gauss-Laguerre quadrature */
    double nodes[4] = {
        0.3225476896193923, 1.7457611011583466,
        4.536620296921127,  9.395070912301133
    };
    double weights[4] = {
        0.6031541043416336, 0.3574186924377997,
        0.0388879085150054, 0.0005392947055614
    };
    double integral = 0.0;
    for (int i = 0; i < 4; i++) {
        double x = nodes[i];
        double E_J = x * kBT;
        double sigma = cross_section(E_J);
        double v = sqrt(2.0 * E_J / reduced_mass_kg);
        double f_M = (2.0 / sqrt(M_PI)) * sqrt(E_J) /
                     pow(kBT, 1.5) * exp(-x);
        integral += weights[i] * sigma * v * f_M * exp(x);
    }
    return integral;
}

/* ================================================================
 * L6: Tritium and Neutron Engineering — 4 functions
 * ================================================================ */

double tritium_burning_fraction(double n_T_burned, double n_T_injected) {
    if (n_T_injected <= 0.0) return 0.0;
    return n_T_burned / n_T_injected;
}

double tritium_breeding_ratio_needed(double f_burn, double eff_extraction) {
    if (f_burn <= 0.0 || eff_extraction <= 0.0) return 1e10;
    return 1.0 / (f_burn * eff_extraction);
}

double neutron_wall_loading(double P_fusion, double fraction_neutron, double S_wall) {
    if (S_wall <= 0.0) return 0.0;
    return fraction_neutron * P_fusion / S_wall;
}

double neutron_energy_multiplication(double E_total_out, double E_neutron_in) {
    if (E_neutron_in <= 0.0) return 1.0;
    return E_total_out / E_neutron_in;
}