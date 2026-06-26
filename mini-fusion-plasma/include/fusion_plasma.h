/**
 * fusion_plasma.h — Fusion Plasma Physics Core Header
 *
 * Refs: Goldston & Rutherford "Introduction to Plasma Physics" (1995)
 *       Wesson "Tokamaks" 4th ed. (2011)
 *       Freidberg "Plasma Physics and Fusion Energy" (2007)
 *       Miyamoto "Plasma Physics and Controlled Nuclear Fusion" (2005)
 *       MIT 22.611/22.615, Princeton PHY 563
 *       Cambridge Part III Plasma Physics / Fusion
 *       Oxford CMT Plasma Theory, ETH 402-0851 Plasma Physics
 *
 * Nine-Level Knowledge Coverage:
 *   L1: PlasmaParameters, FusionFuelType, DebyeSphere, LawsonCriteria,
 *       FusionReactionRate, EnergyBalance, TokamakGeometry,
 *       ConfinementScaling, GradShafranovSolution, HeatingSystem,
 *       TritiumBreedingRatio, DivertorParameters (>5 independent structs)
 *   L2: Debye shielding, Coulomb collision, fusion cross-section,
 *       Lawson criterion, quasi-neutrality, plasma oscillation
 *   L3: Coulomb logarithm, flux coordinates (psi,theta,zeta),
 *       Fokker-Planck operator, Grad-Shafranov operator
 *   L4: Lawson criterion, triple product, ignition condition,
 *       bremsstrahlung, cyclotron radiation, Spitzer resistivity
 *   L5: Bosch-Hale cross-section, Maxwellian reactivity,
 *       IPB98(y,2) scaling, Gauss-Laguerre quadrature
 *   L6: D-T/D-D/D-He3 fusion, tokamak, ITER, DEMO, JET, SPARC
 *   L7: ITER power balance, tritium breeding, neutron wall loading,
 *       divertor heat flux, fusion gain Q
 *   L8: Burning plasma, alpha slowing-down, H-mode pedestal,
 *       neoclassical tearing modes, sawtooth oscillations
 *   L9: Spherical torus, stellarator optimization, liquid metal walls
 */

#ifndef MINI_FUSION_PLASMA_H
#define MINI_FUSION_PLASMA_H

#include <stddef.h>
#include <math.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * L1: Physical Constants (CODATA 2018)
 * ================================================================ */
#define K_BOLTZMANN   1.380649e-23
#define E_CHARGE      1.602176634e-19
#define M_E           9.1093837015e-31
#define M_P           1.67262192369e-27
#define M_D           3.3435837724e-27
#define M_T           5.0073567446e-27
#define M_He3         5.0064118880e-27
#define M_He4         6.6446573357e-27
#define M_N           1.67492749804e-27
#define EPSILON0      8.8541878128e-12
#define MU0           1.25663706212e-6
#define C_LIGHT       299792458.0
#define H_PLANCK      6.62607015e-34
#define HBAR          1.054571817e-34
#define N_AVOGADRO    6.02214076e23
#define SIGMA_SB      5.670374419e-8

/* Fusion-specific derived constants */
#define E_FUSION_DT   17.59e6
#define E_FUSION_DD   3.27e6
#define E_FUSION_DHe3 18.35e6
#define E_ALPHA       3.52e6
#define E_NEUTRON_DT  14.07e6
#define E_FUSION_DT_J (E_FUSION_DT * E_CHARGE)

/* ================================================================
 * L1: Fusion Fuel Types
 * ================================================================ */
typedef enum {
    FUEL_DT   = 0,
    FUEL_DD   = 1,
    FUEL_DHe3 = 2,
    FUEL_pB11 = 3,
    FUEL_CAT_DD = 4,
    FUEL_TT   = 5,
    FUEL_pLi6 = 6
} FusionFuelType;

/* ================================================================
 * L1: Core Physics Type Definitions
 * ================================================================ */

/** PlasmaParameters - plasma state variables
 *  ne,ni [m^-3]; Te,Ti [eV]; B [T]; Ip [A]; R,a [m] */
typedef struct {
    double ne;
    double ni;
    double Te;
    double Ti;
    double B;
    double Ip;
    double R;
    double a;
    double elongation;
    double triangularity;
    double Zeff;
    double V_p;
} PlasmaParameters;

/** DebyeSphere - Debye screening quantities
 *  lambda_D = sqrt(epsilon0 * kB * Te / (ne * e^2))
 *  N_D = (4*pi/3) * ne * lambda_D^3
 *  Ref: Goldston & Rutherford sec 2.4 */
typedef struct {
    double lambda_D;
    double N_D;
    double plasma_param;
    double coulomb_log;
} DebyeSphere;

/** LawsonCriteria - Lawson criterion parameters
 *  triple_product = n * T * tau_E >= 3e21 m^-3 keV s (ignition)
 *  Ref: J.D. Lawson, Proc. Phys. Soc. B 70, 6 (1957) */
typedef struct {
    double density;
    double temperature;
    double tau_E;
    double triple_product;
    double required_min;
    int    satisfied;
} LawsonCriteria;

/** FusionReactionRate - reaction rate data
 *  R = n1*n2*<sigma v> / (1+delta_12)
 *  P_fus = R * E_fus */
typedef struct {
    double n1;
    double n2;
    double sigma_v;
    double reaction_rate;
    double power_density;
} FusionReactionRate;

/** EnergyBalance - plasma power balance
 *  dW/dt = P_fus + P_aux - P_rad - P_transport
 *  Q = P_fus / P_aux */
typedef struct {
    double P_fusion;
    double P_alpha;
    double P_aux;
    double P_radiation;
    double P_transport;
    double P_ohmic;
    double Q;
    int    ignited;
} EnergyBalance;

/** TokamakGeometry - tokamak shape parameters
 *  A = R/a (aspect ratio), kappa = elongation */
typedef struct {
    double R;
    double a;
    double kappa;
    double delta;
    double aspect_ratio;
    double volume;
    double surface_area;
} TokamakGeometry;

/* ================================================================
 * L3: Mathematical Structures
 * ================================================================ */

/** 3D Vector (Euclidean R^3) */
typedef struct { double x, y, z; } Vec3;

/** 3x3 Matrix (diffusion tensor, stress tensor) */
typedef struct { double m[3][3]; } Mat3x3;

/* Vector operations */
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline Vec3 vec3_scale(double s, Vec3 v) {
    return (Vec3){s * v.x, s * v.y, s * v.z};
}
static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
static inline double vec3_norm(Vec3 v) {
    return sqrt(vec3_dot(v, v));
}
static inline Vec3 vec3_normalize(Vec3 v) {
    double n = vec3_norm(v);
    return (n < 1e-30) ? v : vec3_scale(1.0/n, v);
}
static inline Vec3 vec3_zero(void) { return (Vec3){0.0,0.0,0.0}; }

static inline Mat3x3 mat3_identity(void) {
    Mat3x3 I = {{{0}}};
    I.m[0][0] = I.m[1][1] = I.m[2][2] = 1.0;
    return I;
}
static inline Vec3 mat3_mul_vec(Mat3x3 M, Vec3 v) {
    return (Vec3){
        M.m[0][0]*v.x+M.m[0][1]*v.y+M.m[0][2]*v.z,
        M.m[1][0]*v.x+M.m[1][1]*v.y+M.m[1][2]*v.z,
        M.m[2][0]*v.x+M.m[2][1]*v.y+M.m[2][2]*v.z
    };
}

/* ================================================================
 * L3: Flux Coordinates and Derived Types
 * ================================================================ */

/** FluxSurface - flux coordinates (psi, theta, zeta)
 *  psi=poloidal flux, theta=poloidal angle, zeta=toroidal angle
 *  Ref: Wesson sec 3, Freidberg sec 6 */
typedef struct {
    double psi;
    double psi_norm;
    double theta;
    double zeta;
    double q;
    double shear;
} FluxSurface;

/** ConfinementScaling - energy confinement scaling law parameters
 *  tau_E = H * tau_IPB98(y,2)
 *  Ref: ITER Physics Basis, Nucl. Fusion 39, 2137 (1999) */
typedef struct {
    double tau_E_predicted;
    double H_factor;
    double tau_E_empirical;
    double rms_error;
} ConfinementScaling;

/** GradShafranovSolution - GS equation solution
 *  Delta* psi = -mu0*R^2*dp/dpsi - F*dF/dpsi
 *  Ref: Grad & Rubin (1958), Shafranov (1966) */
typedef struct {
    double psi_axis;
    double psi_boundary;
    double beta_p;
    double li;
    double shafranov_shift;
} GradShafranovSolution;

/** HeatingSystem - heating power components
 *  P_total = P_nbi + P_icrf + P_ecrh + P_lh + P_alpha + P_ohmic */
typedef struct {
    double P_nbi;
    double P_icrf;
    double P_ecrh;
    double P_lh;
    double P_alpha;
    double P_ohmic;
    double P_total;
} HeatingSystem;

/** TritiumBreedingRatio - tritium self-sufficiency metrics
 *  TBR = tritium produced / tritium consumed (>1 for self-sufficiency) */
typedef struct {
    double tbr;
    double tbr_required;
    double doubling_time;
    int    self_sufficient;
} TritiumBreedingRatio;

/** DivertorParameters - divertor physics
 *  q_peak: peak heat flux [MW/m^2]
 *  detachment_degree: 0=attached, 1=fully detached */
typedef struct {
    double q_peak;
    double q_parallel;
    double T_e_div;
    double n_e_div;
    double detachment_degree;
} DivertorParameters;

/* ================================================================
 * L2: Fundamental Plasma Parameter Functions
 * ================================================================ */

double debye_length(double Te_eV, double ne);
double plasma_frequency(double ne);
double ion_plasma_frequency(double ni, double mi_kg, double Z);
double coulomb_logarithm(double Te_eV, double ne);
double plasma_parameter(double Te_eV, double ne);
double electron_thermal_velocity(double Te_eV);
double ion_thermal_velocity(double Ti_eV, double mi_kg);
double ion_gyroradius(double Ti_eV, double mi_kg, double B, double Z);
double electron_gyroradius(double Te_eV, double B);
double ion_cyclotron_frequency(double B, double mi_kg, double Z);
double electron_cyclotron_frequency(double B);
double ion_sound_speed(double Te_eV, double Ti_eV, double mi_kg, double Z);
double alfven_speed(double B, double ni, double mi_kg);
double plasma_beta(double n, double T_eV, double B);
double plasma_beta_poloidal(double n, double T_eV, double B_pol);
double safety_factor_cylindrical(double a, double R, double B_phi, double B_theta);
double safety_factor_from_current(double a, double R, double B_phi, double Ip);
double aspect_ratio(double R, double a);
double plasma_volume(double R, double a, double kappa);
double plasma_surface_area(double R, double a, double kappa);
double plasma_cross_section_area(double a, double kappa);
double greenwald_density_limit(double Ip, double a);
double greenwald_fraction(double ne, double n_GW);
double collisionality_nu_star(double R, double nu_ei, double v_the, double epsilon);
double trapped_particle_fraction(double epsilon);
double electron_electron_collision_freq(double ne, double Te_eV, double lnLambda);
double electron_ion_collision_freq(double ni, double Z, double Te_eV, double lnLambda);
double ion_ion_collision_freq(double ni, double Z, double Ti_eV, double mi_kg, double lnLambda);
double energy_equilibration_time(double ne, double Z, double Te_eV, double Ti_eV, double mi_kg, double lnLambda);
double spitzer_resistivity(double Z, double Te_eV, double lnLambda);
double magnetic_pressure(double B);
double thermal_pressure(double n, double T_eV);
double magnetic_reynolds_number(double L, double v, double eta);
double lundquist_number(double L, double B, double ni, double mi_kg, double eta);
double fast_magnetosonic_speed(double c_s, double v_A);
double slow_magnetosonic_speed(double c_s, double v_A);

/* ================================================================
 * L4: Fundamental Fusion Laws
 * ================================================================ */

double lawson_criterion(double n, double T, double tau_E);
double triple_product(double n, double T, double tau_E);
int    ignition_condition(double n, double T, double tau_E);
int    breakeven_condition(double P_fusion, double P_aux);
double fusion_power_density_dt(double nD, double nT, double sigma_v, double Efus);
double fusion_power_density_dd(double nD, double sigma_v);
double fusion_power_density_dhe3(double nD, double nHe3, double sigma_v);
double fusion_power_density_pb11(double np, double nB11, double sigma_v);
double bremsstrahlung_power_density(double ne, double Zeff, double Te_eV);
double recombination_radiation_power(double ne, double Zeff, double Te_eV);
double line_radiation_power(double ne, double nz, double L_z);
double cyclotron_power_density(double ne, double Te_eV, double B);
double bremsstrahlung_relativistic_correction(double Te_eV);
double total_radiation_power(double ne, double Zeff, double Te_eV, double B, double nz, double L_z);
double ohmic_heating_power(double eta, double Ip, double a, double kappa);
double nbi_heating_power_deposited(double P_inj, double f_shine, double f_orbit);
double alpha_power_density(double nD, double nT, double sigma_v);
double alpha_slowing_down_time(double ne, double Te_eV, double lnLambda);
double alpha_birth_velocity(void);
double alpha_critical_energy(double Te_eV);
double fusion_gain_Q(double P_fusion, double P_aux);
double engineering_fusion_gain(double P_elec, double P_wall);
double alpha_heating_fraction(double n, double T, double tau_E);
double energy_confinement_time(double W_plasma, double P_loss);
double plasma_stored_energy(const PlasmaParameters *plasma);
double energy_confinement_iter89p(double Ip_MA, double B, double n_20, double P_MW, double R, double a, double kappa);
double energy_confinement_ipb98y2(double Ip_MA, double B, double n_20, double P_MW, double R, double a, double kappa);

/* ================================================================
 * L5: Bosch-Hale Cross-Section Parameterization
 * ================================================================ */
double bosch_hale_sigma_v_dt(double Ti_keV);
double bosch_hale_sigma_v_dd(double Ti_keV);
double bosch_hale_sigma_v_dhe3(double Ti_keV);
double fusion_reactivity_maxwellian(double Ti_keV, double (*cross_section)(double E_J), double reduced_mass_kg);

/* ================================================================
 * L6: Canonical Systems - Tokamak Operating Scenarios
 * ================================================================ */

double troyon_beta_limit(double beta_percent, double a, double B, double Ip_MA);
double normalized_beta(double beta_percent, double a, double B, double Ip_MA);
int    kruskal_shafranov_limit(double q_edge);
double bootstrap_current_fraction(double epsilon, double beta_p);
double bootstrap_current_density(double epsilon, double B_pol, double dp_dr);
double pfirsch_schluter_factor(double q);
double banana_regime_diffusivity(double epsilon, double q, double rho_i, double nu_ii);
double plateau_regime_diffusivity(double rho_i, double v_thi, double R, double q);
double pfirsch_schluter_diffusivity(double q, double rho_i, double nu_ii);
double neoclassical_ion_heat_diffusivity(double epsilon, double q, double rho_i, double nu_ii, double v_thi, double R);
double bohm_diffusion_coefficient(double Te_eV, double B);
double gyrobohm_diffusion_coefficient(double Te_eV, double B, double rho_i, double a);
double h_mode_power_threshold(double n_20, double B, double S);
double pedestal_pressure(double Ip_MA, double a, double B);
double pedestal_width(double beta_p_ped, double a);
double tearing_mode_delta_prime(double psi_prime_plus, double psi_prime_minus, double psi_resonant);
double ballooning_alpha(double q, double R, double dp_dr, double B);
double mercier_criterion(double q, double shear, double magnetic_well, double pressure_gradient);
double magnetic_shear(double r, double q, double dq_dr);
double magnetic_well_depth(double V_prime_0, double V_prime_psi);
double sawtooth_period(double a, double eta, double v_A);
double sawtooth_mixing_radius(double r_q1);
double disruption_density_limit(double Ip, double a, double margin);
double ntm_threshold_island_width(double rho_i, double beta_p);
double halo_current_fraction(double I_halo, double I_p0);
double tritium_burning_fraction(double n_T_burned, double n_T_injected);
double tritium_breeding_ratio_needed(double f_burn, double eff_extraction);
double neutron_wall_loading(double P_fusion, double fraction_neutron, double S_wall);
double neutron_energy_multiplication(double E_total_out, double E_neutron_in);
double divertor_heat_flux(double P_sol, double alpha_deg, double R, double lambda_q, double f_exp);
double power_exhaust_fraction(double P_rad, double P_heat);

/* L6: Operating point functions */
void iter_operating_point(PlasmaParameters *p);
void demo_operating_point(PlasmaParameters *p);
void sparc_operating_point(PlasmaParameters *p);
void jet_dt_record(PlasmaParameters *p);
void iter_power_balance(const PlasmaParameters *p, EnergyBalance *eb);

#endif /* MINI_FUSION_PLASMA_H */