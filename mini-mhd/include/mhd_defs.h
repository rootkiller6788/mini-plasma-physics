/*
 * mhd_defs.h -- Core Definitions for Magnetohydrodynamics
 *
 * Reference: Goldston & Rutherford (1995), Freidberg "Ideal MHD" (2014)
 *            Biskamp "MHD Turbulence" (2003), Priest "MHD of the Sun" (2014)
 * Course: MIT 22.611, Princeton AST 554, Cambridge Part III MHD
 *
 * Knowledge:
 *   L1 -- Physical constants (CODATA 2018), MHD state structures (8-field model)
 *   L2 -- Core plasma parameters: beta, Mach, Reynolds, Lundquist, Hartmann...
 *   L3 -- Vector algebra operations, conservative/primitive variable transforms
 */

#ifndef MHD_DEFS_H
#define MHD_DEFS_H

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * L1 -- Physical Constants (CODATA 2018)
 * ================================================================ */

#define MHD_MU0         1.25663706212e-6   /* mu_0 vacuum permeability [N/A^2] */
#define MHD_EPS0        8.8541878128e-12   /* epsilon_0 vacuum permittivity [F/m] */
#define MHD_C           2.99792458e8       /* c speed of light [m/s] */
#define MHD_KB          1.380649e-23       /* k_B Boltzmann constant [J/K] */
#define MHD_E           1.602176634e-19    /* e elementary charge [C] */
#define MHD_ME          9.1093837015e-31   /* m_e electron mass [kg] */
#define MHD_MP          1.67262192369e-27  /* m_p proton mass [kg] */
#define MHD_NA          6.02214076e23      /* N_A Avogadro number [1/mol] */
#define MHD_R_GAS       8.314462618        /* R universal gas constant [J/(mol K)] */
#define MHD_SIGMA_SB    5.670374419e-8     /* sigma Stefan-Boltzmann [W/(m^2 K^4)] */
#define MHD_HBAR        1.054571817e-34    /* hbar reduced Planck [J s] */

#define MHD_MU0_INV      7.95774715459e5   /* 1/mu_0 */
#define MHD_HALF_MU0_INV 3.97887357730e5   /* 1/(2*mu_0) */

/* ================================================================
 * L1 -- MHD State Structures
 *
 * The 8-field model (rho, v, B, p) is the standard representation
 * of compressible ideal/resistive MHD. The constraint div(B)=0
 * must be maintained at all times.
 * ================================================================ */

/* Primitive variables: (rho, v, B, p) -- 8 fields */
typedef struct {
    double rho;          /* mass density [kg/m^3] */
    double vx, vy, vz;   /* fluid velocity [m/s] */
    double Bx, By, Bz;   /* magnetic field [T] */
    double p;            /* thermal pressure [Pa] */
} MHDState;

/* Conserved variables: (D, M, B, E) -- for finite-volume schemes */
typedef struct {
    double D;            /* rho -- mass density */
    double Mx, My, Mz;   /* rho*v -- momentum density */
    double Bx, By, Bz;   /* B -- magnetic field */
    double E;            /* total energy density = p/(g-1) + rho*v^2/2 + B^2/(2*mu_0) */
} MHDConserved;

/* MHD flux vector across a face with normal n */
typedef struct {
    double D_flux;                   /* mass flux: rho * v_n */
    double Mx_flux, My_flux, Mz_flux; /* momentum flux tensor row */
    double Bx_flux, By_flux, Bz_flux; /* induction flux */
    double E_flux;                   /* energy flux */
} MHDFlux;

/* Characteristic speeds of the 8x8 MHD flux Jacobian */
typedef struct {
    double c_fast;    /* fast magnetosonic */
    double c_slow;    /* slow magnetosonic */
    double c_alfven;  /* Alfven (intermediate) */
    double c_sound;   /* adiabatic sound */
} MHDWaveSpeeds;

/* Resistive MHD source terms (finite eta) */
typedef struct {
    double ohmic_heating;               /* eta * J^2 [W/m^3] */
    double Bx_diss, By_diss, Bz_diss;   /* eta * Laplacian(B) */
    double E_diss;                      /* total energy dissipation */
} MHDResistiveSource;

/* Supported coordinate systems */
typedef enum {
    MHD_COORD_CARTESIAN   = 0,
    MHD_COORD_CYLINDRICAL = 1,
    MHD_COORD_SPHERICAL   = 2,
    MHD_COORD_TOROIDAL    = 3
} MHDCoordinateSystem;

/* Computational domain geometry */
typedef struct {
    MHDCoordinateSystem coord;
    double Lx, Ly, Lz;       /* domain lengths [m] */
    int    nx, ny, nz;       /* grid cell counts */
    double dx, dy, dz;       /* cell sizes (uniform grid) */
    double x0, y0, z0;       /* domain origin */
} MHDGeometry;

/* Full set of dimensionless plasma parameters */
typedef struct {
    double beta;                 /* plasma beta = 2*mu_0*p/B^2 */
    double mach_alfven;          /* Alfven Mach M_A = v/v_A */
    double mach_sonic;           /* sonic Mach M_s = v/c_s */
    double magnetic_reynolds;    /* R_m = v*L/eta */
    double lundquist;            /* S = v_A*L/eta */
    double magnetic_prandtl;     /* P_m = nu/eta */
    double hartmann;             /* Ha = B*L/sqrt(mu_0*rho*nu*eta) */
    double ion_inertial_length;  /* d_i = c/omega_pi [m] */
    double ion_larmor_radius;    /* rho_i = m_i*v_th/(e*B) [m] */
    double debye_length;         /* lambda_D [m] */
    double elsasser;             /* Lambda = B^2/(mu_0*rho*eta*Omega) */
    double electron_plasma_beta; /* beta_e = 2*mu_0*n*k_B*T_e/B^2 */
    double ion_plasma_beta;      /* beta_i = 2*mu_0*n*k_B*T_i/B^2 */
} PlasmaParameters;

/* ================================================================
 * L2 -- Core MHD Physical Quantities (inline functions)
 *
 * Each function encodes exactly one physical concept
 * with a single, clean mathematical formula.
 * ================================================================ */

/*
 * mhd_sound_speed -- Adiabatic sound speed
 * c_s = sqrt(gamma * p / rho)
 * Speed of compressional (acoustic) waves in the absence of B.
 * For monatomic ideal gas (gamma=5/3): c_s = sqrt(5p/3*rho).
 * Ref: Landau & Lifshitz, Fluid Mechanics, sec.63
 */
static inline double mhd_sound_speed(double p, double rho, double gamma) {
    return sqrt(gamma * p / rho);
}

/*
 * mhd_alfven_speed -- Alfven speed
 * v_A = B / sqrt(mu_0 * rho)
 * Characteristic speed of transverse perturbations along B.
 * Named after Hannes Alfven (Nobel Prize 1970).
 * Ref: Alfven, Nature 150, 405 (1942)
 */
static inline double mhd_alfven_speed(double B_mag, double rho) {
    return B_mag / sqrt(MHD_MU0 * rho);
}

/*
 * mhd_magnetic_pressure -- Magnetic pressure
 * p_mag = B^2 / (2*mu_0)
 * Isotropic part of the Maxwell stress tensor: T_ij^mag = -p_mag * delta_ij + B_i*B_j/mu_0.
 * The diagonal part B^2/(2*mu_0) acts as an isotropic pressure.
 */
static inline double mhd_magnetic_pressure(double B_mag) {
    return (B_mag * B_mag) / (2.0 * MHD_MU0);
}

/*
 * mhd_magnetic_tension -- Magnetic tension
 * T_B = B^2 / mu_0
 * The restoring force due to field line curvature.
 * The Lorentz force decomposes as JxB = -grad_perp(p_mag) + (B^2/mu_0) * n_hat/R_c
 * where n_hat/R_c is the field line curvature vector.
 */
static inline double mhd_magnetic_tension(double B_mag) {
    return (B_mag * B_mag) / MHD_MU0;
}

/*
 * mhd_plasma_beta -- Plasma beta
 * beta = 2 * mu_0 * p / B^2 = p / p_mag
 * The most important dimensionless MHD parameter.
 * beta < 1: magnetically dominated (tokamak ~ 0.01-0.1, corona ~ 0.01)
 * beta > 1: kinetically dominated (solar wind, magnetotail)
 * beta = 1: equipartition
 * High beta is desirable for economical fusion reactors.
 */
static inline double mhd_plasma_beta(double p, double B_mag) {
    if (B_mag < 1e-40) return INFINITY;
    return 2.0 * MHD_MU0 * p / (B_mag * B_mag);
}

/*
 * mhd_mach_alfven -- Alfven Mach number
 * M_A = v / v_A
 * M_A < 1: sub-Alfvenic, M_A > 1: super-Alfvenic (shocks)
 * Solar wind at 1 AU: M_A ~ 10
 */
static inline double mhd_mach_alfven(double v_mag, double B_mag, double rho) {
    double va = mhd_alfven_speed(B_mag, rho);
    if (va < 1e-40) return INFINITY;
    return v_mag / va;
}

/*
 * mhd_mach_sonic -- Sonic Mach number
 * M_s = v / c_s
 * M_s < 1: subsonic (incompressible limit)
 * M_s > 1: supersonic (shocks form)
 */
static inline double mhd_mach_sonic(double v_mag, double p, double rho, double gamma) {
    double cs = mhd_sound_speed(p, rho, gamma);
    if (cs < 1e-40) return INFINITY;
    return v_mag / cs;
}

/*
 * mhd_magnetic_reynolds -- Magnetic Reynolds number
 * R_m = v * L / eta   (where eta = 1/(mu_0*sigma) is magnetic diffusivity)
 * Determines if B is frozen into the fluid.
 * R_m >> 1: advection dominates, frozen-in flux (Alfven theorem)
 * R_m << 1: diffusion dominates, B decays
 * Astrophysics: R_m ~ 10^8-10^20, ideal MHD is excellent.
 * Lab liquid metals: R_m ~ 0.1-10, resistive effects important.
 */
static inline double mhd_magnetic_reynolds(double v_char, double L_char, double eta) {
    if (eta < 1e-40) return INFINITY;
    return v_char * L_char / eta;
}

/*
 * mhd_lundquist -- Lundquist number
 * S = v_A * L / eta
 * Central parameter in magnetic reconnection theory.
 * Sweet-Parker reconnection rate ~ S^{-1/2}.
 * Tokamaks: S ~ 10^7-10^9, Solar corona: S ~ 10^12
 */
static inline double mhd_lundquist(double B_mag, double rho, double L_char, double eta) {
    double va = mhd_alfven_speed(B_mag, rho);
    if (eta < 1e-40) return INFINITY;
    return va * L_char / eta;
}

/*
 * mhd_magnetic_prandtl -- Magnetic Prandtl number
 * P_m = nu / eta = R_m / Re
 * Ratio of viscous to magnetic diffusivity.
 * Liquid metals (Earth core): P_m ~ 10^{-5}
 * Fully ionized plasmas: P_m ~ 10-100
 * ISM: P_m ~ 10^{14}
 * Controls relative roles of kinetic and magnetic energy cascades.
 * Ref: Biskamp, MHD Turbulence (2003)
 */
static inline double mhd_magnetic_prandtl(double nu, double eta) {
    if (eta < 1e-40) return INFINITY;
    return nu / eta;
}

/*
 * mhd_hartmann -- Hartmann number
 * Ha = B * L * sqrt(sigma / (rho * nu)) = B*L / sqrt(mu_0 * rho * nu * eta)
 * Ha^2 = Lorentz force / viscous force
 * Ha >> 1: flow flattened by magnetic damping (Hartmann layer)
 * Classical MHD duct flow problem.
 */
static inline double mhd_hartmann(double B_mag, double L_char,
                                   double rho, double nu, double eta) {
    double denom = MHD_MU0 * rho * nu * eta;
    if (denom < 1e-60) return INFINITY;
    return B_mag * L_char * sqrt(1.0 / denom);
}

/*
 * mhd_ion_inertial_length -- Ion skin depth
 * d_i = c / omega_pi = c * sqrt(epsilon_0 * m_i / (n * e^2))
 * Scale below which electrons and ions decouple (Hall MHD regime).
 * Solar wind: d_i ~ 100 km, tokamaks: d_i ~ 1 cm.
 * MHD valid when d_i/L << 1.
 */
static inline double mhd_ion_inertial_length(double n, double m_i) {
    double omega_pi_sq = n * MHD_E * MHD_E / (MHD_EPS0 * m_i);
    if (omega_pi_sq < 1e-60) return INFINITY;
    return MHD_C / sqrt(omega_pi_sq);
}

/*
 * mhd_ion_larmor_radius -- Ion thermal Larmor radius
 * rho_i = v_{th,perp} / omega_ci = m_i * v_th / (e * B)
 * For MHD validity: rho_i / L << 1 (small gyroradius approximation)
 * Tokamaks: rho_i ~ 1 mm, magnetotail: rho_i ~ 100 km
 */
static inline double mhd_ion_larmor_radius(double T_i, double B_mag, double m_i) {
    double v_th = sqrt(2.0 * MHD_KB * T_i / m_i);
    if (B_mag < 1e-40) return INFINITY;
    return m_i * v_th / (MHD_E * B_mag);
}

/*
 * mhd_debye_length -- Debye screening length
 * lambda_D = sqrt(epsilon_0 * k_B * T_e / (n * e^2))
 * Scale of electrostatic screening. Quasineutrality requires L >> lambda_D.
 * Fusion plasmas: lambda_D ~ 10^{-5} m, solar wind: ~ 10 m
 */
static inline double mhd_debye_length(double T_e, double n) {
    double denom = n * MHD_E * MHD_E;
    if (denom < 1e-60) return INFINITY;
    return sqrt(MHD_EPS0 * MHD_KB * T_e / denom);
}

/*
 * mhd_ion_cyclotron_frequency -- Ion cyclotron frequency
 * omega_ci = e * B / m_i
 * For protons in B=5 T: omega_ci ~ 4.8e8 rad/s.
 */
static inline double mhd_ion_cyclotron_frequency(double B_mag, double m_i) {
    return MHD_E * B_mag / m_i;
}

/*
 * mhd_plasma_frequency -- Electron plasma frequency
 * omega_pe = sqrt(n * e^2 / (epsilon_0 * m_e))
 * For tokamak n~10^20: omega_pe ~ 5.6e11 rad/s.
 * MHD assumes omega << omega_pe.
 */
static inline double mhd_plasma_frequency(double n) {
    return sqrt(n * MHD_E * MHD_E / (MHD_EPS0 * MHD_ME));
}

/*
 * mhd_elsasser -- Elsasser number
 * Lambda = B^2 / (mu_0 * rho * eta * Omega)
 * Ratio of Lorentz to Coriolis force in rotating MHD.
 * Key parameter in geodynamo: Earth ~ 1-10, Jupiter ~ 100.
 * Lambda >> 1: field dominates rotation.
 */
static inline double mhd_elsasser(double B_mag, double rho, double eta, double Omega) {
    double denom = MHD_MU0 * rho * eta * Omega;
    if (denom < 1e-40) return INFINITY;
    return (B_mag * B_mag) / denom;
}

/* ================================================================
 * L3 -- Vector Algebra Operations (inline)
 * ================================================================ */

static inline double mhd_vector_magnitude(double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
}

static inline double mhd_dot_product(double ax, double ay, double az,
                                      double bx, double by, double bz) {
    return ax*bx + ay*by + az*bz;
}

static inline void mhd_cross_product(double ax, double ay, double az,
                                      double bx, double by, double bz,
                                      double *cx, double *cy, double *cz) {
    *cx = ay*bz - az*by;
    *cy = az*bx - ax*bz;
    *cz = ax*by - ay*bx;
}

/* ================================================================
 * L1/L3 -- Energy Density Functions (inline)
 * ================================================================ */

static inline double mhd_kinetic_energy_density(double rho, double v_mag) {
    return 0.5 * rho * v_mag * v_mag;
}

static inline double mhd_magnetic_energy_density(double B_mag) {
    return (B_mag * B_mag) / (2.0 * MHD_MU0);
}

static inline double mhd_internal_energy_density(double p, double gamma) {
    return p / (gamma - 1.0);
}

static inline double mhd_total_energy_density(double rho, double v_mag,
                                               double B_mag, double p,
                                               double gamma) {
    return mhd_internal_energy_density(p, gamma)
         + mhd_kinetic_energy_density(rho, v_mag)
         + mhd_magnetic_energy_density(B_mag);
}

/*
 * mhd_current_from_curlB -- J = curl(B)/mu_0
 * Ampere's law without displacement current (MHD approximation, v << c).
 */
static inline void mhd_current_from_curlB(double cbx, double cby, double cbz,
                                           double *Jx, double *Jy, double *Jz) {
    double inv_mu0 = 1.0 / MHD_MU0;
    *Jx = cbx * inv_mu0;
    *Jy = cby * inv_mu0;
    *Jz = cbz * inv_mu0;
}

/*
 * mhd_lorentz_force -- J x B force density
 * The magnetic force in the momentum equation.
 * Equivalent forms: JxB = curl(B)xB/mu_0 = -grad(B^2/2mu_0) + (B.grad)B/mu_0
 */
static inline void mhd_lorentz_force(double Jx, double Jy, double Jz,
                                      double Bx, double By, double Bz,
                                      double *Fx, double *Fy, double *Fz) {
    mhd_cross_product(Jx, Jy, Jz, Bx, By, Bz, Fx, Fy, Fz);
}

/*
 * mhd_joule_heating -- Ohmic dissipation rate
 * Q_eta = eta * J^2 = eta * |curl(B)|^2 / mu_0^2
 */
static inline double mhd_joule_heating(double J_mag, double eta) {
    return eta * J_mag * J_mag;
}

/* ================================================================
 * Non-inline function declarations (implemented in src/)
 * ================================================================ */

void mhd_compute_all_parameters(const MHDState *state, double gamma,
                                 double L_char, double eta, double nu,
                                 double m_i, double n, double T_e, double T_i,
                                 double Omega, PlasmaParameters *params);

void mhd_primitive_to_conserved(const MHDState *prim, double gamma,
                                 MHDConserved *cons);

void mhd_conserved_to_primitive(const MHDConserved *cons, double gamma,
                                 MHDState *prim);

void mhd_flux_compute(const MHDState *state, double nx, double ny, double nz,
                       double gamma, MHDFlux *flux);

double mhd_max_wavespeed(const MHDState *state, double nx, double ny, double nz,
                          double gamma);

/* Utility functions */
void mhd_state_init_uniform(MHDState *state,
                             double rho, double vx, double vy, double vz,
                             double Bx, double By, double Bz, double p);
void mhd_geometry_init(MHDGeometry *geom, MHDCoordinateSystem coord,
                        double Lx, double Ly, double Lz, int nx, int ny, int nz);
void mhd_state_print(const MHDState *state);
void mhd_params_print(const PlasmaParameters *params);

#endif /* MHD_DEFS_H */
