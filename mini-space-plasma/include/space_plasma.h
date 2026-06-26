/**
 * space_plasma.h — Core Space Plasma Physics Definitions
 *
 * Reference: Kivelson & Russell "Introduction to Space Physics" (1995)
 *            Goldston & Rutherford "Introduction to Plasma Physics" (1995)
 *            MIT 22.611 / 22.616 Space Plasma Physics
 *
 * Covers L1 (Definitions), L2 (Core Concepts), L4 (Fundamental Laws)
 *
 * All units: SI unless otherwise noted.
 */

#ifndef SPACE_PLASMA_H
#define SPACE_PLASMA_H

#include <stddef.h>
#include <math.h>

/* Ensure M_PI is defined (C11 compliance) */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*──────────────────────────────────────────────────────────────────────
 * L1: Fundamental Physical Constants (CODATA 2018)
 *──────────────────────────────────────────────────────────────────────*/

#define SP_KB       (1.380649e-23)    /* Boltzmann constant [J/K]      */
#define SP_HBAR     (1.054571817e-34)  /* Reduced Planck constant [J·s] */
#define SP_C        (299792458.0)      /* Speed of light [m/s]          */
#define SP_ME       (9.1093837015e-31) /* Electron mass [kg]            */
#define SP_MP       (1.67262192369e-27)/* Proton mass [kg]              */
#define SP_EC       (1.602176634e-19)  /* Elementary charge [C]         */
#define SP_EPS0     (8.8541878128e-12) /* Vacuum permittivity [F/m]     */
#define SP_MU0      (1.25663706212e-6) /* Vacuum permeability [N/A²]    */
#define SP_G        (6.67430e-11)      /* Gravitational constant [N·m²/kg²] */
#define SP_NA       (6.02214076e23)    /* Avogadro constant [1/mol]     */
#define SP_SIGMA_SB (5.670374419e-23)  /* Stefan-Boltzmann [W/(m²·K⁴)] */
#define SP_RYDBERG  (13.605693122994)  /* Rydberg energy [eV]          */
#define SP_AU       (1.495978707e11)   /* Astronomical Unit [m]         */
#define SP_RE       (6.371e6)          /* Earth radius [m]              */
#define SP_RSUN     (6.957e8)          /* Solar radius [m]              */
#define SP_MEARTH   (5.9722e24)        /* Earth mass [kg]               */
#define SP_MSUN     (1.98847e30)       /* Solar mass [kg]               */
#define SP_B0_EARTH (3.12e-5)          /* Earth equatorial B [T]       */
#define SP_MU_EARTH (7.78e22)          /* Earth dipole moment [A·m²]   */

/*──────────────────────────────────────────────────────────────────────
 * L1: Core Plasma Parameter Definitions (struct-based)
 *──────────────────────────────────────────────────────────────────────*/

/** @brief Plasma state at a single spatial point */
typedef struct {
    double n_e;          /* Electron number density [m⁻³]             */
    double n_i;          /* Ion number density [m⁻³]                  */
    double T_e;          /* Electron temperature [K]                  */
    double T_i;          /* Ion temperature [K]                       */
    double B[3];         /* Magnetic field vector [T]                 */
    double v[3];         /* Bulk flow velocity [m/s]                  */
    double p_e;          /* Electron pressure [Pa] (derived)          */
    double p_i;          /* Ion pressure [Pa] (derived)               */
    double Z;            /* Ion charge state (dimensionless)          */
    int    is_magnetized;/* Flag: ω_ci τ_ii ≫ 1                       */
} plasma_state_t;

/** @brief Single particle orbit parameters */
typedef struct {
    double r_L;          /* Larmor radius [m]                         */
    double omega_c;      /* Gyrofrequency [rad/s]                     */
    double v_perp;       /* Perpendicular velocity [m/s]              */
    double v_par;        /* Parallel velocity [m/s]                   */
    double mu_m;         /* Magnetic moment (first adiabatic inv)     */
    double pitch_angle;  /* Pitch angle [rad]                         */
    double guiding_center[3]; /* Guiding center position [m]          */
} particle_orbit_t;

/** @brief MHD wave modes characterization */
typedef struct {
    double v_A;          /* Alfven speed [m/s]                        */
    double c_s;          /* Sound speed [m/s]                         */
    double v_fast;       /* Fast magnetosonic speed [m/s]             */
    double v_slow;       /* Slow magnetosonic speed [m/s]             */
    double beta;         /* Plasma beta (2μ₀p/B²)                    */
    double theta;        /* Angle between k and B [rad]              */
} mhd_wave_t;

/** @brief Debye screening parameters */
typedef struct {
    double lambda_D;     /* Electron Debye length [m]                 */
    double lambda_Di;    /* Ion Debye length [m]                      */
    double N_D;          /* Number of particles in Debye sphere       */
    double omega_pe;     /* Electron plasma frequency [rad/s]         */
    double omega_pi;     /* Ion plasma frequency [rad/s]              */
    double tau;          /* Temperature ratio T_e/T_i                 */
} debye_params_t;

/** @brief Magnetospheric region parameters */
typedef struct {
    double R_mp;         /* Magnetopause standoff distance [R_E]      */
    double R_ss;         /* Subsolar distance [R_E]                   */
    double B_dipole[3];  /* Dipole field at location [T]              */
    double L_shell;      /* McIlwain L-parameter                      */
    double pressure_sw;  /* Solar wind dynamic pressure [nPa]        */
    double Dst;          /* Dst index analogue [nT]                   */
} magnetosphere_t;

/** @brief Solar wind state */
typedef struct {
    double v_sw;         /* Solar wind speed [m/s]                    */
    double n_sw;         /* Solar wind density [m⁻³]                  */
    double T_sw;         /* Solar wind temperature [K]                */
    double B_imf[3];     /* Interplanetary magnetic field [T]         */
    double Ma;           /* Alfven Mach number                        */
    double Ms;           /* Sonic Mach number                         */
    double M_ms;         /* Magnetosonic Mach number                  */
    double beta_sw;      /* Solar wind beta                           */
} solar_wind_t;

/** @brief Magnetic reconnection geometry */
typedef struct {
    double B_in;         /* Inflow magnetic field [T]                 */
    double v_in;         /* Inflow velocity [m/s]                     */
    double v_out;        /* Outflow velocity [m/s]                    */
    double delta;        /* Current sheet half-thickness [m]          */
    double L;            /* Current sheet length [m]                  */
    double eta;          /* Resistivity [Ω·m]                         */
    double S;            /* Lundquist number                          */
    double E_rec;        /* Reconnection electric field [V/m]         */
} reconnection_t;

/** @brief Kinetic distribution function (discretized) */
typedef struct {
    double *f;           /* Distribution function values              */
    double  v_min;       /* Minimum velocity [m/s]                    */
    double  v_max;       /* Maximum velocity [m/s]                    */
    double  dv;          /* Velocity resolution [m/s]                 */
    size_t  n_v;         /* Number of velocity bins                   */
    double  n0;          /* Reference density [m⁻³]                   */
    double  T;           /* Temperature [K]                           */
    double  v_th;        /* Thermal velocity [m/s]                    */
} distribution_t;

/** @brief Plasma resistivity models */
typedef enum {
    RESISTIVITY_SPITZER,       /* Spitzer resistivity η ∝ T⁻³ᐟ²      */
    RESISTIVITY_CHROMOSPHERIC, /* Enhanced resistivity for low temp   */
    RESISTIVITY_ANOMALOUS,     /* Anomalous from wave-particle        */
    RESISTIVITY_CONSTANT       /* User-specified constant             */
} resistivity_model_t;

/** @brief Space plasma region identifier */
typedef enum {
    REGION_SOLAR_WIND,
    REGION_BOW_SHOCK,
    REGION_MAGNETOSHEATH,
    REGION_MAGNETOPAUSE,
    REGION_MAGNETOSPHERE,
    REGION_PLASMASPHERE,
    REGION_RING_CURRENT,
    REGION_TAIL_LOBE,
    REGION_PLASMA_SHEET,
    REGION_RADIATION_BELT,
    REGION_IONOSPHERE,
    REGION_HELIOPAUSE
} plasma_region_t;

/**
 * @brief Velocity-space coordinate for distribution function
 *
 * Used in kinetic theory calculations.
 */
typedef struct {
    double v_par;         /* Parallel velocity [m/s]                  */
    double v_perp;        /* Perpendicular velocity magnitude [m/s]   */
    double mu;            /* cos(pitch_angle)                         */
    double E;             /* Kinetic energy [J]                       */
} velocity_coord_t;

/**
 * @brief Dispersion relation result
 */
typedef struct {
    double omega;         /* Real frequency [rad/s]                   */
    double gamma;         /* Growth/damping rate [rad/s]              */
    double k[3];          /* Wave vector [rad/m]                      */
    double v_phase;       /* Phase speed [m/s]                        */
    double v_group;       /* Group speed [m/s]                        */
    int    mode;          /* Wave mode identifier                     */
} dispersion_t;

/*──────────────────────────────────────────────────────────────────────
 * Macro helpers
 *──────────────────────────────────────────────────────────────────────*/

/** @brief Compute scalar pressure from density and temperature */
#define SP_PRESSURE(n, T) ((n) * SP_KB * (T))

/** @brief Compute thermal velocity (most probable) */
#define SP_VTHERM(T, m) (sqrt(2.0 * SP_KB * (T) / (m)))

/** @brief Compute squared magnitude of 3-vector */
#define SP_DOT3(a, b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

/** @brief Compute magnitude of 3-vector */
#define SP_MAG3(a) (sqrt(SP_DOT3(a, a)))

/** @brief Compute cross product z = x × y (3D) */
#define SP_CROSS3(x, y, z) do { \
    (z)[0] = (x)[1]*(y)[2] - (x)[2]*(y)[1]; \
    (z)[1] = (x)[2]*(y)[0] - (x)[0]*(y)[2]; \
    (z)[2] = (x)[0]*(y)[1] - (x)[1]*(y)[0]; \
} while(0)

/** @brief Safe division: returns 0.0 if denominator is 0 */
static inline double sp_safe_div(double num, double den) {
    return (fabs(den) < 1e-300) ? 0.0 : (num / den);
}

#endif /* SPACE_PLASMA_H */
