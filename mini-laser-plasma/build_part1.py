import os

base = r"F:\nano-everything\mini-pure-physics\8. mini-plasma-physics\mini-laser-plasma"

def W(rel, content):
    full = os.path.join(base, rel)
    os.makedirs(os.path.dirname(full), exist_ok=True)
    with open(full, "w", encoding="utf-8") as f:
        f.write(content)
    lines = content.count("\n")
    print(f"  {rel}: {lines} lines")

print("=== Building mini-laser-plasma ===")

# ========== plasma_params.h ==========
W("include/plasma_params.h", """#ifndef PLASMA_PARAMS_H
#define PLASMA_PARAMS_H
#include <math.h>

/*
 * plasma_params.h -- Fundamental plasma parameter calculations
 * References: Goldston & Rutherford (1995), Chen (2016), Kruer (1988)
 * Knowledge: L1 Definitions, L3 Math Structures, L4 Fundamental Laws
 */

typedef struct {
    double ne, Te, Ti, Z, A, B_field, lam_laser, I_laser;
} PlasmaState;

typedef struct {
    double wp, wpi, nc, lambda_D, N_D;
    double v_the, v_thi, cs, nu_ei, ln_Lambda;
    double beta, a0, ne_over_nc;
    double lambda_plasma, skin_depth, omega_laser;
} PlasmaDerived;

/* L1: Core plasma parameters (each function = one knowledge point) */
double plasma_frequency(double ne);
double ion_plasma_frequency(double ne, double Z, double A);
double critical_density(double lambda_m);
double debye_length(double ne, double Te_eV);
double debye_sphere_particles(double ne, double lambda_D);
double electron_thermal_velocity(double Te_eV);
double ion_thermal_velocity(double Ti_eV, double A);
double ion_sound_speed(double Te_eV, double Z, double A);
double skin_depth(double wp);

/* L2: Collision parameters */
double coulomb_logarithm(double ne, double Te_eV, double Z);
double electron_ion_collision_frequency(double ne, double Te_eV, double Z, double lnL);
double collision_mean_free_path(double v_the, double nu_ei);
double spitzer_resistivity(double Te_eV, double Z, double lnL);

/* L3: Dielectric response */
double cold_plasma_permittivity(double omega, double wp);
double collisional_permittivity(double omega, double wp, double nu, double *eps_imag);
double warm_plasma_permittivity(double omega, double k, double wp, double v_the);
double plasma_refractive_index(double ne, double nc);
double laser_wavelength_in_plasma(double lambda_vac, double ne, double nc);
double group_velocity(double ne, double nc);

/* L4: Pressure and beta */
double plasma_beta(double ne, double Te_eV, double Ti_eV, double B_field);
double electron_pressure(double ne, double Te_eV);
double radiation_pressure(double intensity, int reflecting);
int compute_all_derived(const PlasmaState *ps, PlasmaDerived *pd);

#endif
""")

print("Done with plasma_params.h")
