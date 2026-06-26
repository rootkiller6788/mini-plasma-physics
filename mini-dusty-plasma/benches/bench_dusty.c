/**
 * @file  bench_dusty.c
 * @brief Performance benchmarks for dusty plasma computations.
 *
 * Benchmarks key functions to establish baseline performance:
 * - Debye length computation (O(1))
 * - OML current evaluation (O(1))
 * - Floating potential solver (O(iter))
 * - DAW dispersion sweep (O(n_k))
 * - Yukawa lattice sum (O(N^2))
 * - Structure factor computation (O(n_k * n_bins))
 */

#include "dusty_plasma.h"
#include "dusty_charging.h"
#include "dusty_waves.h"
#include "dusty_crystal.h"
#include "dusty_transport.h"
#include <math.h>
#include <stdio.h>
#include <time.h>

#define BENCH_ITER 100000
#define BENCH(name, expr) do { \
    clock_t start = clock(); \
    for (int _i = 0; _i < BENCH_ITER; _i++) { volatile double _r = (expr); (void)_r; } \
    clock_t end = clock(); \
    double ms = 1000.0 * (end - start) / CLOCKS_PER_SEC; \
    printf("  %-30s %8.3f ms (%d iters, %.1f ns/op)\n", name, ms, BENCH_ITER, ms * 1e6 / BENCH_ITER); \
} while(0)

int main(void)
{
    printf("\n=== mini-dusty-plasma Benchmarks ===\n\n");
    printf("Platform: C99/gcc, %d iterations per test\n\n", BENCH_ITER);

    double n_e = 1.0e15, T_e = 3.0 * DUSTY_EV_IN_K;
    double n_i = 1.02e15, T_i = 300.0, m_i = 6.64e-26;
    double a = 5.0e-6, Z_d = 2000.0, m_d = 1.15e-12;

    BENCH("Debye length (electron)", dust_debye_electron(n_e, T_e));
    BENCH("Debye length (ion)", dust_debye_ion(n_i, T_i));
    BENCH("Dust plasma frequency", dust_plasma_freq(1e10, Z_d, m_d));
    BENCH("OML electron current", dust_oml_electron_current(a, n_e, T_e, -3.0));
    BENCH("OML ion current", dust_oml_ion_current(a, n_i, T_i, m_i, -3.0));
    BENCH("Floating potential solve", dust_floating_potential_solve(
        a, n_e, n_i, T_e, T_i, m_i, CHARGE_MODEL_OML,
        0,0,0,0,300,5, -2.5*DUSTY_KB*T_e/DUSTY_EC, 1e-6, 20));
    BENCH("DAW dispersion", dust_acoustic_wave_dispersion(1000, 0.03, 4e-4).omega_r);
    BENCH("Yukawa potential", yukawa_potential(Z_d*DUSTY_EC, Z_d*DUSTY_EC, 2e-4, 4e-4));
    BENCH("Yukawa force magnitude", yukawa_force_magnitude(Z_d*DUSTY_EC, Z_d*DUSTY_EC, 2e-4, 4e-4));
    BENCH("Neutral collision freq", dust_neutral_collision_freq(a, 2.4e21, 300, 6.63e-26, m_d, 0.8));
    BENCH("Diffusion coefficient", dust_diffusion_coefficient(300, m_d, 100));
    BENCH("Einstein frequency", dust_einstein_frequency(Z_d*DUSTY_EC, m_d, 2e-4, 1.0));

    printf("\n=== Done ===\n\n");
    return 0;
}