/**
 * bench_particle_pusher.c — Boris Pusher Performance Benchmark
 *
 * Benchmarks the Boris particle pusher for large ensembles.
 */
#include "plasma_particle.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(void) {
    printf("=== Boris Pusher Benchmark ===\n");

    int N = 1000000;
    printf("Particles: %d\n", N);
    printf("Steps: 100\n\n");

    ParticleEnsemble ensemble;
    ensemble.n_particles = N;
    ensemble.weight = 1.0;
    ensemble.particles = (Particle3D*)malloc(N * sizeof(Particle3D));

    for (int i = 0; i < N; i++) {
        ensemble.particles[i].x = (double)(rand()) / RAND_MAX;
        ensemble.particles[i].y = (double)(rand()) / RAND_MAX;
        ensemble.particles[i].z = 0.0;
        ensemble.particles[i].vx = 1e5 * ((double)rand() / RAND_MAX - 0.5);
        ensemble.particles[i].vy = 1e5 * ((double)rand() / RAND_MAX - 0.5);
        ensemble.particles[i].vz = 1e5 * ((double)rand() / RAND_MAX - 0.5);
        ensemble.particles[i].mass = M_ELECTRON;
        ensemble.particles[i].charge = -E_CHARGE;
        ensemble.particles[i].alive = 1;
    }

    clock_t start = clock();
    double dt = 1e-12;
    for (int step = 0; step < 100; step++) {
        for (int i = 0; i < N; i++) {
            boris_push(&ensemble.particles[i],
                       1e3, 0.0, 0.0,  /* E_x, E_y, E_z */
                       0.0, 0.0, 5.0,  /* B_x, B_y, B_z */
                       dt);
        }
    }
    clock_t end = clock();

    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double pushes_per_sec = (N * 100.0) / elapsed;
    double ns_per_push = elapsed * 1e9 / (N * 100.0);

    printf("Total time: %.3f s\n", elapsed);
    printf("Throughput: %.2e pushes/s\n", pushes_per_sec);
    printf("Per push:   %.2f ns\n", ns_per_push);

    free(ensemble.particles);
    return 0;
}
