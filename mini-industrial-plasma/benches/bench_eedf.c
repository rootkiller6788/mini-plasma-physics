#include "eedf.h"
#include "sheath.h"
#include <stdio.h>
#include <time.h>

int main(void) {
    clock_t start, end;
    int N = 100000;
    start = clock();
    for (int i = 0; i < N; i++) maxwellian_eedf(i*0.01, 3.0);
    end = clock();
    double t_m = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Maxwellian EEDF: %d evals in %.4f s\n", N, t_m);
    start = clock();
    for (int i = 0; i < N; i++) child_langmuir_current_density(200.0, 0.001, 6.63e-26);
    end = clock();
    double t_cl = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Child-Langmuir: %d evals in %.4f s\n", N, t_cl);
    CrossSectionModel cs = {2.5e-20, 25.0, 1.8, 1e-20, 10.0};
    start = clock();
    for (int i = 0; i < N/100; i++) rate_coefficient_maxwellian(&cs, 3.0);
    end = clock();
    double t_k = (double)(end-start)/CLOCKS_PER_SEC;
    printf("Rate coeff (1000pt): %d evals in %.4f s\n", N/100, t_k);
    printf("\nBenchmark complete.\n");
    return 0;
}
