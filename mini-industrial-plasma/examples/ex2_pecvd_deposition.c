#include "sheath.h"
#include "surface.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Example 2: PECVD a-Si:H Deposition ===\n");
    double T_sub = 523.0;
    double n_SiH3 = 1e19;
    double Ea = 0.2 * 1.602e-19;
    double R0 = 1.5e-9;
    double Rdep = pecvd_deposition_rate(T_sub, n_SiH3, Ea, R0);
    printf("Substrate temperature: %.0f K (%.0f C)\n", T_sub, T_sub-273);
    printf("Deposition rate: %.1f nm/min\n", Rdep*1e9*60);
    double H_cont = pecvd_hydrogen_content(T_sub, Rdep);
    printf("Hydrogen content: %.1f at%%\n", H_cont);
    double stress = pecvd_film_stress(T_sub, 573.0, 2.5e-6, 3.5e-6, 80e9, -200e6);
    printf("Film stress: %.0f MPa\n", stress/1e6);
    double theta = surface_coverage_langmuir(1.0, 0.3, T_sub, 0.15);
    printf("Surface H coverage: %.2f\n", theta);
    double S = sticking_probability_arrhenius(T_sub, 0.05*1.602e-19, 0.5);
    printf("SiH3 sticking probability: %.3f\n", S);
    double u[5] = {520, 510, 515, 525, 518};
    double U = pecvd_thickness_uniformity(u, 5, 0.15);
    printf("Thickness uniformity (3-sigma): %.1f%%\n", U);
    printf("\nTypical a-Si:H PECVD: 250C, SiH4/H2=1:5, 50W RF\n");
    return 0;
}
