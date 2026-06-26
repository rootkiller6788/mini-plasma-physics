#include "plasma_types.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

extern void plasma_state_init(PlasmaState*,double,double,double,double,double,double);
extern void plasma_species_init(PlasmaSpecies*,const char*,double,int,double,double,double,PlasmaSpeciesType);
extern int plasma_quasineutrality_check(double,double,double,double);
extern double plasma_coulomb_logarithm(double,double);
extern double plasma_beta(double,double,double,double);
extern double hall_parameter(double,double);
extern double ionization_degree(double,double);
extern void industrial_plasma_defaults(const char*,PlasmaState*,DischargeParams*,DischargeGeometry*);
extern void estimate_plasma_density_rf(double,double,double,double,double,double*,double*);

int main(void) {
    PlasmaState ps; plasma_state_init(&ps, 1e16, 3.0, 300, 300, 5.0, 0.0);
    assert(ps.n_e == 1e16 && ps.lam_D > 0);
    PlasmaSpecies sp; plasma_species_init(&sp, "Ar", 39.948, 0, 0.0, 0.0, 1.64e-30, SPECIES_NEUTRAL);
    assert(strcmp(sp.name, "Ar") == 0);
    assert(plasma_quasineutrality_check(1e16, 3.0, 0.1, 1e6) == 1);
    double lnL = plasma_coulomb_logarithm(1e16, 3.0); assert(lnL > 2);
    double b = plasma_beta(1e16, 3.0, 300, 0.01); assert(b >= 0);
    double h = hall_parameter(0.01, 1e7); assert(h >= 0);
    double al = ionization_degree(1e16, 2.4e20); assert(al > 0 && al < 1);
    double ne,Te; estimate_plasma_density_rf(500,5.0,0.15,0.1,13.56e6,&ne,&Te); assert(ne>0&&Te>0);
    PlasmaState ps2; DischargeParams dp; DischargeGeometry g;
    industrial_plasma_defaults("RIE", &ps2, &dp, &g); assert(dp.power>0);
    printf("All plasma tests passed.\n");
    return 0;
}
