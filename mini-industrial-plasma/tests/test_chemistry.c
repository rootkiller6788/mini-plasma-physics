#include "plasma_types.h"
#include <stdio.h>
#include <assert.h>

extern void setup_argon_chemistry(PlasmaSpecies*,int*,PlasmaReaction*,int*);
extern void setup_sf6_chemistry(PlasmaSpecies*,int*,PlasmaReaction*,int*);
extern void setup_silane_chemistry(PlasmaSpecies*,int*,PlasmaReaction*,int*);
extern void setup_oxygen_chemistry(PlasmaSpecies*,int*,PlasmaReaction*,int*);
extern int rate_matrix_init(RateMatrix*,int,int);
extern void rate_matrix_free(RateMatrix*);
extern double compute_rate_constant_arrhenius(double,double,double,double);
extern void compute_all_rate_constants(RateMatrix*,double,PlasmaReaction*);
extern double euler_chemistry_step(RateMatrix*,PlasmaReaction*,double);

int main(void) {
    double k = compute_rate_constant_arrhenius(3.0, 2.3e-14, 0.59, 15.76); assert(k>0);
    PlasmaSpecies s[12]; PlasmaReaction r[12]; int ns,nr;
    setup_argon_chemistry(s,&ns,r,&nr); assert(ns==5&&nr==7);
    RateMatrix rm; rate_matrix_init(&rm,ns,nr);
    rm.concentrations[0]=2.4e20; rm.concentrations[4]=1e16;
    compute_all_rate_constants(&rm,3.0,r); assert(rm.rate_constants[0]>0);
    double dt=euler_chemistry_step(&rm,r,1e-6); assert(dt>0);
    rate_matrix_free(&rm);
    setup_sf6_chemistry(s,&ns,r,&nr); assert(ns==10&&nr==5);
    setup_silane_chemistry(s,&ns,r,&nr); assert(ns==7&&nr==5);
    setup_oxygen_chemistry(s,&ns,r,&nr); assert(ns==7&&nr==6);
    printf("All chemistry tests passed.\n");
    return 0;
}
