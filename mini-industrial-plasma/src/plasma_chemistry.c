/* plasma_chemistry.c - Plasma Chemistry and Chemical Kinetics
 * Reference: Lieberman 2005, Ch.7-8
 * Course: MIT 22.611, Berkeley EECS 245
 */

#include "plasma_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

extern void plasma_species_init(PlasmaSpecies *sp, const char *name,
    double mass_amu, int charge, double ioniz_eV, double excit_eV,
    double polariz_m3, PlasmaSpeciesType type);
extern void plasma_reaction_init(PlasmaReaction *rxn, const char *eq,
    ReactionType type, double thr, double A, double n,
    double smax, double Eatmax, int *react, int nr, int *prod, int np);

int rate_matrix_init(RateMatrix *rm, int n_species, int n_reactions)
{
    if (!rm || n_species <= 0 || n_reactions <= 0) return -1;
    rm->n_species = n_species;
    rm->n_reactions = n_reactions;
    rm->concentrations = calloc(n_species, sizeof(double));
    rm->production_rates = calloc(n_species, sizeof(double));
    rm->loss_rates = calloc(n_species, sizeof(double));
    rm->rate_constants = calloc(n_reactions, sizeof(double));
    if (!rm->concentrations || !rm->production_rates
        || !rm->loss_rates || !rm->rate_constants) {
        free(rm->concentrations); free(rm->production_rates);
        free(rm->loss_rates); free(rm->rate_constants);
        return -1;
    }
    return 0;
}

void rate_matrix_free(RateMatrix *rm)
{
    if (!rm) return;
    free(rm->concentrations); free(rm->production_rates);
    free(rm->loss_rates); free(rm->rate_constants);
    rm->concentrations = NULL; rm->production_rates = NULL;
    rm->loss_rates = NULL; rm->rate_constants = NULL;
    rm->n_species = 0; rm->n_reactions = 0;
}

double compute_rate_constant_arrhenius(double T_e_eV, double A,
                                        double n, double E_a)
{
    if (T_e_eV <= 0.0) return 0.0;
    if (E_a > 0.0 && T_e_eV < E_a * 0.01) return 0.0;
    return A * pow(T_e_eV, n) * exp(-E_a / T_e_eV);
}

void compute_all_rate_constants(RateMatrix *rm, double T_e_eV,
                                 PlasmaReaction *reactions)
{
    if (!rm || !reactions || T_e_eV <= 0.0) return;
    for (int j = 0; j < rm->n_reactions; j++) {
        PlasmaReaction *rxn = &reactions[j];
        rm->rate_constants[j] = compute_rate_constant_arrhenius(
            T_e_eV, rxn->rate_constant_A,
            rxn->rate_constant_n, rxn->threshold);
    }
}

void compute_rate_matrix(RateMatrix *rm, PlasmaReaction *reactions)
{
    if (!rm || !reactions) return;
    for (int i = 0; i < rm->n_species; i++) {
        rm->production_rates[i] = 0.0;
        rm->loss_rates[i] = 0.0;
    }
    for (int j = 0; j < rm->n_reactions; j++) {
        PlasmaReaction *rxn = &reactions[j];
        double k = rm->rate_constants[j];
        if (k <= 0.0) continue;
        double rp = k;
        for (int ri = 0; ri < rxn->n_reactants; ri++) {
            int id = rxn->reactant_ids[ri];
            if (id >= 0 && id < rm->n_species)
                rp *= rm->concentrations[id];
        }
        for (int ri = 0; ri < rxn->n_reactants; ri++) {
            int id = rxn->reactant_ids[ri];
            if (id >= 0 && id < rm->n_species)
                rm->loss_rates[id] += rp;
        }
        for (int pi = 0; pi < rxn->n_products; pi++) {
            int id = rxn->product_ids[pi];
            if (id >= 0 && id < rm->n_species)
                rm->production_rates[id] += rp;
        }
    }
}

double euler_chemistry_step(RateMatrix *rm, PlasmaReaction *reactions,
                             double dt_max)
{
    if (!rm || !reactions || dt_max <= 0.0) return 0.0;
    compute_rate_matrix(rm, reactions);
    double dt = dt_max;
    for (int i = 0; i < rm->n_species; i++) {
        double net_rate = rm->production_rates[i] - rm->loss_rates[i];
        double abs_rate = fabs(net_rate);
        double n_i = rm->concentrations[i];
        if (abs_rate > 0.0 && n_i > 1.0) {
            double mdt = 0.1 * n_i / abs_rate;
            if (mdt < dt) dt = mdt;
        }
    }
    if (dt < 1e-15) dt = 1e-15;
    for (int i = 0; i < rm->n_species; i++) {
        double net_rate = rm->production_rates[i] - rm->loss_rates[i];
        rm->concentrations[i] += dt * net_rate;
        if (rm->concentrations[i] < 0.0) rm->concentrations[i] = 0.0;
    }
    return dt;
}

void setup_argon_chemistry(PlasmaSpecies *s, int *ns, PlasmaReaction *r, int *nr)
{
    if(!s||!ns||!r||!nr)return;
    *ns=5;
    plasma_species_init(&s[0],"Ar",39.948,0,0.0,0.0,1.64e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[1],"Ar+",39.948,1,27.63,0.0,1.0e-30,SPECIES_ION_POS);
    plasma_species_init(&s[2],"Ar*",39.948,0,4.21,0.0,2.0e-30,SPECIES_METASTABLE);
    plasma_species_init(&s[3],"Ar**",39.948,0,2.64,11.55,3.0e-30,SPECIES_EXCITED);
    plasma_species_init(&s[4],"e",5.486e-4,-1,0.0,0.0,0.0,SPECIES_ELECTRON);
    *nr=7;
    {int ra[]={4,0};int pr[]={4,4,1};
     plasma_reaction_init(&r[0],"e+Ar->2e+Ar+",REACTION_IONIZATION,15.76,2.3e-14,0.59,3e-20,50.0,ra,2,pr,3);}
    {int ra[]={4,0};int pr[]={4,2};
     plasma_reaction_init(&r[1],"e+Ar->e+Ar*",REACTION_EXCITATION,11.55,5.0e-15,0.5,1e-20,20.0,ra,2,pr,2);}
    {int ra[]={4,2};int pr[]={4,4,1};
     plasma_reaction_init(&r[2],"e+Ar*->2e+Ar+",REACTION_IONIZATION,4.21,6.8e-14,0.5,5e-20,15.0,ra,2,pr,3);}
    {int ra[]={4,2};int pr[]={4,3};
     plasma_reaction_init(&r[3],"e+Ar*->e+Ar**",REACTION_EXCITATION,2.1,3.0e-14,0.3,3e-20,10.0,ra,2,pr,2);}
    {int ra[]={3};int pr[]={2};
     plasma_reaction_init(&r[4],"Ar**->Ar*+hv",REACTION_RECOMBINATION,0.0,3.2e7,0.0,0.0,0.0,ra,1,pr,1);}
    {int ra[]={2,2};int pr[]={1,0,4};
     plasma_reaction_init(&r[5],"Ar*+Ar*->Ar++Ar+e",REACTION_PENNING,0.0,6.2e-16,0.0,0.0,0.0,ra,2,pr,3);}
    {int ra[]={1,4,4};int pr[]={0,4};
     plasma_reaction_init(&r[6],"Ar++2e->Ar+e",REACTION_RECOMBINATION,0.0,8.75e-39,-4.5,0.0,0.0,ra,3,pr,2);}
}

void setup_sf6_chemistry(PlasmaSpecies *s, int *ns, PlasmaReaction *r, int *nr)
{
    if(!s||!ns||!r||!nr)return;
    *ns=10;
    plasma_species_init(&s[0],"SF6",146.06,0,0.0,0.0,6.5e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[1],"SF4",108.06,0,0.0,0.0,5.0e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[2],"SF2",70.06,0,0.0,0.0,4.0e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[3],"F",19.00,0,0.0,0.0,0.5e-30,SPECIES_RADICAL);
    plasma_species_init(&s[4],"F2",38.00,0,0.0,0.0,1.4e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[5],"SiF4",104.08,0,0.0,0.0,3.0e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[6],"e",5.486e-4,-1,0.0,0.0,0.0,SPECIES_ELECTRON);
    plasma_species_init(&s[7],"SF5+",127.06,1,15.0,0.0,3.0e-30,SPECIES_ION_POS);
    plasma_species_init(&s[8],"SF3+",89.06,1,12.0,0.0,2.5e-30,SPECIES_ION_POS);
    plasma_species_init(&s[9],"F-",19.00,-1,0.0,0.0,1.0e-30,SPECIES_ION_NEG);
    *nr=5;
    {int ra[]={6,0};int pr[]={6,1,3,3};plasma_reaction_init(&r[0],"e+SF6->e+SF4+2F",REACTION_DISSOCIATION,4.0,1e-14,0.3,1e-20,10.0,ra,2,pr,4);}
    {int ra[]={6,1};int pr[]={6,2,3,3};plasma_reaction_init(&r[1],"e+SF4->e+SF2+2F",REACTION_DISSOCIATION,3.0,2e-14,0.3,2e-20,8.0,ra,2,pr,4);}
    {int ra[]={6,0};int pr[]={6,6,7,3};plasma_reaction_init(&r[2],"e+SF6->SF5++F+2e",REACTION_IONIZATION,15.29,4e-15,0.5,5e-21,30.0,ra,2,pr,4);}
    {int ra[]={6,0};int pr[]={9};plasma_reaction_init(&r[3],"e+SF6->F-+SF5",REACTION_ATTACHMENT,0.0,2.2e-13,-1.0,1e-19,0.3,ra,2,pr,1);}
    {int ra[]={3,3,3,3};int pr[]={5};plasma_reaction_init(&r[4],"4F->SiF4",REACTION_RADICAL_RECOMB,0.0,1e-16,0.0,0.0,0.0,ra,4,pr,1);}
}

void setup_silane_chemistry(PlasmaSpecies *s, int *ns, PlasmaReaction *r, int *nr)
{
    if(!s||!ns||!r||!nr)return;
    *ns=7;
    plasma_species_init(&s[0],"SiH4",32.12,0,0.0,0.0,4.6e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[1],"SiH3",31.12,0,0.0,0.0,5.0e-30,SPECIES_RADICAL);
    plasma_species_init(&s[2],"SiH2",30.11,0,0.0,0.0,5.5e-30,SPECIES_RADICAL);
    plasma_species_init(&s[3],"H",1.008,0,0.0,0.0,0.67e-30,SPECIES_RADICAL);
    plasma_species_init(&s[4],"H2",2.016,0,0.0,0.0,0.8e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[5],"Si2H6",62.22,0,0.0,0.0,8.0e-30,SPECIES_DIMER);
    plasma_species_init(&s[6],"e",5.486e-4,-1,0.0,0.0,0.0,SPECIES_ELECTRON);
    *nr=5;
    {int ra[]={6,0};int pr[]={6,1,3};plasma_reaction_init(&r[0],"e+SiH4->e+SiH3+H",REACTION_DISSOCIATION,8.75,4.0e-14,0.3,8e-21,15.0,ra,2,pr,3);}
    {int ra[]={6,0};int pr[]={6,2,3,3};plasma_reaction_init(&r[1],"e+SiH4->e+SiH2+2H",REACTION_DISSOCIATION,10.5,1.5e-14,0.4,3e-21,18.0,ra,2,pr,4);}
    {int ra[]={6,4};int pr[]={6,3,3};plasma_reaction_init(&r[2],"e+H2->e+2H",REACTION_DISSOCIATION,8.9,1.0e-14,0.3,1e-20,15.0,ra,2,pr,3);}
    {int ra[]={3,0};int pr[]={1,4};plasma_reaction_init(&r[3],"H+SiH4->SiH3+H2",REACTION_CHARGE_EXCHANGE,0.0,2.8e-18,0.0,0.0,0.0,ra,2,pr,2);}
    {int ra[]={1,1};int pr[]={5};plasma_reaction_init(&r[4],"SiH3+SiH3->Si2H6",REACTION_RADICAL_RECOMB,0.0,1.5e-16,0.0,0.0,0.0,ra,2,pr,1);}
}

void setup_oxygen_chemistry(PlasmaSpecies *s, int *ns, PlasmaReaction *r, int *nr)
{
    if(!s||!ns||!r||!nr)return;
    *ns=7;
    plasma_species_init(&s[0],"O2",32.00,0,0.0,0.0,1.6e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[1],"O",16.00,0,0.0,0.0,0.8e-30,SPECIES_RADICAL);
    plasma_species_init(&s[2],"O3",48.00,0,0.0,0.0,2.8e-30,SPECIES_NEUTRAL);
    plasma_species_init(&s[3],"O2+",32.00,1,12.06,0.0,1.0e-30,SPECIES_ION_POS);
    plasma_species_init(&s[4],"O+",16.00,1,13.62,0.0,0.5e-30,SPECIES_ION_POS);
    plasma_species_init(&s[5],"O-",16.00,-1,0.0,0.0,1.0e-30,SPECIES_ION_NEG);
    plasma_species_init(&s[6],"e",5.486e-4,-1,0.0,0.0,0.0,SPECIES_ELECTRON);
    *nr=6;
    {int ra[]={6,0};int pr[]={1,5};plasma_reaction_init(&r[0],"e+O2->O+O-",REACTION_ATTACHMENT,4.2,8.8e-17,-1.0,1e-22,6.7,ra,2,pr,2);}
    {int ra[]={6,0};int pr[]={6,6,3};plasma_reaction_init(&r[1],"e+O2->2e+O2+",REACTION_IONIZATION,12.06,9.0e-16,0.5,2.5e-20,31.0,ra,2,pr,3);}
    {int ra[]={6,0};int pr[]={6,1,1};plasma_reaction_init(&r[2],"e+O2->e+2O",REACTION_DISSOCIATION,6.0,4.2e-15,0.3,1e-20,10.0,ra,2,pr,3);}
    {int ra[]={5,3};int pr[]={1,0};plasma_reaction_init(&r[3],"O-+O2+->O+O2",REACTION_RECOMBINATION,0.0,2.0e-13,-0.5,0.0,0.0,ra,2,pr,2);}
    {int ra[]={5,4};int pr[]={1,1};plasma_reaction_init(&r[4],"O-+O+->2O",REACTION_RECOMBINATION,0.0,2.0e-13,-0.5,0.0,0.0,ra,2,pr,2);}
    {int ra[]={1,0,0};int pr[]={2,0};plasma_reaction_init(&r[5],"O+O2+O2->O3+O2",REACTION_RADICAL_RECOMB,0.0,6.0e-46,0.0,0.0,0.0,ra,3,pr,2);}
}