#ifndef __tsp_
#define __tsp_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//configuration
const int L = 32; //number of cities
double position_x[L], position_xy[L][2];
int cities_x[L], cities_xy[L];;
int population_x[L], population_xy[L];
int new_population_x[L], new_population_xy[L];
int old_population_x[L], old_population_xy[L];
int new_n = 0;
int offspring_father_x[L], offspring_mother_x[L];
int offspring_father_xy[L], offspring_mother_xy[L];

//simulation
double F_x, F_xy;
int ngen = 1000;
int n_father_x, n_mother_x, n_father_xy, n_mother_xy;
int father[L], mother[L];
double walker;
double temp = 0.8;
double beta = 1/temp;
double ene_new_x, ene_old_x = 0, ene_new_xy, ene_old_xy = 0;
double weight_x, weight_xy;
int MCstep = 8000;
double acc_x, att_x;
double acc_xy, att_xy;

/*
// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x,stima_g;
double err_u, err_c, err_m, err_x, err_g;*/


//functions
void Input(void);
void InputGibbs(void);
void Restart(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void Measure(void);
void Creation(void);
void Cost(void);
void Check(void);
void Metropolis(void);
void Mutation(void);
void Crossover(void);
void Shuffle(int*,int);
int Pbc(int);
double Error(double,double,int);

#endif

