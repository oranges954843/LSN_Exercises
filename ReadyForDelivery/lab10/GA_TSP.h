#ifndef __tsp_
#define __tsp_

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//configuration
const int L = 32; //number of cities
const int N = 1000; //number of elements in the population
const int Nmigr = 5;
double position_x[L], position_xy[L][2];
int cities_x[L], cities_xy[L];;
int population_x[N][L], population_xy[N][L];
int new_population_x[N][L], new_population_xy[N][L];
int new_n = 0;
int casual_n = 0;
int offspring_father_x[L], offspring_mother_x[L];
int offspring_father_xy[L], offspring_mother_xy[L];

//simulation
double F_x[N] = {0}, F_xy[N] = {0};
double fitness_x[N] = {0}, fitness_xy[N] = {0};
int ngen = 250;
int n_father_x, n_mother_x, n_father_xy, n_mother_xy;
int father[L], mother[L];
double walker;
int highest_nxy, lowest_nxy;
int order[N];

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
void Selection(void);
void Mutation(void);
void Crossover(void);
void Shuffle(int*,int);
int Pbc(int);
double Error(double,double,int);

#endif

