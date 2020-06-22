#ifndef __NVT__
#define __NVT__


#include "random.h"
int seed[4];
Random rnd;

//functions
double PsiQGS(double);
double PsiQ2P(double,double);
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
double Error(double,double,int);

// simulation
int nstep, nblk;
const int n_props = 2;
int iGS, i2P;
double delta;
double walker[n_props];

// averages
double blk_av[n_props], blk_norm;
double acc_GS, att_GS, acc_2P, att_2P;
double glob_av[n_props], glob_av2[n_props];
double stima_GS, stima_2P, err_GS, err_2P;
double xGS, yGS, zGS;
double x2P, y2P, z2P;


#endif

