
#ifndef __NVT__
#define __NVT__

#include <string>

using namespace std;

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, ie;
double bin_size,nbins,sd;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_etot, err_etot;
double acc_par, att_par;
double p;

//configuration
const int m_part=108;
double x;

// thermodynamical state
double mu, sigma;
double beta, temp;
double ene_old, ene_new;
double mu_old, sigma_old;

// simulation
int nstep, nblk;
double delta, eqstep;
bool minimize = false;

//functions
void Input(void);
void Reset(int);
void Restart(void);
void Accumulate(void);
void Averages(int);
void Move(void);
void MoveParameters(void);
void ConfFinal(void);
void ConfStart(void);
void FillConfig(void);
double Boltzmann(double,double);
double Weight(double,double);
void Measure(void);
double Error(double,double,int);

#endif

