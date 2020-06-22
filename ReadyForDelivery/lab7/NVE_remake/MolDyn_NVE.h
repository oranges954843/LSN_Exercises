#include <cstring>

using namespace std;

//parameters, observables
const int m_props=400;
int n_props;
int iv,ik,it,ie;
int igofr;
double stima[m_props];

const double pi=3.1415927;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint, seed, restart;
int nsim = 0, nbins;
double delta;
string phase[3] = {"solid","liquid","gas"};
int index_ph = 0;
double eqstep, fs;
double bin_size;

// blocking average
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double err_e, err_v, err_k, err_t, err_gdir;
double stima_e, stima_v, stima_k, stima_t;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void WriteIstantValues(void);
void Restart(void);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);
void Reset(int);
void Accumulate(void);
void Averages(int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
