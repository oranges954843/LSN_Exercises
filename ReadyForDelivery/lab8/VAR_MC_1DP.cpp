#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

#include "VAR_MC_1DP.h"

using namespace std;

#define mass 1.0
#define hbar 1.0
#define M 500

int main(){ 

	Input();
	ConfFinal();
	
	do{
		Restart();
		cout << "ene_old = " << ene_old << endl;
		for(int istep=1; istep<=M; istep++){
			MoveParameters();
		}
		ConfFinal();
		cout << "Acceptance rate for mu and sigma = " << acc_par/att_par << endl << endl;
		temp = temp - 0.01;

	}while(temp > 0.0);
	
	minimize = true;
	
	cout << "Simulation with best parameters....." << endl;
	cout << "mu = " << mu << ", sigma = " << sigma << endl;
	cout << "ene best = " << ene_old << endl;
	for(int iblk=1; iblk<=nblk; ++iblk){
		Reset(iblk);
		for(int istep=1; istep<=nstep; istep++){
			Move();
			FillConfig();
			Measure();
			Accumulate();
		}
		Averages(iblk);
		}
	ConfFinal();
	cout << "Simulation ended" << endl;
				
	return 0;
}


void Input(void){
	ifstream ReadInput, ReadConf;
	int p1, p2;

	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
		
	cout << "----------------------------" << endl << endl;
	cout << "Single quantum particle in 1D" << endl;
	cout << "Monte Carlo simulation" << endl;
	cout << "External potential v(r) = x^4 - 2.5 * x^2" << endl;
	cout << "Wave function proportional to ";
	cout << "exp(-(x-mu)^2/(2*sigma^2)) + exp(-(x+mu)^2/(2*sigma^2))" << endl;
	cout << "The program uses hbar = 1 and mass = 1 units " << endl << endl;

	ReadInput.open("initial_values.dat");

	ReadInput >> mu;

	cout << "Initial mu = " << mu << endl;

	ReadInput >> sigma;
	cout << "Initial sigma = " << sigma << endl;

	ReadInput >> temp;
	cout << "Initial temperature = " << temp << endl;

	ReadInput >> x;
	cout << "Input position = " << x << endl;

	ReadInput >> delta; 
	ReadInput >> nblk; 
	ReadInput >> nstep;
	ReadInput >> eqstep; 

	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in a block = " << nstep << endl;
	cout << "Number of steps to equilibrate the system = " << eqstep << endl;
	ReadInput.close();

	Measure();
	cout << "Input energy = " << walker[ie] << endl << endl;

	ie = 0;

	n_props = 1;

	attempted = 0;
	accepted = 0;
	cout << "Performing equilibration......" << endl;
	for(int istep=1; istep <= eqstep; ++istep) Move();
	cout << "Equilibration ended" << endl;
	cout << "initial acceptance rate : ";
	cout << accepted/attempted << endl;

	for(int iblk=1; iblk<=nblk; iblk++){
		Reset(iblk);
		for(int istep=1; istep <= nstep; ++istep){
			Move();
			Measure();
			Accumulate();
		}
		Averages(iblk);
	}
	ene_old = glob_av[ie] / (double)nblk;
	cout << "Initial energy = " << ene_old << endl;
	cout << "Initial position = " << x << endl << endl;
}


void Restart(){
	ifstream ReadConf("config.final");

	ReadConf >> x;
	ReadConf >> mu_old;
	ReadConf >> sigma_old;
	ReadConf.close();

	beta = 1 / temp;
	cout << "Current temperature = " << temp << endl;
	acc_par = 0;
	att_par = 0;
}


void MoveParameters(void){

	ofstream WriteEne("istant_etot.out",ios::app);
	
	WriteEne << ene_old << endl;
	WriteEne.close();

	mu = rnd.Rannyu(-1,1);
	sigma = rnd.Rannyu(0.6,1.4);
	

	for(int iblk=1; iblk<=nblk; iblk++){
		Reset(iblk);
		for(int istep=1; istep<=nstep; istep++){
			Move();
			Measure();
			Accumulate();
		}
		Averages(iblk);
	}
	ene_new = glob_av[ie] / (double)nblk;
	
	p = Boltzmann(ene_old,ene_new);	
	if(p >= rnd.Rannyu()){
		
		acc_par++;
		ene_old = ene_new;
	}
	else{
		mu = mu_old;
		sigma = sigma_old;
	} 
	att_par++;
}

	
double Boltzmann(double Eold, double Enew){
	
	double r;
	return r = exp(-beta * (Enew - Eold));
}


void Move(void){

	double p, xold, xnew;

	xold = x;
	
	xnew = x + delta*(rnd.Rannyu() - 0.5);

	p = Weight(xold,xnew);
	if(p >= rnd.Rannyu()){

		x = xnew;
		accepted = accepted + 1.0;
	}
	attempted = attempted + 1.0;
}


double Weight(double xold, double xnew){

	double delta, rnew, rold;

	rnew = xnew * mu / pow(sigma,2);
	rold = xold * mu / pow(sigma,2);
	delta = (xold * xold - xnew * xnew) / pow(sigma,2);

	return exp(delta)*pow(cosh(rnew),2)/pow(cosh(rold),2);
}


void Measure(){

	double ekin, epot;
	double argo;
	double lambda = hbar * hbar / (2 * mass);

	argo = x * mu / pow(sigma,2);

	ekin = (x*x + mu*mu - sigma*sigma) * cosh(argo) - 2*x*mu * sinh(argo);
	ekin = - lambda * ekin / (pow(sigma,4) * cosh(argo));

	epot = pow(x,4) - 2.5 * pow(x,2);
	walker[ie] = ekin + epot;
}


void Reset(int iblk){

	if(iblk == 1){
		for(int i=0; i<n_props; ++i){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i){
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void){

	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){
    
	stima_etot = blk_av[ie]/blk_norm; 
	glob_av[ie] += stima_etot;

	if(minimize == true){

		glob_av2[ie] += stima_etot*stima_etot;
		err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

		ofstream WriteEtot;
		WriteEtot.open("ave_etot.final",ios::app);
		WriteEtot << iblk << "   " << glob_av[ie]/(double)iblk;
		WriteEtot << "   " << err_etot << endl;
		WriteEtot.close();
	}
}


void FillConfig(void){
	
	ofstream WriteSampledPos("sampled_positions.out",ios::app);
	WriteSampledPos << x << endl;
	WriteSampledPos.close();
}

 
void ConfFinal(void){

	ofstream WriteConf;

	WriteConf.open("config.final");
	WriteConf << x << endl;
	WriteConf << mu << endl;
	WriteConf << sigma << endl;
	WriteConf.close();

	rnd.SaveSeed();
}


double Error(double sum, double sum2, int iblk){

	if( iblk == 1 ) return 0.0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
