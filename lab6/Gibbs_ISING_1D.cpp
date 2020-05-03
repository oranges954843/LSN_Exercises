#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main() { 
	Input(); //Inizialization
accepted = 0;
attempted = 0;
	for(int iblk=1; iblk <= nblk; ++iblk) { //Simulation
		Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep) {
			Move();
			Measure();
			Accumulate(); //Update block averages
			//cout << "attempted = " << attempted << endl;
//cout << "accepted = " << accepted << endl;
		}
		Averages(iblk);   //Print results for current block
	}
	ConfFinal(); //Write final configuration
accepted = 0;
attempted = 0;
		temp = temp - 0.1;
	do {
		Restart();
		for (int iblk=1; iblk <= nblk; ++iblk) {
			Reset(iblk);
			for (int istep=1; istep <= nstep; ++istep) {
				Move();
				Measure();
				Accumulate();
		//cout << "attempted = " << attempted << endl;
		//cout << "accepted = " << accepted << endl;
			}
			Averages(iblk);
		}
		ConfFinal();
		temp = temp - 0.1;
	} while(temp > 0.4);

	return 0;
}


void Input(void) {
	ifstream ReadInput;

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

	//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;   
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
	//Read input informations
	ReadInput.open("input.dat");

	ReadInput >> temp;

	beta = 1.0/temp;
	cout << "Temperature = " << temp << endl;

	ReadInput >> nspin;
	cout << "Number of spins = " << nspin << endl;

	ReadInput >> J;
	cout << "Exchange interaction = " << J << endl;

	ReadInput >> h;
	cout << "External field = " << h << endl << endl;

	ReadInput >> nblk;

	ReadInput >> nstep;

	cout << "The program perform Gibbs moves" << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();


	//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
 
	n_props = 4; //Number of observables

	//initial configuration
	for (int i=0; i<nspin; ++i) {
		if(rnd.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
  
	//Evaluate energy etc. of the initial configuration
	Measure();

	//Print initial values for the potential energy and virial
	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
	cout << "Initial heat capacity = " << walker[ic]/(double)nspin << endl;
	cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
	cout << "Initial susceptibility = " << walker[ix]/(double)nspin << endl;

	//Equilibration
	for (int i=0; i<10000; i++) Move();  
}


void Move() {
	int o;
	double p, prob_up, prob_down;
	double energy_up, energy_down;
	//pick a spin randomly in the array
	o = (int)(rnd.Rannyu()*nspin);
	energy_up = Boltzmann(-2,o);
	energy_down = Boltzmann(2,o);
	prob_up = pow(1 + exp(-beta * energy_up),-1);
	prob_down = pow(1 + exp(-beta * energy_down),-1);
	p = rnd.Rannyu();
	if(prob_up == prob_down ){
		if(p < 0.5) s[o] = 1;
		else s[o] = -1;
	}else if(prob_up != prob_down){
			if(p < prob_up) s[0] = 1;
			else s[o] = -1;
	}
}


double Boltzmann(double sm, int ip) {
	return -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) ;//- h * sm;
}

void Measure() {
	/*int bin*/;
	double u = 0.0, u_free = 0.0, m = 0.0, chi = 0.0;
	double sum_spin = 0.0;
//cout << endl << "-----------------" << endl;
	//cycle over spins
	for (int i=0; i<nspin; ++i) {
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		u_free += -J * s[i] * s[Pbc(i+1)];
		//cout << s[i] << "   ";
		sum_spin += s[i];
		m += s[i] * exp(-beta*u);
		chi += s[i]*s[i];
	}
	m = sum_spin;
	chi = beta * pow(sum_spin,2); 
	walker[iu] = u_free;
	walker[ic] = u_free * u_free;
	walker[im] = m;
	walker[ix] = chi;
	//cout << endl << "energy = " << u_free << endl;
}


void Reset(int iblk) { //Reset block averages

	if(iblk == 1) {
		for(int i=0; i<n_props; ++i) {
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void) { //Update block averages

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) { //Print results for current block
	ofstream Ene, Heat, Mag, Chi;
	double c;

    
	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy per spin
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u=Error(glob_av[iu],glob_av2[iu],iblk);

	c = blk_av[ic]/blk_norm - pow(stima_u * nspin,2);
	stima_c = beta * beta * c / (double)nspin; //Heat capacity
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c=Error(glob_av[ic],glob_av2[ic],iblk);

	stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m=Error(glob_av[im],glob_av2[im],iblk);

	stima_x = blk_av[ix]/blk_norm/(double)nspin; //Susceptibility per spin 
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x=Error(glob_av[ix],glob_av2[ix],iblk);

	//Writing measurments of one block
	Ene.open("ene_gibbs.out",ios::app);
	Heat.open("heat_gibbs.out",ios::app);
	Mag.open("mag_gibbs.out",ios::app);
	Chi.open("chi_gibbs.out",ios::app);

	Ene << iblk <<  "   " << stima_u << "   " << glob_av[iu]/(double)iblk <<  "   " << err_u << endl;
	Heat << iblk <<  "   " << stima_c <<  "   " << glob_av[ic]/(double)iblk <<  "   " << err_c << endl;
	Mag << iblk <<  "   " << stima_m <<  "   " << glob_av[im]/(double)iblk <<  "   " << err_m << endl;
	Chi << iblk <<  "   " << stima_x <<  "   " << glob_av[ix]/(double)iblk <<  "   " << err_x << endl;

	Ene.close();
	Heat.close();
	Mag.close();
	Chi.close();
}


void ConfFinal(void) {
	ofstream WriteConf, WriteEne, WriteHeat, WriteMag, WriteChi;

	WriteConf.open("config_gibbs.final");
	for (int i=0; i<nspin; ++i) {
		WriteConf << s[i] << endl;
	}
	WriteConf.close();

	//Writing final results
	WriteEne.open("eneVStemp_gibbs.out",ios::app);
	WriteHeat.open("heatVStemp_gibbs.out",ios::app);
	WriteMag.open("magVStemp_gibbs.out",ios::app);
	WriteChi.open("chiVStemp_gibbs.out",ios::app);

	WriteEne << temp << "   ";
	WriteEne << glob_av[iu]/(double)nblk << "   ";
	WriteEne << err_u << endl;

	WriteHeat << temp << "   ";
	WriteHeat << glob_av[ic]/(double)nblk << "   ";
	WriteHeat << err_c << endl;

	WriteMag << temp << "   ";
	WriteMag << glob_av[im]/(double)nblk << "   ";
	WriteMag << err_m << endl;

	WriteChi << temp << "   ";
	WriteChi << glob_av[ix]/(double)nblk << "   ";
	WriteChi << err_c << endl;

	WriteEne.close();
	WriteHeat.close();
	WriteMag.close();
	WriteChi.close();

	rnd.SaveSeed();
}


int Pbc(int i) { //Algorithm for periodic boundary conditions
	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;
	return i;
}


double Error(double sum, double sum2, int iblk) {
	if(iblk == 0) return 0;
	if(sum2/(double)iblk >= pow(sum/(double)iblk,2)){
		return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
	}else return 0;
}


void Restart(void) {
	//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.out");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
	//New input informations
	beta = 1.0/temp;
	//New initial values
	Measure();
}

// aggiungere anche controllo sulla radie q
///stimare in averages la capacita termica con l energia/////////
/*
//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
	energy_old += -J * s[i] * s[Pbc(i+1)]
	o = (int)(rnd.Rannyu()*nspin);
	if(metro==1) {
		Boltzmann(o,i);
	} else {//Gibbs sampling
		gibbs();
	}
}
*/
