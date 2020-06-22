#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

#define eqstep 10000

int main() { 

	Input(); //Inizialization

	cout << "Equilibration with ";
	cout << eqstep << " steps....." << endl;
	for(int i=0; i<eqstep; i++){
		Move();
	}
	cout << "Equilibration ended" << endl << endl;
	cout << "Initial acceptance rate: " << accepted/attempted << endl << endl;

	for(int iblk=1; iblk <= nblk; iblk++){ //Simulation
		Reset(iblk);
		for(int istep=1; istep <= nstep; istep++){
			Move();
			Measure();
			Accumulate();
		}
		Averages(iblk);
	}
	ConfFinal(); //Write final configuration

	temp = temp - 0.1;
	do{
		Restart();
		for(int iblk=1; iblk <= nblk; ++iblk){
			Reset(iblk);
			for(int istep=1; istep <= nstep; ++istep){
				Move();
				Measure();
				Accumulate();
			}
			Averages(iblk);
		}
		ConfFinal();
		temp = temp - 0.1;
	} while(temp > 0.4);
	
	cout << "Turning on the interaction" << endl << endl;
	turn_on = true;
	for (int i=0; i<nspin; ++i) {
		if(rnd.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
	cout << "Re-equilibration....." << endl;
	for (int i=0; i<eqstep; i++) Move();
	cout << "Equilibration ended" << endl << endl;
	temp = 2.0;

	for(int iblk=1; iblk <= nblk; iblk++){ //Simulation
		Reset(iblk);
		for(int istep=1; istep <= nstep; istep++){
			Move();
			Measure();
			Accumulate();
		}
		Averages(iblk);
	}
	ConfFinal(); //Write final configuration
	
	temp = temp - 0.1;
	do{
		Restart();
		for(int iblk=1; iblk <= nblk; ++iblk){
			Reset(iblk);
			for(int istep=1; istep <= nstep; ++istep){
				Move();
				Measure();
				Accumulate();
			}
			Averages(iblk);
		}
		ConfFinal();
		temp = temp - 0.1;
	}while(temp > 0.4);

	return 0;
}


void Input(void) {
	ifstream ReadInput;

	cout << endl << "Classic 1D Ising model" << endl;
	cout << "Monte Carlo simulation" << endl;
	cout << "Nearest neighbour interaction" << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T" << endl;
	cout << "The program uses k_B=1 and mu_B=1 units" << endl << endl;

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
	cout << "The program perform Metropolis moves" << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();


	//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility

	n_props = 4; //Number of observables

	iup = 0;
	idown = 1;
	iup_h = 2;
	idown_h = 3;

	delta_ene[iup] = exp(4 * beta * J);
	delta_ene[idown] = exp(- 4 * beta * J);
	delta_ene[iup_h] = exp(2 * beta * h);
	delta_ene[idown_h] = exp(- 2 * beta * h);

	//initial configuration
	for (int i=0; i<nspin; ++i) {
		if(rnd.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
	cout << "Turning off the interaction" << endl << endl;
	turn_on = false;

	Measure();

	//Print initial values for the potential energy and virial
	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
	cout << "Initial heat capacity = " << walker[ic]/(double)nspin << endl;
	cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
	cout << "Initial susceptibility = " << walker[ix]/(double)nspin << endl << endl;
}


void Restart(void) {
	
	//initial configuration
	ifstream ReadConf("config.metro");
	for (int i=0; i<nspin; ++i) {
		ReadConf >> s[i];
	}
	ReadConf.close();

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
	delta_ene[iup] = exp(4 * beta * J);
	delta_ene[idown] = exp(- 4 * beta * J);
	if(turn_on == true){
		delta_ene[iup_h] = exp(2 * beta * h);
		delta_ene[idown_h] = exp(- 2 * beta * h);
	}	
}


void Move(void) {
	int o;
	double r, weight, snew;
	
	o = (int)(rnd.Rannyu()*nspin);
	snew = -1.0 * s[o];
	weight = Boltzmann(s[o],o);
	if(turn_on == true){
		if(snew - s[o] == 2)
			weight = weight * delta_ene[iup_h];
		else if(snew - s[o] == -2)
			weight = weight * delta_ene[idown_h];
	}
	r = rnd.Rannyu();
	if(r < weight){
		s[o] = snew;
		accepted++;
	}
	attempted++;
}


double Boltzmann(double sm, int ip) {
	double nn_ene = sm * (s[Pbc(ip+1)] + s[Pbc(ip-1)]);
	if(nn_ene == 2) return delta_ene[idown];
	else if(nn_ene == 0) return 1;
	else return delta_ene[iup];
}


void Measure() {
	double u = 0;
	double sum = 0;

	if(turn_on == false){
		for (int i=0; i<nspin; ++i) {
			u += - J * s[i] * s[Pbc(i+1)];
			sum += s[i];
		}
		walker[iu] = u;
		walker[ic] = u * u;
		walker[ix] = beta * sum * sum;
		walker[im] = sum;
	}
	else if(turn_on == true){
		for (int i=0; i<nspin; ++i) {
			u += - J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
			sum += s[i];
		}
		walker[im] = sum;
	}
}


void Reset(int iblk) {

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


void Accumulate(void) {

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) {

	if(turn_on == false){
		ofstream WriteEne, WriteHeat, WriteChi;

		stima_c = blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm ,2);
		stima_c = beta * beta * stima_c;

		stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
		glob_av[iu]  += stima_u;
		glob_av2[iu] += stima_u*stima_u;
		err_u=Error(glob_av[iu],glob_av2[iu],iblk);

		stima_c = stima_c/(double)nspin; //Heat capacity
		glob_av[ic]  += stima_c;
		glob_av2[ic] += stima_c*stima_c;
		err_c=Error(glob_av[ic],glob_av2[ic],iblk);

		stima_x = blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
		glob_av[ix]  += stima_x;
		glob_av2[ix] += stima_x*stima_x;
		err_x=Error(glob_av[ix],glob_av2[ix],iblk);

		//Writing measurments of one block
		WriteEne.open("ene.metro",ios::app);
		WriteHeat.open("heat.metro",ios::app);
		WriteChi.open("chi.metro",ios::app);

		WriteEne << iblk << "   " << glob_av[iu]/(double)iblk << "   " << err_u << endl;
		WriteHeat << iblk << "   " << glob_av[ic]/(double)iblk << "   " << err_c << endl;
		WriteChi << iblk << "   " << glob_av[ix]/(double)iblk << "   " << err_x << endl;

		WriteEne.close();
		WriteHeat.close();
		WriteChi.close();
	}
	else if(turn_on == true){
		ofstream WriteMag;

		stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
		glob_av[im]  += stima_m;
		glob_av2[im] += stima_m*stima_m;
		err_m=Error(glob_av[im],glob_av2[im],iblk);

		WriteMag.open("mag.metro",ios::app);

		WriteMag << iblk << "   " << glob_av[im]/(double)iblk << "   " << err_m << endl;

		WriteMag.close();
	}
}


void ConfFinal(void) {
	ofstream WriteConf;

	WriteConf.open("config.metro");
	for (int i=0; i<nspin; ++i) {
		WriteConf << s[i] << endl;
	}
	WriteConf.close();

	if(turn_on == false){//Writing final results
		ofstream WriteEne, WriteHeat, WriteChi;

		WriteEne.open("eneVStemp.metro",ios::app);
		WriteHeat.open("heatVStemp.metro",ios::app);
		WriteChi.open("chiVStemp.metro",ios::app);

		WriteEne << temp << "   ";
		WriteEne << glob_av[iu]/(double)nblk << "   ";
		WriteEne << err_u << endl;

		WriteHeat << temp << "   ";
		WriteHeat << glob_av[ic]/(double)nblk << "   ";
		WriteHeat << err_c << endl;

		WriteChi << temp << "   ";
		WriteChi << glob_av[ix]/(double)nblk << "   ";
		WriteChi << err_c << endl;

		WriteEne.close();
		WriteHeat.close();
		WriteChi.close();
	}
	else if(turn_on == true){
		ofstream WriteMag;

		WriteMag.open("magVStemp.metro",ios::app);

		WriteMag << temp << "   ";
		WriteMag << glob_av[im]/(double)nblk << "   ";
		WriteMag << err_m << endl;

		WriteMag.close();
	}

	rnd.SaveSeed();
}


int Pbc(int i) {
	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;
	return i;
}


double Error(double sum, double sum2, int iblk) {
	if(iblk == 1) return 0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
