#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

#include "Monte_Carlo_NVT.h"

using namespace std;


int main(){ 

	index_ph = 0;
	do{
		Input(); //Equilibration
		ConfFinal(); //Write final configuration

		Input(); //Inizialization
		for(int istep=1; istep <= nstep*nblk; istep++){
			if(istep%100000 == 0) cout << "Number of MC steps " << istep << endl;
			Move();
			Measure();
		}
		ConfFinal(); //Write final configuration

		index_ph++;
		equil = false;

	}while(index_ph<3);

	return 0;
}


void Input(void){
	ifstream ReadInput, ReadConf;

	int p1, p2;

	if(equil == false){
		ifstream Primes("Primes");
		Primes >> p1 >> p2 ;
		Primes.close();

		ifstream input("seed.in");
		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		rnd.SetRandom(seed,p1,p2);
		input.close();
 
		//Read input informations
		ReadInput.open("input." + phase[index_ph]);

		ReadInput >> temp;
	
		beta = 1.0/temp;
		cout << "----------------------------" << endl << endl;
		cout << "Temperature = " << temp << endl;

		ReadInput >> npart;
		cout << "Number of particles = " << npart << endl;

		ReadInput >> rho;
		cout << "Density of particles = " << rho << endl;
		vol = (double)npart/rho;
		box = pow(vol,1.0/3.0);
		cout << "Volume of the simulation box = " << vol << endl;
		cout << "Edge of the simulation box = " << box << endl;

		ReadInput >> rcut;
		cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
		//Tail corrections for potential energy and pressure
		vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
		ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3)); 
		cout << "Tail correction for the potential energy = " << vtail << endl;
		cout << "Tail correction for the virial = " << ptail << endl << endl; 

		ReadInput >> delta; 
		ReadInput >> nblk; 
		ReadInput >> nstep;
		ReadInput >> eqstep; 

		cout << "The program perform Metropolis moves with uniform translations" << endl;
		cout << "Moves parameter = " << delta << endl;
		cout << "Number of MC steps = " << nstep*nblk << endl;
		cout << "Number of steps to equilibrate the system = " << eqstep << endl << endl;
		ReadInput.close();


		//Prepare arrays for measurements
		iv = 0; //Potential energy
		iw = 1; //Virial
 
		n_props = 2; //Number of observables

		//Read initial configuration
		cout << "Read initial configuration from file config.init " << endl << endl;
		ReadConf.open("config.init");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = Pbc( x[i] * box );
			y[i] = Pbc( y[i] * box );
			z[i] = Pbc( z[i] * box );

		}
		ReadConf.close();

		attempted = 0;
		accepted = 0;
		cout << "Performing equilibration......" << endl;
		for(int istep=1; istep <= eqstep; ++istep) Move();
		equil = true;
		cout << "Equilibration ended" << endl << endl;
		cout << "Acceptance rate for the " << phase[index_ph] << " phase: ";
		cout << accepted/attempted << endl;
	}
	else{
		ifstream Primes("Primes");
		Primes >> p1 >> p2 ;
		Primes.close();

		ifstream input("seed.out");
		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		rnd.SetRandom(seed,p1,p2);
		input.close();

		cout << "Read initial configuration from file config.final " << endl << endl;
		ReadConf.open("config.final");
		for (int i=0; i<npart; ++i){
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = Pbc( x[i] * box );
			y[i] = Pbc( y[i] * box );
			z[i] = Pbc( z[i] * box );

		}
		ReadConf.close();
	}
}


void Move(void){
	int o;
	double p, energy_old, energy_new;
	double xold, yold, zold, xnew, ynew, znew;


	for(int i=0; i<npart; ++i){
		o = (int)(rnd.Rannyu()*npart);

		xold = x[o];
		yold = y[o];
		zold = z[o];

		energy_old = Boltzmann(xold,yold,zold,o);
 
		if(phase[index_ph] == "solid"){
			xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
			ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
			znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );
		}
		else if(phase[index_ph] == "liquid"){
			xnew = Pbc( x[o] + delta*(rnd.Rannyu(-0.8,0.8) ));
			ynew = Pbc( y[o] + delta*(rnd.Rannyu(-0.8,0.8) ));
			znew = Pbc( z[o] + delta*(rnd.Rannyu(-0.8,0.8) ));
		}
		else{
			xnew = Pbc( x[o] + delta*(rnd.Rannyu(-50,50) ));
			ynew = Pbc( y[o] + delta*(rnd.Rannyu(-50,50) ));
			znew = Pbc( z[o] + delta*(rnd.Rannyu(-50,50) ));
		}

		energy_new = Boltzmann(xnew,ynew,znew,o);

		p = exp(beta*(energy_old-energy_new));
		if(p >= rnd.Rannyu()){
			x[o] = xnew;
			y[o] = ynew;
			z[o] = znew;
    
			accepted = accepted + 1.0;
		}
		attempted = attempted + 1.0;
	}
}


double Boltzmann(double xx, double yy, double zz, int ip){
	double ene=0.0;
	double dx, dy, dz, dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
		
			dx = Pbc(xx - x[i]);
			dy = Pbc(yy - y[i]);
			dz = Pbc(zz - z[i]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut){
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}

	return 4.0*ene;
}


void Measure(){
	double v = 0.0, w = 0.0;
	double vij, wij;
	double dx, dy, dz, dr;
	double p, u;
	ofstream WritePres, WriteEpot;


	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

				v += vij;
				w += wij;	
			}
 		}       
	}

	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
	p = rho * temp + (walker[iw] + ptail * npart) / vol;
	u = walker[iv]/(double)npart + vtail;

	WritePres.open(phase[index_ph] + "_pres.istant",ios::app);
	WriteEpot.open(phase[index_ph] + "_epot.istant",ios::app);
	WritePres << p << endl;
	WriteEpot << u << endl;
	WritePres.close();
	WriteEpot.close();
}


void ConfFinal(void){
	ofstream WriteConf;

	cout << endl << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}


double Pbc(double r){ //Algorithm for periodic boundary conditions with side L=box

	return r - box * rint(r/box);
}
