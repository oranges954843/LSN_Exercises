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
		Input();
		ConfFinal();

		Input(); //Inizialization
		for(int iblk=1; iblk <= nblk; ++iblk){ //Simluation
			Reset(iblk);
			for(int istep=1; istep <= nstep; istep++){
				Move();
				Measure();
				Accumulate();
			}
			Averages(iblk);
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
		
		cout << "----------------------------" << endl << endl;
		cout << "Classic Lennard-Jones fluid        " << endl;
		cout << "Monte Carlo simulation             " << endl << endl;
		cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
		cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
		cout << "The program uses Lennard-Jones units " << endl;
 
		//Read input informations
		ReadInput.open("input." + phase[index_ph]);

		ReadInput >> temp;
	
		beta = 1.0/temp;
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
		cout << "Tail correction for the virial = " << ptail << endl; 

		ReadInput >> delta; 
		ReadInput >> nblk; 
		ReadInput >> nstep;
		ReadInput >> eqstep; 

		cout << "The program perform Metropolis moves with uniform translations" << endl;
		cout << "Moves parameter = " << delta << endl;
		cout << "Number of blocks = " << nblk << endl;
		cout << "Number of steps in one block = " << nstep << endl;
		cout << "Number of steps to equilibrate the system = " << eqstep << endl << endl;
		ReadInput.close();


		//Prepare arrays for measurements
		iv = 0; //Potential energy
		iw = 1; //Virial
 
		n_props = 2; //Number of observables

		//measurement of g(r)
		igofr = 2;//indice di g(r) dobbiamo calcolare un istogramma
		nbins = 100;
		n_props = n_props + nbins;
		bin_size = (box/2.0)/(double)nbins;//calcoliamo g(r) da 0 a L/2

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
		Measure();

		//Print initial values for the potential energy and virial
		cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
		cout << "Initial virial (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
		cout << "Initial pressure (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
	}
}


void Move(void){//Metropolis algorithm
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

		//Metropolis test
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
	int bin;
	double v = 0.0, w = 0.0;
	double vij, wij;
	double dx, dy, dz, dr;

	//reset the hystogram of g(r)
	for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){//doppio ciclo che cicla una volta sola su tutte le coppie del sistema

			// distance i-j in pbc
			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			//update of the histogram of g(r)
			//inserire ogni volta che troviamo una coppia ad una certa distanza sommiamo 2 al bin corrispondente
			bin = 2;
			for(double r=0; r < box/2.0; r += bin_size){
				if((dr >= r) && (dr < r + bin_size)){
					walker[bin] += 2;
				}
				bin++;
			}
			if(dr < rcut){
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

				// contribution to energy and virial
				v += vij;
				w += wij;	
			}
 		}       
	}

	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
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

	double r, gdir, vol_bin;
	ofstream WriteGofr, WriteGave, WriteEpot, WritePres;
    
	cout << "Block number " << iblk << endl;
	cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
	WriteEpot.open(phase[index_ph] + "_output.epot.0",ios::app);
	WritePres.open(phase[index_ph] + "_output.pres.0",ios::app);
	WriteGofr.open(phase[index_ph] + "_output.gofr.0",ios::app);
	WriteGave.open(phase[index_ph] + "_output.gave.0",ios::app);
    
	stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy per particle
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
	stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
	glob_av[iw] += stima_pres;
	glob_av2[iw] += stima_pres*stima_pres;
	err_press=Error(glob_av[iw],glob_av2[iw],iblk);		

	WriteEpot << iblk << "   " << stima_pot << "   " << glob_av[iv]/(double)iblk << "   " << err_pot << endl;
	WritePres << iblk << "   " << stima_pres << "   " << glob_av[iw]/(double)iblk << "   " << err_press << endl;

	r = 0;
	for(int k=igofr; k < igofr + nbins; k++){
		vol_bin = (4 * pi / 3) * (pow(r+bin_size,3) - pow(r,3));
		gdir = blk_av[k]/blk_norm/(rho * (double)npart * vol_bin);
		glob_av[k] += gdir;
		glob_av2[k] += gdir*gdir;
		err_gdir = Error(glob_av[k],glob_av2[k],iblk);
		WriteGofr << iblk << "   " << r << "   " << r + bin_size << "   ";
		WriteGofr << gdir << endl;
		if(iblk == nblk){
			WriteGave << iblk << "   " << r << "   " << r + bin_size << "   ";
			WriteGave << glob_av[k]/(double)iblk << "   " << err_gdir << endl;
		}
		r = r + bin_size;
	}

	WriteEpot.close();
	WritePres.close();
	WriteGofr.close();
	WriteGave.close();
}


void ConfFinal(void){
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}


double Pbc(double r){
	return r - box * rint(r/box);
}


double Error(double sum, double sum2, int iblk){
	if( iblk == 1 ) return 0.0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
