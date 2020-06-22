#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>

#include "main.h"

using namespace std;

#define a0 1.0
#define N 100
#define M 1000000
#define eqstep 2000

int main() {
	
	ofstream WriteR100, WriteR210;
	WriteR100.open("r100_eq.out");
	WriteR210.open("r210_eq.out");

	Input();
	cout << "Equilibration....." << endl << endl;
	for(int i=0; i<eqstep; i++){ //Equlibration
		Move();
		WriteR100 << i << "  " << walker[iGS] << endl;
		WriteR210 << i << "  " << walker[i2P] << endl;
	}
	WriteR100.close();
	WriteR210.close();
	cout << "Equilibration ended" << endl << endl;
	cout << "Acceptance rate for GS " << acc_GS/att_GS << endl;
	cout << "Acceptance rate for 2P " << acc_2P/att_2P << endl << endl;

	cout << "Simulation....." << endl << endl;
	for (int iblk=1; iblk<=nblk; iblk++) {
		Reset(iblk);
		for (int istep=0; istep<nstep; istep++) {
			Move();
			Accumulate();
		}
		Averages(iblk);
	}
	rnd.SaveSeed();
	cout << "Simulation ended" << endl << endl;

	return 0;
}


void Input(void){
	
	cout << endl << "This program implements the Metropolis algorithm" << endl;
	cout << "applied to the estimation of the expectation values of the radius" << endl;
	cout << "in the ground state and in one of the three 2P excited states" << endl;
	cout << "of the Hydrogen Atom" << endl << endl;

	nblk = N;
	nstep = (int)M/N;

	cout << "Number of blocks " << nblk << endl;
	cout << "Number of steps in a block " << nstep << endl;

	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();

	iGS = 0;
	i2P = 1;

	xGS = 0;
	yGS = 0;
	zGS = 0;

	x2P = 1;
	y2P = 1;
	z2P = 1;

	cout << "Starting point for the ground state (";
	cout <<  xGS << "," << yGS << "," << zGS << ")" << endl;
	cout << "Starting point for the excited state (";
	cout <<  x2P << "," << y2P << "," << z2P << ")" << endl << endl;

}


void Move(void){
	double probGS, prob2P;
	double xold, yold, zold, xnew, ynew, znew;
	double rold, rnew;

	//Ground State
	xold = xGS;
	yold = yGS;
	zold = zGS;

	rold = sqrt(xold*xold + yold*yold + zold*zold);
	walker[iGS] = rold;

	xnew = xold + rnd.Rannyu(-1.2*a0,1.2*a0);
	ynew = yold + rnd.Rannyu(-1.2*a0,1.2*a0);
	znew = zold + rnd.Rannyu(-1.2*a0,1.2*a0);

	rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);

	probGS = PsiQGS(rnew)/PsiQGS(rold);
	//cout << "probGS = " << probGS << endl;
	if (probGS >= rnd.Rannyu()) {
	
		xGS = xnew;
		yGS = ynew;
		zGS = znew;
			
		walker[iGS] = rnew;
		acc_GS = acc_GS + 1;
	}
	att_GS = att_GS + 1;

	//Excited State
	xold = x2P;
	yold = y2P;
	zold = z2P;

	rold = sqrt(xold*xold + yold*yold + zold*zold);
	walker[i2P] = rold;

	xnew = xold + rnd.Rannyu(-3*a0,3*a0);
	ynew = yold + rnd.Rannyu(-3*a0,3*a0);
	znew = zold + rnd.Rannyu(-3*a0,3*a0);

	rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);

	prob2P = PsiQ2P(rnew,znew/rnew)/PsiQ2P(rold,zold/rold);
	//cout << "prob2P = " << prob2P << endl;
	if (prob2P >= rnd.Rannyu()) {
	
		x2P = xnew;
		y2P = ynew;
		z2P = znew;
			
		walker[i2P] = rnew;
		acc_2P = acc_2P + 1;
	}
	att_2P = att_2P + 1;
}


double PsiQGS(double r) {
	return exp(-2*r/a0);
}


double PsiQ2P(double r, double cos_t) {
	return r*r*exp(-r/a0)*pow(cos_t,2);
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
	att_GS = 0;
	acc_GS = 0;
	att_2P = 0;
	acc_2P = 0;
}


void Accumulate(void){

	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){

	ofstream WriteGS, Write2P;
    	
	WriteGS.open("output_GS.dat",ios::app);
	Write2P.open("output_2P.dat",ios::app);
    
	stima_GS = blk_av[iGS]/blk_norm;
	glob_av[iGS] += stima_GS;
	glob_av2[iGS] += stima_GS*stima_GS;
	err_GS = Error(glob_av[iGS],glob_av2[iGS],iblk);

	stima_2P = blk_av[i2P]/blk_norm;
	glob_av[i2P] += stima_2P;
	glob_av2[i2P] += stima_2P*stima_2P;
	err_2P = Error(glob_av[i2P],glob_av2[i2P],iblk);

	WriteGS << iblk << "   " << glob_av[iGS]/(double)iblk;
	WriteGS << "   " << err_GS << endl;
	Write2P << iblk << "   " << glob_av[i2P]/(double)iblk;
	Write2P << "   " << err_2P << endl;


	WriteGS.close();
	Write2P.close();
}


double Error(double sum, double sum2, int iblk){
	if( iblk == 1 ) return 0.0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}	 
