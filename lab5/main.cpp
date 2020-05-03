#include "random.h"
#include "main.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstring>

using namespace std;

#define a0 1.0
#define N 200
#define M 1000000

int main() {

double xnew, ynew, znew;
double x = 1, y = 1, z = 1;
double a, r, rnew, sum;
double r_ave[N] = {0}, r_ave2[N] = {0}, r_ave_prog[N] = {0}, r_ave2_prog[N] = {0};
int seed[4];
int accepted, rejected, p1, p2;
int L = M/N;
ofstream WriteR100, WriteR210;
Random rnd;

ifstream Primes("Primes");
if (Primes.is_open()){
        Primes >> p1 >> p2 ;
} else cerr << "PROBLEM: Unable to open Primes" << endl;
Primes.close();
ifstream input("seed.in");
string property;
if (input.is_open()){
        while ( !input.eof() ){
                input >> property;
                if( property == "RANDOMSEED" ){
                        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                        rnd.SetRandom(seed,p1,p2);
                }
        }
        input.close();
} else cerr << "PROBLEM: Unable to open seed.in" << endl;


///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
//sampling psi_r100
accepted = 0;
rejected = 0;
x = 3;
y = 4;
z = 5;
do { //equilibration of the system
	xnew = rnd.Rannyu(x-a0,x+a0);
	ynew = rnd.Rannyu(y-a0,y+a0);
	znew = rnd.Rannyu(z-a0,z+a0);
	r = sqrt(x*x + y*y + z*z);
	rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
	a = minA(PsiQ100(rnew),PsiQ100(r));
	if (rnd.Rannyu() < a) {
		x = xnew;
		y = ynew;
		z = znew;
		accepted++;
	} else rejected ++;
} while (rejected < 100);
cout << "accepted = " << accepted << endl;
cout << "rejected = " << rejected << endl;

for (int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) {
		xnew = rnd.Rannyu(x-a0,x+a0);
		ynew = rnd.Rannyu(y-a0,y+a0);
		znew = rnd.Rannyu(z-a0,z+a0);
		r = sqrt(x*x + y*y + z*z);
		rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
		a = minA(PsiQ100(rnew),PsiQ100(r));
		if (rnd.Rannyu() < a) {
			x = xnew;
			y = ynew;
			z = znew;
			sum += rnew;
		} else sum += r;
	}
	r_ave[i] = sum/L;
	r_ave2[i] = r_ave[i]*r_ave[i];
}
WriteR100.open("r100_ave.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++) {
		r_ave_prog[i] += r_ave[j];
		r_ave2_prog[i] += r_ave2[j];
	}
	r_ave_prog[i] /= (i+1); 
	r_ave2_prog[i] /= (i+1);
	if (i != 0) {
		WriteR100 << i*L << "   ";
		WriteR100 << r_ave_prog[i] << "   ";
		WriteR100 << sqrt((r_ave2_prog[i] - pow(r_ave_prog[i],2))/i) << endl;
	} else {
		WriteR100 << i*L << "   ";
		WriteR100 << r_ave_prog[i] << "   ";
		WriteR100 << 0 << endl;
	}
}
WriteR100.close();

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
//sampling psi_r210
for (int i=0; i<N; i++)  {
	r_ave_prog[i] = 0;
	r_ave2_prog[i] = 0;
}
accepted = 0;
rejected = 0;
x = 3;
y = 4;
z = 5;
do { //equilibration of the system
	xnew = rnd.Rannyu(x-a0,x+a0);
	ynew = rnd.Rannyu(y-a0,y+a0);
	znew = rnd.Rannyu(z-a0,z+a0);
	r = sqrt(x*x + y*y + z*z);
	rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
	a = minA(PsiQ210(rnew,znew/rnew),PsiQ210(r,z/r));
	if (rnd.Rannyu() < a) {
		x = xnew;
		y = ynew;
		z = znew;
		accepted++;
	} else rejected ++;
} while (rejected < 1000);
cout << "accepted r210 = " << accepted << endl;
cout << "rejected r210 = " << rejected << endl;

for (int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) {
		xnew = rnd.Rannyu(x-a0,x+a0);
		ynew = rnd.Rannyu(y-a0,y+a0);
		znew = rnd.Rannyu(z-a0,z+a0);
		r = sqrt(x*x + y*y + z*z);
		rnew = sqrt(xnew*xnew + ynew*ynew + znew*znew);
		a = minA(PsiQ210(rnew,znew/rnew),PsiQ210(r,z/r));
		if (rnd.Rannyu() < a) {
			x = xnew;
			y = ynew;
			z = znew;
			sum += rnew;
		} else sum += r;
	}
	r_ave[i] = sum/L;
	r_ave2[i] = r_ave[i]*r_ave[i];
}
WriteR210.open("r210_ave.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++) {
		r_ave_prog[i] += r_ave[j];
		r_ave2_prog[i] += r_ave2[j];
	}
	r_ave_prog[i] /= (i+1); 
	r_ave2_prog[i] /= (i+1);
	if (i != 0) {
		WriteR210 << i*L << "   ";
		WriteR210 << r_ave_prog[i] << "   ";
		WriteR210 << sqrt((r_ave2_prog[i] - pow(r_ave_prog[i],2))/i) << endl;
	} else {
		WriteR210 << i*L << "   ";
		WriteR210 << r_ave_prog[i] << "   ";
		WriteR210 << 0 << endl;
	}
}
WriteR210.close();

rnd.SaveSeed();

return 0;

}



/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double minA(double PsiNew, double Psi) {
	if (PsiNew < Psi) {
		return PsiNew/Psi;
	} else return 1.0;
}

double PsiQ100(double r) {
	return exp(-2*r/a0);
}

double PsiQ210(double r, double cos_t) {
	return r*r*exp(-r/a0)*pow(cos_t,2);
}
	 
