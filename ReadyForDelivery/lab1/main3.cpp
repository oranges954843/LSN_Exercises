#include <iostream>
#include <fstream>
#include "random.h"
#include <cmath>

using namespace std;

#define l 0.5 //needle length 
#define d 1.0 //spacing between straight lines
#define N 400 //number of blocks
#define Nthr 40000 //number of throws
 
int main(int argc, char *argv[]) {

Random rnd;
int seed[4];
int p1, p2, Nhit;
int L = Nthr/N;
double x; //x component of the needle
double o; //distance of the origin of the needle from the straight lines
double theta; 
double myPi[N] = {0}, myPi2[N] = {0};
double myEstimate[N] = {0}, myEstimate2[N] = {0};
ofstream WritePi;


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



/*simulation of the Buffon's experiment with implementation 
of block averaging*/
for (int i=0; i<N; i++) {
	Nhit = 0;
	for (int j=0; j<L; j++) {
		theta = acos(rnd.Rannyu(1,-1)); //angle uniformly distributed in [0,PI)
		x = l * cos(rnd.Rannyu(0,2) * theta);
		o = rnd.Rannyu(0,d);
		if ((x+o >= d) || (o+x <= 0.0)) {
			Nhit++;
		} 
	}
	myEstimate[i] = (2*l*L)/(Nhit*d);
	myEstimate2[i] = myEstimate[i]*myEstimate[i];
}
WritePi.open("my_pi.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++) {
		myPi[i] += myEstimate[j];
		myPi2[i] += myEstimate2[j];
	}
	myPi[i] = myPi[i]/(i+1);
	myPi2[i] = myPi2[i]/(i+1);
	if (i == 0) {
		WritePi << i*L << "   ";
		WritePi << myPi[i] << "   ";
		WritePi << 0 << endl;
	} 
	else {
		WritePi << i*L << "   ";
		WritePi << myPi[i] << "   ";
		WritePi << sqrt((myPi2[i]-myPi[i]*myPi[i])/i) << endl;
	}
}
WritePi.close();
			
rnd.SaveSeed();

return 0;
}






