#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;

#define M 10000 //# of realizations
#define decay_rate 1
#define mu 0
#define gamma 1

int main(int argc, char *argv[]) {


int N[4] = {1,2,10,100};
double sum_standard, sum_exponential, sum_lorentzian;
Random rnd;
int seed[4];
int p1, p2;
ofstream WriteStan, WriteExpo, WriteLore;

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



/*standard, exponential and lorentzian dice
with N=1*/
WriteStan.open("standard_N1.out");
WriteExpo.open("exponential_N1.out");
WriteLore.open("lorentzian_N1.out");
for (int i=0; i<M; i++) {
	sum_standard = 0;
	sum_exponential = 0;
	sum_lorentzian = 0;	
	for (int j=0; j<N[0]; j++) {
		sum_standard += rnd.Rannyu();
		sum_exponential += rnd.Exp(decay_rate);
		sum_lorentzian += rnd.Cauchy_Lorentz(mu,gamma);	
	}
	WriteStan << sum_standard/N[0] << endl;
	WriteExpo << sum_exponential/N[0] << endl;
	WriteLore << sum_lorentzian/N[0] << endl;
}
WriteStan.close();
WriteExpo.close();
WriteLore.close();



/*standard, exponential and lorentzian dice
with N=2*/
WriteStan.open("standard_N2.out");
WriteExpo.open("exponential_N2.out");
WriteLore.open("lorentzian_N2.out");
for (int i=0; i<M; i++) {
	sum_standard = 0;
	sum_exponential = 0;
	sum_lorentzian = 0;
	for (int j=1; j<N[1]; j++) {
		sum_standard += rnd.Rannyu();
		sum_exponential += rnd.Exp(decay_rate);
		sum_lorentzian += rnd.Cauchy_Lorentz(mu,gamma);
	}
	WriteStan << sum_standard/N[1] << endl;
	WriteExpo << sum_exponential/N[1] << endl;
	WriteLore << sum_lorentzian/N[1] << endl;
}
WriteStan.close();
WriteExpo.close();
WriteLore.close();



/*standard, exponential and lorentzian dice
with N=10*/
WriteStan.open("standard_N10.out");
WriteExpo.open("exponential_N10.out");
WriteLore.open("lorentzian_N10.out");
for (int i=0; i<M; i++) {
	sum_standard = 0;
	sum_exponential = 0;
	sum_lorentzian = 0;
	for (int j=0; j<N[2]; j++) {
		sum_standard += rnd.Rannyu();
		sum_exponential += rnd.Exp(decay_rate);
		sum_lorentzian += rnd.Cauchy_Lorentz(mu,gamma);
	}
	WriteStan << sum_standard/N[2] << endl;
	WriteExpo << sum_exponential/N[2] << endl;
	WriteLore << sum_lorentzian/N[2] << endl;
}
WriteStan.close();
WriteExpo.close();
WriteLore.close();



/*standard, exponential and lorentzian dice
with N=100*/
WriteStan.open("standard_N100.out");
WriteExpo.open("exponential_N100.out");
WriteLore.open("lorentzian_N100.out");
for (int i=0; i<M; i++) {
	sum_standard = 0;
	sum_exponential = 0;
	sum_lorentzian = 0;
	for (int j=0; j<N[3]; j++) {
		sum_standard += rnd.Rannyu();
		sum_exponential += rnd.Exp(decay_rate);
		sum_lorentzian += rnd.Cauchy_Lorentz(mu,gamma);
	}
	WriteStan << sum_standard/N[3] << "   " << (i+1)*N[3] << endl;
	WriteExpo << sum_exponential/N[3] << endl;
	WriteLore << sum_lorentzian/N[3] << endl;
}
WriteStan.close();
WriteExpo.close();
WriteLore.close();

rnd.SaveSeed();

return 0;

}
