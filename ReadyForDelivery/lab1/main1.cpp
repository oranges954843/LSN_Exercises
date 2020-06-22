#include <iostream>
#include <fstream>
#include <cmath>

#include "random.h"

using namespace std;

#define M 10000 //# of throws
#define N 50  //# of blocks

int main(int argc, char *argv[]) {


int L = M/N;
double ave[N] = {0}, ave2[N] = {0};
double ave_prog[N] = {0}, ave2_prog[N] = {0};
double r[M] = {0}; //it will contain M pseudo-random numbers
Random rnd;
int seed[4];
int p1, p2;
int obs;
double sum = 0;
ofstream WriteAve, WriteAveError, WriteSigma2, WriteChi2;



//generating M pseudo-random numbers in [0,1)
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
for (int i=0; i<M; i++) {
	r[i] = rnd.Rannyu();
}



/*computing mean value and its statistical uncertainty and writing 
the output in file ave.out*/
for (int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) {	
		sum += r[j+i*L];
	}
	ave[i] = sum/L;
	ave2[i] = ave[i]*ave[i];
}
WriteAve.open("ave.out");
if (WriteAve.is_open()) {
	for (int i=0; i<N; i++) {
		for (int j=0; j<i+1; j++) {
			ave_prog[i] += ave[j];
			ave2_prog[i] += ave2[j];
		}
		ave_prog[i] = ave_prog[i]/(i+1);
		ave2_prog[i] = ave2_prog[i]/(i+1);
		if (i == 0) {
			WriteAve << i*L << "   ";
			WriteAve << ave_prog[i] << "   ";
 			WriteAve << 0 << endl;
		} 
		else {
			WriteAve << i*L << "   ";
			WriteAve << ave_prog[i] << "   ";
			WriteAve << sqrt((ave2_prog[i]-pow(ave_prog[i],2))/i) << endl;
		}
	}
	WriteAve.close();
} else cerr << "PROBLEM: Unable to open ave.out" << endl;



//cleaning vectors
for (int i=0; i<N; i++) {
	ave[i] = 0;
	ave2[i] = 0;
	ave_prog[i] = 0;
	ave2_prog[i] = 0;
}



/*computing sigma2 and its statistical uncertainty and writing
the output in file sigma2.out*/
for (int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) {	
		sum += pow(r[j+i*L]-0.5,2);
	}
	ave[i] = sum/L;
	ave2[i] = ave[i]*ave[i];
}
WriteSigma2.open("sigma2.out");
if (WriteSigma2.is_open()) {
	for (int i=0; i<N; i++) {
		for (int j=0; j<i+1; j++) {
			ave_prog[i] += ave[j];
			ave2_prog[i] += ave2[j];
		}
		ave_prog[i] = ave_prog[i]/(i+1);
		ave2_prog[i] = ave2_prog[i]/(i+1);
		if (i == 0) {
			WriteSigma2 << i*L << "   ";
			WriteSigma2 << ave_prog[i] << "   ";
			WriteSigma2 << 0 << endl;
		} 
		else {
			WriteSigma2 << i*L << "   ";
			WriteSigma2 << ave_prog[i] << "   ";
			WriteSigma2 << sqrt((ave2_prog[i]-pow(ave_prog[i],2))/i) << endl;
		}
	}
	WriteSigma2.close();
} else cerr << "PROBLEM: Unable to open sigma2.out" << endl;



/*computing chi2 and its statistical uncertainty and writing
the output in file chi2.out*/
WriteChi2.open("chi2.out");
for (int i=0; i<100; i++) {
	sum = 0;
	for (int j=0; j<M; j++) {
		r[j] = rnd.Rannyu();
	}
	for (double x=0; x<1; x += 0.01) {
		obs = 0;
		for (int k=0; k<M; k++) {
			if (r[k] < x+0.01 && r[k] >= x) {
				obs = obs + 1;
			}
		}
		sum += pow(obs-100.,2)/100;
	}
	WriteChi2 << i+1 << "   ";
	WriteChi2 << sum << endl;
}
WriteChi2.close();


rnd.SaveSeed();

return 0;


}
