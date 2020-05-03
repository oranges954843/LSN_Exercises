#include "random.h"
#include <iostream>
#include <cmath>
#include <fstream> 
#include <cstring>


using namespace std;

#define N 100
#define M 10000
#define r 0.1 //risk-free interest rate
#define K 100.0 //strike price
#define sigma 0.25 //volatility
#define S0 100.0 //initial asset price
#define T 1.0 //delivery time

int main(int argc, char *argv[]) {


int L = M/N;
Random rnd;
double S[100] = {0};
double ST;
double PutAve[N] = {0}, CallAve[N] = {0}; //put and call option price
double PutAve2[N] = {0}, CallAve2[N] = {0};
double PutAveProg[N] = {0}, CallAveProg[N] = {0};
double PutAve2Prog[N] = {0}, CallAve2Prog[N] = {0}; 
int seed[4];
int p1, p2;
ofstream WriteCall, WritePut;


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



S[0] = K; //strike price 100.0
for (int i=0; i<N; i++) {
	for (int j=0; j<L; j++) {
		for (int k=1; k<100; k++) {
			S[k] = S[k-1]*exp((r-sigma*sigma*0.5)*0.01+sigma*rnd.Gauss(0,1)*sqrt(0.01));
		}
		ST = S[99];  
		if (ST == K) {
			PutAve[i] += 0.0;
			CallAve[i] += 0.0;
		}
		else {
			if (ST > K) {
				CallAve[i] += exp(-r*T)*(ST-K);
				PutAve[i] += 0.0;
			} 
			else {
				CallAve[i] += 0.0;
				PutAve[i] += exp(-r*T)*(K-ST);
			}
		}
	}
	CallAve[i] = CallAve[i]/L;
	PutAve[i] = PutAve[i]/L;
	CallAve2[i] = CallAve[i]*CallAve[i];
	PutAve2[i] = PutAve[i]*PutAve[i];
}
WriteCall.open("call_GBM.out");
WritePut.open("put_GBM.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++)  {
			PutAveProg[i] += PutAve[j];
			CallAveProg[i] += CallAve[j];
			PutAve2Prog[i] += PutAve2[j];
			CallAve2Prog[i] += CallAve2[j];
		}
	CallAveProg[i] = CallAveProg[i]/(i+1);
	PutAveProg[i] = PutAveProg[i]/(i+1);
	CallAve2Prog[i] = CallAve2Prog[i]/(i+1);
	PutAve2Prog[i] = PutAve2Prog[i]/(i+1);
	if (i != 0) {
		WriteCall << i*L << "   ";
		WriteCall << CallAveProg[i] << "   ";
		WriteCall << sqrt((CallAve2Prog[i] - pow(CallAveProg[i],2))/i) << endl;
		WritePut << i*L << "   ";
		WritePut << PutAveProg[i] << "   ";
		WritePut << sqrt((PutAve2Prog[i] - pow(PutAveProg[i],2))/i) << endl;
	} else {
		WriteCall << i*L << "   ";
		WriteCall << CallAveProg[i] << "   ";
		WriteCall << 0 << endl;
		WritePut << i*L << "   ";
		WritePut << PutAveProg[i] << "   ";
		WritePut << 0 << endl;
	}
}
WriteCall.close();
WritePut.close();

rnd.SaveSeed();


	
return 0;

}
