#include "integrale.h"
#include "random.h"
#include "integranda.h"
#include "FunzioneBase.h"

#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

#define N 60
#define M 12000
#define nstep 100

int main(int argc, char *argv[]) {


FunzioneBase * f = new integranda();
Integrale * I = new Integrale(0,1,f);
double sum;
double ave[N] = {0};
double ave2[N] = {0};
double ave_prog[N] = {0};
double ave2_prog[N] = {0};
int L = M/N;
ofstream WriteUnif, WriteNonUnif; 


//uniform distribution
for(int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) { //L estimates of the integral
		sum += I->Media(nstep);
	}
	ave[i] = sum/L;
	ave2[i] = ave[i]*ave[i];
}
WriteUnif.open("Uniform.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++)  {
		ave_prog[i] += ave[j];
		ave2_prog[i] += ave2[j];
	}
	ave_prog[i] = ave_prog[i]/(i+1);
	ave2_prog[i] = ave2_prog[i]/(i+1);
	if (i != 0) {
		WriteUnif << i*L << "   ";
		WriteUnif << ave_prog[i] << "   ";
		WriteUnif << sqrt((ave2_prog[i]-pow(ave_prog[i],2))/i) << endl;
	} else {
		WriteUnif << i*L << "   ";
		WriteUnif << ave_prog[i] << "   ";
		WriteUnif << 0 << endl;
	}
}
WriteUnif.close();


//cleaning vectors
for (int i=0; i<N; i++) {
	ave_prog[i] = 0;
	ave2_prog[i] = 0;
}


//non-uniform distribution
for(int i=0; i<N; i++) {
	sum = 0;
	for (int j=0; j<L; j++) {
		sum += I->MediaNonU(nstep);
	}
	ave[i] = sum/L;
	ave2[i] = ave[i]*ave[i];
}
WriteNonUnif.open("NonUniform.out");
for (int i=0; i<N; i++) {
	for (int j=0; j<i+1; j++)  {
		ave_prog[i] += ave[j];
		ave2_prog[i] += ave2[j];
	}
	ave_prog[i] = ave_prog[i]/(i+1);
	ave2_prog[i] = ave2_prog[i]/(i+1);
	if (i != 0) {
		WriteNonUnif << i*L << "   ";
		WriteNonUnif << ave_prog[i] << "   ";
		WriteNonUnif << sqrt((ave2_prog[i]-pow(ave_prog[i],2))/i) << endl;
	} else {
		WriteNonUnif << i*L << "   ";
		WriteNonUnif << ave_prog[i] << "   ";
		WriteNonUnif << 0 << endl;
	}
}
WriteNonUnif.close();

delete f;
delete I;

return 0;

}
	
	
