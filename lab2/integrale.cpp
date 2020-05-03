#include "integrale.h"
#include "FunzioneBase.h"
#include "random.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

Integrale::Integrale(double a, double b, FunzioneBase * f) {
	
_integranda = f;
_a = min(a,b);
_b = max(a,b);
_generatore = new Random();
int seed[4];
int p1, p2;

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
                        _generatore->SetRandom(seed,p1,p2);
                }
        }
        input.close();
} else cerr << "PROBLEM: Unable to open seed.in" << endl;


if (a > b) _sign = -1;
else _sign = 1;

}

Integrale::~Integrale() { 
_generatore->SaveSeed();
delete _generatore; }


double Integrale::Media(int nstep) {
	
double x;
double sum = 0;
for (int i = 0; i<nstep; i++) {
	x = _generatore->Rannyu(_a,_b);
	sum += _integranda->Eval(x);
}

return (sum/nstep)*(_b -_a);
}


double Integrale::MediaNonU(int nstep) {
	
double x;
double sum = 0;
for (int i = 0; i<nstep; i++) {
	x = _generatore->NonU();
	sum += _integranda->EvalNonU(x);
}

return (sum/nstep)*(_b -_a);
}
	

