#include "random.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

#define N 100 // number of blocks
#define M 10000 //number of random walks
#define a 1.0 //length of a step in the lattice
#define T 100 //time steps

int main(int argc, char *argv[]) {


double ave[N] = {0};
double ave2[N] = {0};
double ave_prog = 0;
double ave2_prog = 0;
int L = M/N;
double theta, phi;
double sum, coin; //it tells me which sign along a specific direction
Random rnd;
double x, y, z;
ofstream WriteZ3, WriteR3;
int seed[4];
int p1, p2;

ifstream Primes("Primes");
if (Primes.is_open()){
	Primes >> p1 >> p2 ;
} else cerr << "PROBLEM: Unable to open Primes" << endl;
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
} else cerr << "PROBLEM: Unable to open seed.in" << endl;


//Z3 cubic lattice
WriteZ3.open("Z3_lattice.out");
//initial position: (0,0,0)
WriteZ3 << 0 << "   ";
WriteZ3 << 0 << "   ";
WriteZ3 << 0 << endl;
for (int t=1; t<=T; t++) {
	ave_prog = 0;
	ave2_prog = 0;
	for (int i=0; i<N; i++) { //every cicle is a simulation of M random walks with t time-steps diveded in N blocks
		sum = 0;
		for (int j=0; j<L; j++) {
			x = 0;
			y = 0;
			z = 0;
			for (int n=1; n<=t; n++) { //un random walk di t+1 passi
				coin = rnd.Rannyu(0,1); //sign of the displacement along x
				if (coin<0.5) {
					x += 2 * a * 0.5; 
				} else x -= 2 * a * 0.5;
				coin = rnd.Rannyu(0,1); //sign of the displacement along y
				if (coin<0.5) {
					y += 2 * a * 0.5;
				} else y -= 2 * a * 0.5;
				coin = rnd.Rannyu(0,1); //sign of the displacement along z
				if (coin<0.5) {
					z += 2 * a * 0.5;
				} else z -= 2 * a * 0.5;
			}
			sum += sqrt(x*x+y*y+z*z);
		}
		ave[i] = sum/L;
		ave2[i] = ave[i]*ave[i];
	}
	for (int i=0; i<N; i++)  {
		ave_prog += ave[i];
		ave2_prog += ave2[i];
	}
	ave_prog = ave_prog/N;
	ave2_prog = ave2_prog/N;
	WriteZ3 << t << "   ";
	WriteZ3 << ave_prog << "   ";
	WriteZ3 << sqrt((ave2_prog - pow(ave_prog,2))/(N-1)) << endl;
}
WriteZ3.close();
////////////////////////
////////////////////////


//R3 cubic lattice
WriteR3.open("R3_lattice.out");
//initial position: (0,0,0)
WriteR3 << 0 << "   ";
WriteR3 << 0 << "   ";
WriteR3 << 0 << endl;
for (int t=1; t<=T; t++) {
	ave_prog = 0;
	ave2_prog = 0;
	for (int i=0; i<N; i++) { //every cicle is a simulation of M random walks with t time-steps diveded in N blocks
		sum = 0;
		for (int j=0; j<L; j++) {
			x = 0;
			y = 0;
			z = 0;
			for (int n=1; n<=t; n++) { //un random walk di t+1 passi
				theta = rnd.Rannyu(0,M_PI);
				phi = rnd.Rannyu(0,2*M_PI);
				x += a * sin(theta) * cos(phi);
				y += a * sin(theta) * cos(phi);
				z += a * sin(theta) * cos(phi);
				
			}
			sum += sqrt(x*x+y*y+z*z);
		}
		ave[i] = sum/L;
		ave2[i] = ave[i]*ave[i];
	}
	for (int i=0; i<N; i++)  {
		ave_prog += ave[i];
		ave2_prog += ave2[i];
	}
	ave_prog = ave_prog/N;
	ave2_prog = ave2_prog/N;
	WriteR3 << t << "   ";
	WriteR3 << ave_prog << "   ";
	WriteR3 << sqrt((ave2_prog - pow(ave_prog,2))/(N-1)) << endl;
}
WriteR3.close();
////////////////////////
////////////////////////


rnd.SaveSeed();

return 0;

}
