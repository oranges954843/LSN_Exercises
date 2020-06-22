#include "random.h"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;

#define N 20 // number of blocks
#define M 10000 //number of random walks
#define a 1.0 //length of a step in the lattice
#define T 100 //time steps

int main(int argc, char *argv[]) {


double ave, ave2;
double theta, phi;
double sum, p;
int L = M / N;
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
WriteZ3 << 0 << "   " << 0 << "   " << 0 << endl;
for (int t=1; t<=T; t++) {
	ave = 0, ave2 = 0;
	for (int i=0; i<N; i++){ //every cicle is a simulation of L random walks with t time-steps 

		sum = 0;
		for(int j=0; j<L; j++){

			x = 0;
			y = 0;
			z = 0;
			for (int n=1; n<=t; n++) {

				p = rnd.Rannyu(0,6);
				if(p < 1){
					x += 2 * a * 0.5; 
				} 
				else if((p >= 1)&&(p < 2)){
					x -= 2 * a * 0.5;
				}
				else if((p >= 2)&&(p < 3)){
					y += 2 * a * 0.5;
				}
				else if((p >= 3)&&(p < 4)){
					y -= 2 * a * 0.5;
				}
				else if((p >= 4)&&(p < 5)){
					z += 2 * a * 0.5;
				}
				else{
					z -= 2 * a * 0.5;
				}
			}
			sum += x*x+y*y+z*z;
		}
		ave += sum/L;
		ave2 += pow(sum/L,2);
	}
	ave = sqrt(ave / N);
	ave2 = sqrt(ave2 / N);
	WriteZ3 << t << "   " << ave << "   ";
	WriteZ3 << sqrt((ave2 - pow(ave,2))/(N-1)) << endl;
}
WriteZ3.close();


//R3 cubic lattice
WriteR3.open("R3_lattice.out");
WriteR3 << 0 << "   " << 0 << "   " << 0 << endl;
for(int t=1; t<=T; t++){
	ave = 0, ave2 = 0;
	for(int i=0; i<N; i++){//every cicle is a simulation of L random walks with t time-steps 
		sum = 0;
		for(int j=0; j<L; j++){

			x = 0;
			y = 0;
			z = 0;
			for (int n=1; n<=t; n++) {

				theta = rnd.Rannyu(0,M_PI);
				phi = rnd.Rannyu(0,2*M_PI);
				x += a * sin(theta) * cos(phi);
				y += a * sin(theta) * sin(phi);
				z += a * cos(theta);
			}
			sum += x*x+y*y+z*z;
		}
		ave += sum/L;
		ave2 += pow(sum/L,2);
	}
	ave = sqrt(ave / N);
	ave2 = sqrt(ave2 / N);
	WriteR3 << t << "   " << ave << "   ";
	WriteR3 << sqrt((ave2 - pow(ave,2))/(N-1)) << endl;
}
WriteR3.close();


rnd.SaveSeed();

return 0;

}
