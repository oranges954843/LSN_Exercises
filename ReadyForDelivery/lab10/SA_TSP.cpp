#include <iostream>
#include <fstream>
#include <cmath>


#include "SA_TSP.h"

using namespace std;

int main(){

	
	Creation();

	do{
		cout << "-----------------------------------" << endl;
		cout << "Temperature = " << temp << endl;
		acc_x = 0, att_x = 0;
		acc_xy = 0, att_xy = 0;

		int p1, p2;
        	ifstream Primes("Primes");
		ifstream input("seed.out");
        	Primes >> p1 >> p2 ;
        	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        	rnd.SetRandom(seed,p1,p2);
		Primes.close();
		input.close();

		for(int istep = 1; istep <= MCstep; istep++){
		
			Metropolis();
			if(istep%50 == 0)
				Averages(istep);		
		}
		
		cout << "accepatance rate on the circle = " << acc_x / att_x << endl;
		cout << "accepatance rate in the square = " << acc_xy / att_xy << endl << endl;
		temp = temp - 0.05;
		beta = 1 / temp;
	}

	while(temp >= 0.05);

	ConfFinal();

	
	return 0;
}


void Creation(void){

	int p1, p2;
        ifstream Primes("Primes");
        Primes >> p1 >> p2 ;
        Primes.close();

        ifstream input("seed.in");
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed,p1,p2);
        input.close();

	for(int i = 0; i < L; i++){
		position_x[i] = rnd.Rannyu();
		cities_x[i] = i;
		for(int j = 0; j < 2; j++)
			position_xy[i][j] = rnd.Rannyu();
		cities_xy[i] = i;
	}

	population_x[0] = cities_x[0];
	population_xy[0] = cities_xy[0];
	for(int l = 1; l < L; l++){
		Shuffle(cities_x,l);
		Shuffle(cities_xy,l);
		population_x[l] = cities_x[l];
		population_xy[l] = cities_xy[l];
	}

	Check();
	Cost();

	ene_old_x = F_x;
	ene_old_xy = F_xy;
}


void Cost(void){

	F_x = 0;
	F_xy = 0;
	for(int l = 0; l < L; l++){
		F_x += pow(position_x[population_x[l]] - position_x[population_x[Pbc(l+1)]],2);
		for(int j = 0; j < 2; j++)
			F_xy += pow(position_xy[population_xy[l]][j] - position_xy[population_xy[Pbc(l+1)]][j],2);
	}
}


void Check(void){
	
	int storage_x[L], storage_xy[L];
		
	for(int l = 0; l < L; l++){
		if((l == 0) && (population_x[l] != cities_x[0]))
			cerr << "Error on the circle!! The first city is not 1 " << endl;
		if((l == L-1) && (population_x[Pbc(l+1)] != cities_x[0]))
			cerr << "Error on the circle!! The p.b.c are not satisfied " << endl;
			storage_x[l] = population_x[l];
		for(int i = 0; i <= l; i++){
			if((storage_x[i] == storage_x[l]) && (i != l))
				cerr << "Error on the circle!! There is a copy " << endl;
		}
		if((l == 0) && (population_xy[l] != cities_xy[0]))
			cerr << "Error on the square!! The first city is not 1 " << endl;
		if((l == L-1) && (population_xy[Pbc(l+1)] != cities_xy[0]))
			cerr << "Error on the square!! The p.b.c are not satisfied " << endl;
		storage_xy[l] = population_xy[l];
		for(int i = 0; i <= l; i++){
			if((storage_xy[i] == storage_xy[l]) && (i != l))
				cerr << "Error on the square!! There is a copy " << endl;
		}
	}
}


void Metropolis(void){

	
	for(int l = 0; l < L; l++){
		old_population_x[l] = population_x[l];
		old_population_xy[l] = population_xy[l];
	}

	for(int i = 0; i < 4; i++){
		Mutation(); 
	}

	Check();
	Cost();

			
	ene_new_x = F_x;
	ene_new_xy = F_xy;

	weight_x = exp( -beta * (ene_new_x - ene_old_x));
	weight_xy = exp( -beta * (ene_new_xy - ene_old_xy));

	if(weight_x > rnd.Rannyu()){
		ene_old_x = ene_new_x;
		acc_x++;
	}
	else{
		for(int l = 0; l < L; l++){
			population_x[l] = old_population_x[l];
		}
	}
	att_x++;

	if(weight_xy > rnd.Rannyu()){
		ene_old_xy = ene_new_xy;
		acc_xy++;
	}
	else{
		for(int l = 0; l < L; l++){
			population_xy[l] = old_population_xy[l];
		}
	}
	att_xy++;
}
 

void Mutation(){
	
	int l = (int)(rnd.Rannyu()*(L-1) + 1);
	Shuffle(population_x,l);
	Shuffle(population_xy,l);
}


void Averages(int istep){
	
	ofstream WriteFX("best_CostFunction_circle.out",ios::app);
	ofstream WriteFXY("best_CostFunction_square.out",ios::app);
	
	WriteFX << istep << "   " << F_x << endl;

	WriteFXY << istep << "   " << F_xy << endl;

	WriteFX.close();
	WriteFXY.close();

}
void ConfFinal(void){

	ofstream WritePathX("best_path_circle.out");
	ofstream WritePathXY("best_path_square.out");
	
	
	for(int l = 0; l < L; l++){
		WritePathX << population_x[l] << "   ";
		WritePathX << position_x[population_x[l]] << endl;
		WritePathXY << population_xy[l] << "   ";
		for(int j = 0; j < 2; j++)
			WritePathXY << position_xy[population_xy[l]][j] << "  ";
		WritePathXY << endl;
	}

	WritePathX.close();
	WritePathXY.close();

	rnd.SaveSeed();
}


void Shuffle(int * local_population, int l){
	int i, appo;
	i = (int)(rnd.Rannyu()*(L-l) + l);
	appo = local_population[l];
	local_population[l] = local_population[i];
	local_population[i] = appo;
}


int Pbc(int i){

	if(i >= L) i = i - L;
	else if(i < 0) i = i + L;
	return i;
}
