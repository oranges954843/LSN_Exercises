#include <iostream>
#include <fstream>
#include <cmath>

#include "GA_TSP.h"
using namespace std;

int main(){

	Creation();

	for(int igen = 1; igen <= ngen; igen++){

		if(igen%50 == 0)
			cout << igen << "-th generation" << endl;

		new_n = 0;
		for(int istep = 1; istep <= N/2; istep++){
			Selection();
			Crossover();
			Mutation();

			for(int l = 0; l < L; l++){
				new_population_x[new_n][l] = offspring_father_x[l];
				new_population_x[new_n+1][l] = offspring_mother_x[l];
				new_population_xy[new_n][l] = offspring_father_xy[l];
				new_population_xy[new_n+1][l] = offspring_mother_xy[l];
			}
			new_n = new_n + 2;
		}
		
		for(int n = 0; n < N; n++){
			for(int l = 0; l < L; l++){
				population_x[n][l] = new_population_x[n][l];
				population_xy[n][l] = new_population_xy[n][l];
			}
		}

		Check();
		Cost();
		Averages(igen);
	}

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

	for(int n = 0; n < N; n++){
		population_x[n][0] = cities_x[0];
		population_xy[n][0] = cities_xy[0];
		for(int l = 1; l < L; l++){
			Shuffle(cities_x,l);
			Shuffle(cities_xy,l);
			population_x[n][l] = cities_x[l];
			population_xy[n][l] = cities_xy[l];
		}
	}

	Check();
	Cost();
}


void Cost(void){

	double sum_x = 0, sum_xy = 0;

	for(int n = 0; n < N; n++){
		F_x[n] = 0;
		F_xy[n] = 0;
		for(int l = 0; l < L; l++){
			F_x[n] += pow(position_x[population_x[n][l]] - position_x[population_x[n][Pbc(l+1)]],2);
			for(int j = 0; j < 2; j++)
				F_xy[n] += pow(position_xy[population_xy[n][l]][j] - position_xy[population_xy[n][Pbc(l+1)]][j],2);
		}
		sum_x += 1.0 / F_x[n];
		sum_xy += 1.0 / F_xy[n];
	}

	double sumP_x = 0, sumP_xy = 0;
	double max = 1.0 - 1E-8;

	for(int n = 0; n < N; n++){
		fitness_x[n] = 1.0 / F_x[n] / sum_x;
		fitness_xy[n] = 1.0 / F_xy[n] / sum_xy;
		sumP_x += fitness_x[n];
		sumP_xy += fitness_xy[n];
	}

	if(sumP_x < max){
		cerr << "Error!! Probability is not conserved on the circle" << endl;
		cout << "Total probability = " << sumP_x << endl;
	}
	if(sumP_xy < max){
		cerr << "Error!! Probability is not conserved on the square" << endl;
		cout << "Total probability = " << sumP_xy << endl;
	}
}


void Check(void){
	
	int storage_x[L], storage_xy[L];
	for(int n = 0; n < N; n++){
		for(int l = 0; l < L; l++){
			if((l == 0) && (population_x[n][l] != cities_x[0]))
				cerr << "Error on the circle!! The first city is not 1 " << endl;
			if((l == L-1) && (population_x[n][Pbc(l+1)] != cities_x[0]))
				cerr << "Error on the circle!! The p.b.c are not satisfied " << endl;
			storage_x[l] = population_x[n][l];
			for(int i = 0; i <= l; i++){
				if((storage_x[i] == storage_x[l]) && (i != l))
					cerr << "Error on the circle!! There is a copy " << endl;
			}
			if((l == 0) && (population_xy[n][l] != cities_xy[0]))
				cerr << "Error on the square!! The first city is not 1 " << endl;
			if((l == L-1) && (population_xy[n][Pbc(l+1)] != cities_xy[0]))
				cerr << "Error on the square!! The p.b.c are not satisfied " << endl;
			storage_xy[l] = population_xy[n][l];
			for(int i = 0; i <= l; i++){
				if((storage_xy[i] == storage_xy[l]) && (i != l))
					cerr << "Error on the square!! There is a copy " << endl;
			}
		}
	}
}


void Selection(void){

	double f, m;
	double partition_unity_x[N+1] = {0}, partition_unity_xy[N+1] = {0};
	double sum_x = 0, sum_xy = 0;

	for(int n = 0; n < N; n++){
		sum_x += fitness_x[n];
		sum_xy += fitness_xy[n];
		partition_unity_x[n+1] = sum_x;
		partition_unity_xy[n+1] = sum_xy;
	}
	
	f = rnd.Rannyu();
	m = rnd.Rannyu();
			
	for(int n = 0; n < N; n++){
		if((f >= partition_unity_x[n]) && (f < partition_unity_x[n+1])) 
			n_father_x = n;
		if((m >= partition_unity_x[n]) && (m < partition_unity_x[n+1])) 
			n_mother_x = n;
	}
	for(int n = 0; n < N; n++){
		if((f >= partition_unity_xy[n]) && (f < partition_unity_xy[n+1])) 
			n_father_xy = n;
		if((m >= partition_unity_xy[n]) && (m < partition_unity_xy[n+1])) 
			n_mother_xy = n;
	}
}


void Crossover(void){
	
	double Pc = 0.4;
	int cut, k_f = 0, k_m = 0;;

	if(Pc > rnd.Rannyu()){
		cut = (int)(rnd.Rannyu()*(L-1)) + 1;
		int father_end_x[L-cut], mother_end_x[L-cut];

		for(int l = 0; l < L; l++){
			father[l] = population_x[n_father_x][l];
			mother[l] = population_x[n_mother_x][l];
			offspring_father_x[l] = population_x[n_father_x][l];
			offspring_mother_x[l] = population_x[n_mother_x][l];
		}

		for(int l = 0; l < L-cut; l++){
			father_end_x[l] = population_x[n_father_x][l+cut];
			mother_end_x[l] = population_x[n_mother_x][l+cut];
		}

		for(int l = 1; l < L; l++){
			for(int i = 0; i < L-cut; i++){
				if(mother[l] == father_end_x[i]){
					offspring_father_x[k_f+cut] = mother[l];
					k_f++;
				}
				if(father[l] == mother_end_x[i]){
					offspring_mother_x[k_m+cut] = father[l];
					k_m++;
				}
			}
		}
	}
	else{
		for(int l = 0; l < L; l++){
			offspring_father_x[l] = population_x[n_father_x][l];
			offspring_mother_x[l] = population_x[n_mother_x][l];
		}
	}
	
	if(Pc > rnd.Rannyu()){
		cut = (int)(rnd.Rannyu()*(L-1)) + 1;
		int father_end_xy[L-cut], mother_end_xy[L-cut];

		for(int l = 0; l < L; l++){
			father[l] = population_xy[n_father_xy][l];
			mother[l] = population_xy[n_mother_xy][l];
			offspring_father_xy[l] = population_xy[n_father_xy][l];
			offspring_mother_xy[l] = population_xy[n_mother_xy][l];
		}
		for(int l = 0; l < L-cut; l++){
			father_end_xy[l] = population_xy[n_father_xy][l+cut];
			mother_end_xy[l] = population_xy[n_mother_xy][l+cut];
		}

		k_f = 0, k_m = 0;
		for(int l = 1; l < L; l++){
			for(int i = 0; i < L-cut; i++){
				if(mother[l] == father_end_xy[i]){
					offspring_father_xy[k_f+cut] = mother[l];
					k_f++;
				}
				if(father[l] == mother_end_xy[i]){
					offspring_mother_xy[k_m+cut] = father[l];
					k_m++;
				}
			}
		}
	}
	else{
		for(int l = 0; l < L; l++){
			offspring_father_xy[l] = population_xy[n_father_xy][l];
			offspring_mother_xy[l] = population_xy[n_mother_xy][l];
		}
	}
}
 

void Mutation(void){
	
	double Pm = 0.08;

	for(int l = 1; l < L; l++){
		if(Pm > rnd.Rannyu()){
			Shuffle(offspring_father_x,l);
			Shuffle(offspring_mother_x,l);
		}
		if(Pm > rnd.Rannyu()){
			Shuffle(offspring_father_xy,l);
			Shuffle(offspring_mother_xy,l);
		}
	}
}


void Averages(int igen){
	
	ofstream WriteFX("ave_CostFunction_circle.out",ios::app);
	ofstream WriteFXY("ave_CostFunction_square.out",ios::app);
	ofstream WriteFXBest("best_CostFunction_circle.out",ios::app);
	ofstream WriteFXYBest("best_CostFunction_square.out",ios::app);

	double ave = 0;
	double min_x, min_xy, appo_x, appo_xy;

	for(int n1 = 0; n1 < N; n1++){
		min_x = F_x[n1];
		min_xy = F_xy[n1];
		for(int n2 = n1+1; n2 < N; n2++){
			if(F_x[n2] < min_x){
				appo_x = F_x[n2];
				F_x[n2] = min_x;
				min_x = appo_x;
			}
			if(F_xy[n2] < min_xy){
				appo_xy = F_xy[n2];
				F_xy[n2] = min_xy;
				min_xy = appo_xy;
			}
		}
		F_x[n1] = min_x;
		F_xy[n1] = min_xy;
	}	

	for(int n = 0; n < N/2; n++)
		ave += F_x[n];
	ave = ave / (double)(0.5 * N);

	WriteFX << igen << "   " << ave << endl;
	WriteFXBest << igen << "   " << F_x[0] << endl;

	ave = 0;
	for(int n = 0; n < N/2; n++)
		ave += F_xy[n];
	ave = ave / (double)(0.5 * N);

	WriteFXY << igen << "   " << ave << endl;
	WriteFXYBest << igen << "   " << F_xy[0] << endl;

	WriteFX.close();
	WriteFXY.close();
	WriteFXBest.close();
	WriteFXYBest.close();
}


void ConfFinal(void){

	ofstream WritePathX("best_path_circle.out");
	ofstream WritePathXY("best_path_square.out");
	
	for(int l = 0; l < L; l++){
		WritePathX << population_x[0][l] << "   ";
		WritePathX << position_x[population_x[0][l]] << endl;
		WritePathXY << population_xy[0][l] << "   ";
		for(int j = 0; j < 2; j++)
			WritePathXY << position_xy[population_xy[0][l]][j] << "  ";
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
