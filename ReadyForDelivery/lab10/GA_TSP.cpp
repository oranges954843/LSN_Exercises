#include <iostream>
#include <fstream>
#include <cmath>

#include "GA_TSP.h"
#include "mpi.h"

using namespace std;
int size, rank;

int main(int argc, char* argv[]){

	
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat1, stat2, stat3, stat4;
	MPI_Request req;

	int* imesg1 = new int[L];
	int* imesg2 = new int[L]; 
	int* imesg3 = new int[L];
	int* imesg4 = new int[L];
	
	int itag1=0, itag2=1,itag3=2, itag4=3;

	if(size != 4){
		cout << "need 4 processes not " << size << endl;
		return 1;
	}
	
	Creation();


	for(int igen = 1; igen <= ngen; igen++){
		
		if(igen%50 == 0)
			cout << "I'm " << rank << " at the " << igen << "-th generation" << endl;

		new_n = 0;
		for(int istep = 1; istep <= N/2; istep++){
			Selection();
			Crossover();
			Mutation();

			for(int l = 0; l < L; l++){
				new_population_xy[new_n][l] = offspring_father_xy[l];
				new_population_xy[new_n+1][l] = offspring_mother_xy[l];
			}
			new_n = new_n + 2;
		}
		
		for(int n = 0; n < N; n++){
			for(int l = 0; l < L; l++){
				population_xy[n][l] = new_population_xy[n][l];
			}
		}

		Check();
		Cost();
		Averages(igen);
		if((igen%Nmigr == 0) && (igen%5 == 0)){
			for(int imigr = 0; imigr < 25; imigr++){
				if(rank == 1){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg1[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg1[0],L,MPI_INTEGER,0,itag1,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg2[0],L,MPI_INTEGER,0,itag2,MPI_COMM_WORLD,&stat2);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg2[l];
				}
				else if(rank == 0){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg2[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg2[0],L,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg1[0],L,MPI_INTEGER,1,itag1,MPI_COMM_WORLD,&stat1);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg1[l];
				}
				if(rank == 3){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg3[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg3[0],L,MPI_INTEGER,2,itag3,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg4[0],L,MPI_INTEGER,2,itag4,MPI_COMM_WORLD,&stat4);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg4[l];
				}
				else if(rank == 2){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg4[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg4[0],L,MPI_INTEGER,3,itag4,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg3[0],L,MPI_INTEGER,3,itag3,MPI_COMM_WORLD,&stat3);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg3[l];
				}
			}
		}
		if((igen%Nmigr == 0) && (igen%10 == 0)){
			for(int imigr = 0; imigr < 35; imigr++){
				if(rank == 1){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg1[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg1[0],L,MPI_INTEGER,2,itag1,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg2[0],L,MPI_INTEGER,2,itag2,MPI_COMM_WORLD,&stat2);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg2[l];
				}
				else if(rank == 2){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg2[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg2[0],L,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg1[0],L,MPI_INTEGER,1,itag1,MPI_COMM_WORLD,&stat1);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg1[l];
				}
				if(rank == 3){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg3[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg3[0],L,MPI_INTEGER,0,itag3,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg4[0],L,MPI_INTEGER,0,itag4,MPI_COMM_WORLD,&stat4);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg4[l];
				}
				else if(rank == 0){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg4[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg4[0],L,MPI_INTEGER,3,itag4,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg3[0],L,MPI_INTEGER,3,itag3,MPI_COMM_WORLD,&stat3);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg3[l];
				}
			}	
		}
		if((igen%Nmigr == 0) && (igen%15 != 0)){
			for(int imigr = 0; imigr < 45; imigr++){
				if(rank == 1){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg1[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg1[0],L,MPI_INTEGER,3,itag1,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg2[0],L,MPI_INTEGER,3,itag2,MPI_COMM_WORLD,&stat2);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg2[l];
				}
				else if(rank == 3){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg2[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg2[0],L,MPI_INTEGER,1,itag2,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg1[0],L,MPI_INTEGER,1,itag1,MPI_COMM_WORLD,&stat1);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg1[l];
				}
				if(rank == 2){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg3[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg3[0],L,MPI_INTEGER,0,itag3,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg4[0],L,MPI_INTEGER,0,itag4,MPI_COMM_WORLD,&stat4);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg4[l];
				}
				else if(rank == 0){
					casual_n = int(pow(rnd.Rannyu(),4)*N);
					for(int l = 0; l < L; l++) imesg4[l] = population_xy[order[casual_n]][l];
					MPI_Isend(&imesg4[0],L,MPI_INTEGER,2,itag4,MPI_COMM_WORLD,&req); 
					MPI_Recv(&imesg3[0],L,MPI_INTEGER,2,itag3,MPI_COMM_WORLD,&stat3);
					casual_n = int(pow(rnd.Rannyu(),0.2)*N);
					for(int l = 0; l < L; l++) population_xy[order[casual_n]][l] = imesg3[l];
				}
			}	
		}
	}

	ConfFinal();

	
	MPI_Finalize();
	return 0;
}


void Creation(void){

	int p[8] = {0};
        ifstream Primes("Primes");
	int p1, p2;
	Primes >> p1 >> p2;
	for(int i = 0; i < 8; i++)
		Primes >> p[i];

	ifstream input("seed.in");
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);

        Primes.close(); 
        input.close();

	for(int i = 0; i < L; i++){
		for(int j = 0; j < 2; j++)
			position_xy[i][j] = rnd.Rannyu();
		cities_xy[i] = i;
	}

	if(rank == 0)
		rnd.SetRandom(seed,p[0],p[1]);
	if(rank == 1)
		rnd.SetRandom(seed,p[2],p[3]);
	if(rank == 2)
		rnd.SetRandom(seed,p[4],p[5]);
	if(rank == 3)
		rnd.SetRandom(seed,p[6],p[7]);

	for(int n = 0; n < N; n++){
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

	double sum_xy = 0;

	lowest_nxy = 0, highest_nxy = 0;
	for(int n = 0; n < N; n++){

		F_xy[n] = 0;
		for(int l = 0; l < L; l++){

			for(int j = 0; j < 2; j++)
				F_xy[n] += pow(position_xy[population_xy[n][l]][j] - position_xy[population_xy[n][Pbc(l+1)]][j],2);
		}
		if(F_xy[n] <= F_xy[lowest_nxy])
			lowest_nxy = n;
		if(F_xy[n] > F_xy[highest_nxy])
			highest_nxy = n;

		sum_xy += 1.0 / F_xy[n];
	}

	for(int n = 0; n < N; n++){
		
		fitness_xy[n] = 1.0 / F_xy[n] / sum_xy;
	}

}


void Check(void){
	
	int storage_xy[L];
	for(int n = 0; n < N; n++){
		for(int l = 0; l < L; l++){
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
	double partition_unity_xy[N+1] = {0};
	double sum_xy = 0;

	for(int n = 0; n < N; n++){

		sum_xy += fitness_xy[n];
		partition_unity_xy[n+1] = sum_xy;
	}
	
	f = rnd.Rannyu();
	m = rnd.Rannyu();
			
	for(int n = 0; n < N; n++){
		if((f >= partition_unity_xy[n]) && (f < partition_unity_xy[n+1])) 
			n_father_xy = n;
		if((m >= partition_unity_xy[n]) && (m < partition_unity_xy[n+1])) 
			n_mother_xy = n;
	}
}


void Crossover(void){
	
	double Pc = 0.4;
	int cut, k_f = 0, k_m = 0;
	
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
			Shuffle(offspring_father_xy,l);
			Shuffle(offspring_mother_xy,l);
		}
	}
}


void Averages(int igen){
	
	ofstream WriteFXY, WriteFXYBest;
	if(rank == 0){
		WriteFXY.open("ave_CostFunction_MPI_1.out",ios::app);
		WriteFXYBest.open("best_CostFunction_MPI_1.out",ios::app);
	}
	else if(rank == 1){
		WriteFXY.open("ave_CostFunction_MPI_2.out",ios::app);
		WriteFXYBest.open("best_CostFunction_MPI_2.out",ios::app);
	}
	else if(rank == 2){
		WriteFXY.open("ave_CostFunction_MPI_3.out",ios::app);
		WriteFXYBest.open("best_CostFunction_MPI_3.out",ios::app);
	}
	else{
		WriteFXY.open("ave_CostFunction_MPI_4.out",ios::app);
		WriteFXYBest.open("best_CostFunction_MPI_4.out",ios::app);
	}

	double ave = 0;
	double appo_xy, min_xy = 0;
		
	
	for(int n = 0; n < N; n++) order[n] = n;
	for(int n1 = 0; n1 < N; n1++){
	
		min_xy = F_xy[n1];
		lowest_nxy = order[n1];
		for(int n2 = n1+1; n2 < N; n2++){
			if(F_xy[n2] < min_xy){
				appo_xy = F_xy[n2];
				F_xy[n2] = min_xy;
				min_xy = appo_xy;
				order[n1] = order[n2];
				order[n2] = lowest_nxy;
				lowest_nxy = order[n1];
			}
		}
		order[n1] = lowest_nxy;
		F_xy[n1] = min_xy;
	}	

	for(int n = 0; n < N/2; n++)
		ave += F_xy[n];
	ave = ave / (double)(0.5 * N);

	WriteFXY << igen << "   " << ave << endl;
	WriteFXYBest << igen << "   " << F_xy[0] << endl;

	WriteFXY.close();
	WriteFXYBest.close();
}


void ConfFinal(void){

	ofstream WritePathXY("best_path_MPI.out",ios::app);
	
	for(int l = 0; l < L; l++){
	
		WritePathXY << population_xy[0][l] << "   ";
		for(int j = 0; j < 2; j++)
			WritePathXY << position_xy[population_xy[0][l]][j] << "  ";
		WritePathXY << endl;
	}

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
