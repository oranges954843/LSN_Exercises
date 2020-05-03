/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

#define nblocks 100
#define M 4000 

using namespace std;

int main(){
 
ofstream WriteEtot, WriteEpot, WriteEkin, WriteTemp;
int nconf = 1;
int L = M/nblocks;
cout << "nblocks = " << nblocks << endl;
cout << "M = " << M << endl;
cout << "L = " << L << endl;
int k = 0;
double v_epot[M] = {0}, v_ekin[M] = {0}, v_etot[M] = {0}, v_temp[M] = {0};
double ave_etot[nblocks] = {0}, ave_epot[nblocks] = {0}, ave_ekin[nblocks] = {0}, ave_temp[nblocks] = {0};
double ave2_etot[nblocks] = {0}, ave2_epot[nblocks] = {0}, ave2_ekin[nblocks] = {0}, ave2_temp[nblocks] = {0};
double ave_prog_etot[nblocks] = {0}, ave_prog_epot[nblocks] = {0}, ave_prog_ekin[nblocks] = {0}, ave_prog_temp[nblocks] = {0};
double ave2_prog_etot[nblocks] = {0}, ave2_prog_epot[nblocks] = {0}, ave2_prog_ekin[nblocks] = {0}, ave2_prog_temp[nblocks] = {0};
double sum_etot, sum_epot, sum_ekin, sum_temp;

Input();             //Inizialization
for(int istep=1; istep <= nstep; ++istep){
	Move();           //Move particles with Verlet algorithm
	if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
	if(istep%10 == 0){
        	Measure();     //Properties measurement
        	ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        	nconf += 1;
		v_etot[k] = stima_etot;
		v_epot[k] = stima_pot;
		v_ekin[k] = stima_kin;
		v_temp[k] = stima_temp;
		k++;
     	}
}
ConfFinal();         //Write final configuration to restart


cout << " k = " << k << endl;

for(int i=0; i<3; i++) {
	cout << endl << "///////////////////////////////" << endl;
	cout << "Simulazione " << i+2 << endl;
	NewInput();
      	for(int istep=1; istep <= nstep; ++istep){
		Move();           //Move particles with Verlet algorithm
        	if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        	if(istep%10 == 0){
          		Measure();     //Properties measurement
          		ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          		nconf += 1;
			v_epot[k] = stima_pot;
			v_ekin[k] = stima_kin;
			v_etot[k] = stima_etot;
			v_temp[k] = stima_temp;
			k++;
        	}

	}
    	NewConfFinal();         //Write final configuration to restart
}


//block averaging
for (int i=0; i<nblocks; i++){
	sum_etot = 0;
	sum_epot = 0;
	sum_ekin = 0;
	sum_temp = 0;
	for (int j=0; j<L; j++){
		sum_etot += v_etot[j+i*L];
		sum_epot += v_epot[j+i*L];
		sum_ekin += v_ekin[j+i*L];
		sum_temp += v_temp[j+i*L];
	}
	ave_etot[i] = sum_etot/L;
	ave_epot[i] = sum_epot/L;
	ave_ekin[i] = sum_ekin/L;
	ave_temp[i] = sum_temp/L;
	ave2_etot[i] = ave_etot[i]*ave_etot[i];
	ave2_epot[i] = ave_epot[i]*ave_epot[i];
	ave2_ekin[i] = ave_ekin[i]*ave_ekin[i];
	ave2_temp[i] = ave_temp[i]*ave_temp[i];
}

WriteEtot.open("ave_etot.out");
WriteEpot.open("ave_epot.out");
WriteEkin.open("ave_ekin.out");
WriteTemp.open("ave_temp.out");
for (int i=0; i<nblocks; i++){
	for (int j=0; j<i+1; j++){
		ave_prog_etot[i] += ave_etot[j];
		ave_prog_epot[i] += ave_epot[j];
		ave_prog_ekin[i] += ave_ekin[j];
		ave_prog_temp[i] += ave_temp[j];
		ave2_prog_etot[i] += ave2_etot[j];
		ave2_prog_epot[i] += ave2_epot[j];
		ave2_prog_ekin[i] += ave2_ekin[j];
		ave2_prog_temp[i] += ave2_temp[j];
	}
	ave_prog_etot[i] /= (i+1);
	ave_prog_epot[i] /= (i+1);
	ave_prog_ekin[i] /= (i+1);
	ave_prog_temp[i] /= (i+1);
	ave2_prog_etot[i] /= (i+1);
	ave2_prog_epot[i] /= (i+1);
	ave2_prog_ekin[i] /= (i+1);
	ave2_prog_temp[i] /= (i+1);
	WriteEtot << ave_prog_etot[i] << "   " << sqrt(ave2_prog_etot[i] - pow(ave_prog_etot[i],2)) << endl; 
	WriteEpot << ave_prog_epot[i] << "   " << sqrt(ave2_prog_epot[i] - pow(ave_prog_epot[i],2)) << endl;
	WriteEkin << ave_prog_ekin[i] << "   " << sqrt(ave2_prog_ekin[i] - pow(ave_prog_ekin[i],2)) << endl;
	WriteTemp << ave_prog_temp[i] << "   " << sqrt(ave2_prog_temp[i] - pow(ave_prog_temp[i],2)) << endl;
}
WriteEtot.close();
WriteEpot.close();
WriteEkin.close();
WriteTemp.close();


return 0;

}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void NewInput(void){ //Prepare all stuff for the simulation
  ifstream ReadConf,ReadOldConf;
  double t,newfs, xnew, ynew, znew;
  double fx[m_part], fy[m_part], fz[m_part];

  cout << endl;
  cout << "An other iteration " << endl;

//Read initial configuration
  cout << "Read initial configuration from file old.0 and old.final " << endl << endl;
  ReadConf.open("old.0"); //r(t)
  ReadOldConf.open("old.final"); //r(t-dt)
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    ReadOldConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();
  ReadOldConf.close();


  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }
  t = 0.0;
  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - x[i])/(delta);
    vy[i] = Pbc(ynew - y[i])/(delta);
    vz[i] = Pbc(znew - z[i])/(delta);

    t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  }
  t /= (double)npart;
  newfs = sqrt(3 * temp / (2 * t));

  for (int i=0; i<npart; ++i){

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2 * delta);
    vz[i] = Pbc(znew - zold[i])/(2 * delta);

    vx[i] *= newfs;
    vy[i] *= newfs;
    vz[i] *= newfs;
   
    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
    
    xold[i] = Pbc(x[i] - delta * vx[i]);
    yold[i] = Pbc(y[i] - delta * vy[i]);
    zold[i] = Pbc(z[i] - delta * vz[i]);
  }

  return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void ConfFinal(void){ //Write final configuration
  ofstream WriteConf,WriteOldConf,WriteOld0;

  cout << endl;
  cout << "Print final configuration to file config.final " << endl;
  WriteConf.open("config.final");
  cout << "Print final configuration to file old.final and old.0 " << endl << endl;
  WriteOldConf.open("old.final");
  WriteOld0.open("old.0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOld0 << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteOldConf.close();
  WriteOld0.close();

  return;
}
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
void NewConfFinal(void){ //Write final configuration
  ofstream WriteOldConf,WriteOld0;

  cout << endl;
  cout << "Print final configuration to file old.final and old.0 " << endl << endl;
  WriteOldConf.open("old.final");
  WriteOld0.open("old.0");
  
  for (int i=0; i<npart; ++i){
    WriteOld0 << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteOldConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteOldConf.close();
  WriteOld0.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
