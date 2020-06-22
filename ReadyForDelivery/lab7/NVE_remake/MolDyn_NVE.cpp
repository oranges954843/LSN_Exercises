#include <stdlib.h>     
#include <iostream> 
#include <fstream> 
#include <cmath>      
#include "MolDyn_NVE.h"


using namespace std;

int main(){

	do{
		Input();
		do{
			cout << "Simulation " << nsim + 1 << endl;
			for(int istep=1; istep<=eqstep; istep++){

				Move();
				if(istep%10 == 0){
					Measure();
					WriteIstantValues();
				}
			}
			cout << "Estimated temperature = " << stima[it] << endl;
			ConfFinal();
			nsim++;
			if(nsim < 5){
				Restart();
			}
		}while(nsim < 5);

		cout << "Measurments........... " << endl;
		for(int iblk=1; iblk<=nblk; iblk++){
			Reset(iblk);
			for(int istep=1; istep<=nstep; istep++){
				Move();
				if(istep%10 == 0){
        				Measure();
					Accumulate();
     				}
			}
			Averages(iblk);
		}
		cout << "End of " + phase[index_ph] + " phase experiment" << endl << endl;
		nsim = 0;
		index_ph++;

	}while(index_ph<3);

	return 0;
}



void Input(void){

	ifstream ReadInput,ReadConf;
	double sumv2;
	
	iv = 0;
	ik = 1;
	ie = 2;
	it = 3;
	n_props = 4;


	cout << "-----------------------------" << endl;
	cout << "Classic Lennard-Jones fluid ";
	cout << "in " + phase[index_ph] + " phase" << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl;
	cout << "The program uses Lennard-Jones units " << endl << endl;

	seed = 1;
	srand(seed);
  
	ReadInput.open("input.Ar_" + phase[index_ph]);

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;


	igofr = 4;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box/2.0)/(double)nbins;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput >> eqstep;

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps in a block = " << nstep << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of equilibration steps = " << eqstep << endl;
	ReadInput.close();


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
	sumv2 = 0;
	for (int i=0; i<npart; ++i){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];

		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double)npart;
	
	fs = sqrt(3 * temp / sumv2);
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;

		xold[i] = Pbc(x[i] - vx[i] * delta);
		yold[i] = Pbc(y[i] - vy[i] * delta);
		zold[i] = Pbc(z[i] - vz[i] * delta);
		
	}

	cout << "Equilibrating " << phase[index_ph] << " phase" << endl << endl;
}


void Restart(void){
		
	ifstream ReadConf,ReadOldConf;
	double sumv2, xnew, ynew, znew;
	double fx[m_part], fy[m_part], fz[m_part];

	ReadConf.open("rt.final"); //r(t)
	ReadOldConf.open("old_rt.final"); //r(t-dt)
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
	
	sumv2 = 0.0;
	double sumv[3] = {0.0, 0.0, 0.0};
	for(int i=0; i<npart; ++i){ //Verlet integration scheme
	
		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		//v(t+dt/2)
		vx[i] = Pbc(xnew - x[i])/(delta);
		vy[i] = Pbc(ynew - y[i])/(delta);
		vz[i] = Pbc(znew - z[i])/(delta);
		
		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;

	}
	for(int idim=0; idim<3; idim++) sumv[idim] /= (double)npart;
	for(int i=0; i<npart; i++){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double)npart;
	fs = sqrt(3.0 * temp / sumv2); //new scale factor

	for(int idim=0; idim<3; idim++) sumv[idim] = 0.0;
	for (int i=0; i<npart; ++i){
		//v(t)
		vx[i] = Pbc(x[i] - xold[i])/(2 * delta);
		vy[i] = Pbc(y[i] - yold[i])/(2 * delta);
		vz[i] = Pbc(z[i] - zold[i])/(2 * delta);

		sumv[0] += vx[i];
		sumv[1] += vy[i];
		sumv[2] += vz[i];
	}
	for(int idim=0; idim<3; idim++) sumv[idim] /= (double)npart;
	for(int i=0; i<npart; i++){
		vx[i] = vx[i] - sumv[0];
		vy[i] = vy[i] - sumv[1];
		vz[i] = vz[i] - sumv[2];
		
		vx[i] *= fs;
		vy[i] *= fs;
		vz[i] *= fs;
    
		xold[i] = Pbc(x[i] - delta * vx[i]); //xnew(t) = x(t+dt) - vx*(t)dt
		yold[i] = Pbc(y[i] - delta * vy[i]);
		zold[i] = Pbc(z[i] - delta * vz[i]);
	}
}


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


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
    		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));
			}
		}
	}

	return f;
}


void Measure(){ //Properties measurement

	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	
	v = 0.0;
	t = 0.0;

	for(int k=igofr; k<igofr+nbins; k++) stima[k] = 0.0;
	
  	for (int i=0; i<npart-1; ++i){
	    	for (int j=i+1; j<npart; ++j){

			dx = Pbc( x[i] - x[j] );
			dy = Pbc( y[i] - y[j] );
			dz = Pbc( z[i] - z[j] );

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			bin = 4;
			for(double r=0; r < box/2.0; r += bin_size){
				if((dr >= r) && (dr < r + bin_size)){
					stima[bin] += 2;
				}
				bin++;
			}
			if(dr < rcut){
				vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
				v += vij;
			}
		}          
	}

	//Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	stima[iv] = v/(double)npart; //Potential energy per particle
	stima[ik] = t/(double)npart; //Kinetic energy per particle
	stima[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
	stima[ie] = (t+v)/(double)npart; //Total energy per particle

	return;
}


void WriteIstantValues(void){

	ofstream Epot, Ekin, Etot, Temp;

	Epot.open(phase[index_ph] + "_epot.istant",ios::app);
	Ekin.open(phase[index_ph] + "_ekin.istant",ios::app);
	Temp.open(phase[index_ph] + "_temp.istant",ios::app);
	Etot.open(phase[index_ph] + "_etot.istant",ios::app);

	Epot << stima[iv]  << endl;
	Ekin << stima[ik]  << endl;
	Temp << stima[it] << endl;
	Etot << stima[ie] << endl;

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
}


void ConfFinal(void){
	ofstream WriteConf,WriteRt,WriteOldRt;

	WriteConf.open("config_Ar_" + phase[index_ph] + ".final");
	WriteRt.open("rt.final");
	WriteOldRt.open("old_rt.final");

	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
		WriteRt << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
		WriteOldRt << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
	}
	WriteConf.close();
	WriteRt.close();
	WriteOldRt.close();
}


double Pbc(double r){ 
	return r - box * rint(r/box);
}


void Reset(int iblk) {
        if(iblk == 1) {
                for(int i=0; i<n_props; ++i) {
                        glob_av[i] = 0;
                        glob_av2[i] = 0;
                }
        }
        for(int i=0; i<n_props; ++i) {
                blk_av[i] = 0;
        }
       	blk_norm = 0;
}


void Accumulate(void) {
        for(int i=0; i<n_props; ++i) {
                blk_av[i] = blk_av[i] + stima[i];
        }
        blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) {

	ofstream WriteEtot, WriteEpot, WriteEkin, WriteTemp;
	double r, gdir, vol_bin;
	ofstream WriteGofr, WriteGave;
	
        stima_e = blk_av[ie]/blk_norm; 
        glob_av[ie]  += stima_e;
        glob_av2[ie] += stima_e*stima_e;
        err_e=Error(glob_av[ie],glob_av2[ie],iblk);

	stima_v = blk_av[iv]/blk_norm;
        glob_av[iv]  += stima_v;
        glob_av2[iv] += stima_v*stima_v;
        err_v=Error(glob_av[iv],glob_av2[iv],iblk);

	stima_k = blk_av[ik]/blk_norm;
        glob_av[ik]  += stima_k;
        glob_av2[ik] += stima_k*stima_k;
        err_k=Error(glob_av[ik],glob_av2[ik],iblk);

	stima_t = blk_av[it]/blk_norm;
        glob_av[it]  += stima_t;
        glob_av2[it] += stima_t*stima_t;
        err_t=Error(glob_av[it],glob_av2[it],iblk);

	//Writing final results
	WriteEtot.open(phase[index_ph] + "_NVE_etot.out",ios::app);
        WriteEpot.open(phase[index_ph] + "_NVE_epot.out",ios::app);
        WriteEkin.open(phase[index_ph] + "_NVE_ekin.out",ios::app);
        WriteTemp.open(phase[index_ph] + "_NVE_temp.out",ios::app);
	WriteGofr.open(phase[index_ph] + "_NVE_gofr.out",ios::app);
	WriteGave.open(phase[index_ph] + "_NVE_gave.out",ios::app);

	r = 0;
	for(int k=igofr; k < igofr + nbins; k++){
		vol_bin = (4 * pi / 3) * (pow(r+bin_size,3)-pow(r,3));
		gdir = blk_av[k]/blk_norm/(rho * (double)npart * vol_bin);
		glob_av[k] += gdir;
		glob_av2[k] += gdir*gdir;
		err_gdir = Error(glob_av[k],glob_av2[k],iblk);
		
		WriteGofr << iblk << "   " << r << "   " << r + bin_size << "   ";
		WriteGofr << gdir << endl;
		
		if(iblk == nblk){
			WriteGave << iblk << "   " << r << "   " << r + bin_size << "   ";
			WriteGave << glob_av[k]/(double)iblk << "   " << err_gdir << endl;
		}
		r = r + bin_size;
	}

	WriteEtot << iblk <<  "   ";
	WriteEtot << glob_av[ie]/(double)iblk << "   ";
	WriteEtot << err_e << endl;
	
	WriteEpot << iblk <<  "   ";
	WriteEpot << glob_av[iv]/(double)iblk << "   ";
	WriteEpot << err_v << endl;

	WriteEkin << iblk <<  "   ";
	WriteEkin << glob_av[ik]/(double)iblk << "   ";
	WriteEkin << err_k << endl;

	WriteTemp << iblk <<  "   ";
	WriteTemp << glob_av[it]/(double)iblk << "   ";
	WriteTemp << err_t << endl;

        WriteEtot.close();
        WriteEkin.close();
        WriteEpot.close();
        WriteTemp.close();
	WriteGofr.close();
	WriteGave.close();
}


double Error(double sum, double sum2, int iblk) {
        if(iblk == 1) return 0;
	else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
