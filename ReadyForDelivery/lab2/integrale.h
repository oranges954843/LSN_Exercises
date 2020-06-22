#ifndef _Integrale_H_
#define _Integrale_H_

#include "FunzioneBase.h"
#include "random.h"

class Integrale {
		
	public:
		Integrale(double a, double b, FunzioneBase * f);
		~Integrale();
		double Media(int nstep);
		double MediaNonU(int nstep);
		//double HitorMiss(int nstep, double M, double m);

	private:
		double _a, _b;
		double _sum;
		double _h;
		double _integrale;
		FunzioneBase * _integranda;
		int _sign;
		Random * _generatore;
};

#endif
