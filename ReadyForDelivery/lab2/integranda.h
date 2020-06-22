#ifndef _INTEGRANDA_H_
#define _INTEGRANDA_H_

#include "FunzioneBase.h"

class integranda : public FunzioneBase {

	public:
		double Eval(double x) const;
		double EvalNonU(double x) const;

};


#endif 
