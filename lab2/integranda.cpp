#include "integranda.h"
#include <iostream>
#include <cmath>

using namespace std;


double integranda::Eval(double x) const {
	return 0.5*M_PI*cos(0.5*M_PI*x);
}
double integranda::EvalNonU(double x) const {
	if (x == 0) {
		return 0;
	} 
	else { 
		return (M_PI*M_PI)*x*sin(0.5*M_PI*x)/(12*x*x);
	}
}

