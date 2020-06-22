#ifndef _FunzioneBase_h_
#define _FunzioneBase_h_


class FunzioneBase {

	public:
		virtual double Eval(double x) const = 0;
		virtual double EvalNonU(double x) const = 0; 
};	

#endif
