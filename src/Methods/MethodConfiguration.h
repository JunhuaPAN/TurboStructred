#ifndef TurboStructured_Methods_MethodConfiguration
#define TurboStructured_Methods_MethodConfiguration

#include "EquationsOfState/EquationsOfState.h"

//Method configuration class
struct MethodConfiguration {
	int RungeKuttaOrder;
	double CFL;
	double Eps;
	double OperatingPresure;

	//for general eos cases
	bool UseExactPressureDerivative;
	GeneralEOS* eos;
};

#endif
