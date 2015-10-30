#ifndef TurboStructured_Utility_Integration
#define TurboStructured_Utility_Integration

//Functions operating geometrical objects

#include <cmath>
#include <string>
#include <iostream>
#include <vector>
#include "Vector.h"

// compute integral of real value function f(x) - 1D
double ComputeIntegrall(const std::valarray<double>& f, const std::valarray<double>& hx) {
	if ((f.size() == 0) || (hx.size() == 0)) return 0;
	double res = f[0] * hx[0] + f[f.size() - 1] * hx[hx.size() - 1];
	for (int i = 1; i < hx.size(); i++) {
		res += 0.5 * f[i] * (0.5 * hx[i + 1] + 0.5 * hx[i - 1] + hx[i]);
	};

	return res;
};


#endif
