#ifndef TurboStructured_utility_NumericQuadrature
#define TurboStructured_utility_NumericQuadrature

#include <valarray>
#include <vector>
#include "Vector.h"
#include "grid.h"


// class for numerical integration over the cell by Gause quadrature formula

struct IntegrInfo {
	// points number
	int nx{ 1 };
	int ny{ 1 };
	int nz{ 1 };
};

class NumericQuadrature{
private:
	// integration order
	int n;

	// dimension number
	int Dims;

	// points number
	IntegrInfo Pn;

	// Jacobian value
	inline double ComputeJacobian(CellInfo& cell, int nDims) {
		double J = 0.5 * cell.hx;
		if(nDims > 1) J *= 0.5 * cell.hy;
		if (nDims > 2) J *= 0.5 * cell.hz;

		return J;
	};

	// weights for 1D integration
	std::valarray<double> wx{ 1 };
	std::valarray<double> wy{ 1 };
	std::valarray<double> wz{ 1 };

	// Legandr's polinome roots
	std::valarray<double> x{ 0 };
	std::valarray<double> y{ 0 };
	std::valarray<double> z{ 0 };

	// fill weights and Legandr's polinome roots
	void FillCoefficients();

public:
	// constructor
	NumericQuadrature(int _n, int _Dims) : n(_n), Dims(_Dims) {
		Pn.nx = _n;
		if (_Dims > 1) Pn.ny = _n;
		if (_Dims > 2) Pn.nz = _n;
		FillCoefficients();
	};

	// main function
	template< typename F >
	double Integrate(CellInfo& cell, F f ) {

		// integration result
		double res = 0;
		double rx, ry, rz;
		for (int k = 0; k < Pn.nz; k++) {
			for (int j = 0; j < Pn.ny; j++) {
				for (int i = 0; i < Pn.nx; i++) {
					rx = cell.x + x[i] * 0.5 * cell.hx;
					ry = cell.y + y[i] * 0.5 * cell.hy;
					rz = cell.z + z[i] * 0.5 * cell.hz;
					res += wx[i] * wy[j] * wz[k] * f(Vector(rx, ry, rz));
				};
			};
		};

		// multiply by jacobian
		res *= ComputeJacobian(cell, Dims);

		return res;
	};

	template< typename F >
	std::valarray<double> IntegrateGroap(CellInfo& cell, F f, size_t nvars) {

		// integration result
		std::valarray<double> res( 0.0, nvars);
		double rx, ry, rz;
		double w;

		// Gauss's integration method
		for (int k = 0; k < Pn.nz; k++) {
			for (int j = 0; j < Pn.ny; j++) {
				for (int i = 0; i < Pn.nx; i++) {
					// Compute exact variablse in a Gauss point
					rx = cell.x + x[i] * 0.5 * cell.hx;
					ry = cell.y + y[i] * 0.5 * cell.hy;
					rz = cell.z + z[i] * 0.5 * cell.hz;
					auto vars = f(Vector(rx, ry, rz));

					// Compute total weight
					w = wx[i] * wy[j] * wz[k];

					// Calculate contribution of given point
					for (int nv = 0; nv < nvars; nv++) {
						res[nv] += w * vars[nv];
					};
				};
			};
		};

		// multiply by jacobian
		res *= ComputeJacobian(cell, Dims);

		return res;
	};

};

void NumericQuadrature::FillCoefficients() {
	// create two arrays for data containing
	std::valarray<double> r;	// Legandre's roots
	std::valarray<double> w;	// weights

	// create Legandre's roots of n-th degree polinome and weght coefficients
	switch (n) {
	// odd number first
	case 1:
		r = { 0 };
		w = { 1 };
		break;
	case 3:
		/* n = 3 */
		r = { 0.0000000000000000000000000,0.7745966692414833770358531 };
		w = { 0.8888888888888888888888889,0.5555555555555555555555556 };
		break;
	case 5:
		/* n = 5 */
		r = { 0.0000000000000000000000000,0.5384693101056830910363144,0.9061798459386639927976269 };
		w = { 0.5688888888888888888888889,0.4786286704993664680412915,0.2369268850561890875142640 };
		break;
	case 7:
		r = { 0.0000000000000000000000000,0.4058451513773971669066064,0.7415311855993944398638648,0.9491079123427585245261897 };
		w = { 0.4179591836734693877551020,0.3818300505051189449503698,0.2797053914892766679014678,0.1294849661688696932706114 };
		break;

	// even number of points
	case 2:
		/* n = 2 */
		r = { 0.5773502691896257645091488 };
		w = { 1.0000000000000000000000000 };
		break;
	case 4:
		/* n = 4 */
		r = { 0.3399810435848562648026658,0.8611363115940525752239465 };
		w = { 0.6521451548625461426269361,0.3478548451374538573730639 };
		break;
	case 6:
		/* n = 6 */
		r = { 0.2386191860831969086305017,0.6612093864662645136613996,0.9324695142031520278123016 };
		w = { 0.4679139345726910473898703,0.3607615730481386075698335,0.1713244923791703450402961 };
		break;
	case 8:
		/* n = 8 */
		r = { 0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609 };
		w = { 0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314 };
		break;
	default:
		std::cerr << "Bad polinome's order in numerical integration, n = " << n << std::endl;
	};

	// fill both arrays for 1D case
	wx.resize(n);
	x.resize(n);

	size_t ws{ w.size() };
	for (int i = 0; i < ws; i++) {
		wx[i] = w[ws - 1 - i];
		wx[n - 1 - i] = wx[i];
		x[n - 1 - i] = r[ws - 1 - i];
		x[i] = -x[n - 1 - i];
	};

	// repeat for other directions
	if (Dims > 1) {
		wy = wx;
		y = x;
	};
	if (Dims > 2) {
		wz = wx;
		z = x;
	};

	return;
};

#endif