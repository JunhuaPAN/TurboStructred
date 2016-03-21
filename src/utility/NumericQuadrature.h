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
					// Compute exact variables in a Gauss point
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
	// http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/

	switch (n) {
	// odd number first
	case 1:
		r = { 0 };
		w = { 2 };
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
		/* n = 7 */
		r = { 0.0000000000000000000000000,0.4058451513773971669066064,0.7415311855993944398638648,0.9491079123427585245261897 };
		w = { 0.4179591836734693877551020,0.3818300505051189449503698,0.2797053914892766679014678,0.1294849661688696932706114 };
		break;
	case 9:
		/* n = 9 */
		r = { 0.0000000000000000000000000,0.3242534234038089290385380,0.6133714327005903973087020,0.8360311073266357942994298,0.9681602395076260898355762 };
		w = { 0.3302393550012597631645251,0.3123470770400028400686304,0.2606106964029354623187429,0.1806481606948574040584720,0.0812743883615744119718922 };
		break;
	case 11:
		/* n = 11 */
		r = { 0.0000000000000000000000000,0.2695431559523449723315320,0.5190961292068118159257257,0.7301520055740493240934163,0.8870625997680952990751578,0.9782286581460569928039380 };
		w = { 0.2729250867779006307144835,0.2628045445102466621806889,0.2331937645919904799185237,0.1862902109277342514260976,0.1255803694649046246346943,0.0556685671161736664827537 };
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
	case 10:
		/* n = 10 */
		r = { 0.1488743389816312108848260,0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640 };
		w = { 0.2955242247147528701738930,0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688 };
	case 12:
		/* n = 12 */
		r = { 0.1252334085114689154724414,0.3678314989981801937526915,0.5873179542866174472967024,0.7699026741943046870368938,0.9041172563704748566784659,0.9815606342467192506905491 };
		w = { 0.2491470458134027850005624,0.2334925365383548087608499,0.2031674267230659217490645,0.1600783285433462263346525,0.1069393259953184309602547,0.0471753363865118271946160 };
		break;
	case 32:
		/* n = 32 */
		r = { 0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811 };
		w = { 0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071 };
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