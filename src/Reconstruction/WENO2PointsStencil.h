#ifndef TurboStructured_Reconstruction_WENO2PointsStencil
#define TurboStructured_Reconstruction_WENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/IReconstruction.h>


//Reconstruction by convex combination of 2 points stencil reconstructions
class WENO2PointsStencil : public IReconstruction {
public:
	
	//constructor
	WENO2PointsStencil() {	};
	WENO2PointsStencil(CellReconstruction& _recons, int _nDims, int _nValues) :
		IReconstruction(_recons, _nDims, _nValues) {};
};

template<>
WENO2PointsStencil ComputeReconstruction<WENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim, double gamma) {
	//compute Gradients
	auto size = value.size();
	std::valarray< Vector > gradientsL(size), gradientsR(size);		//for 1D case only
	std::vector< std::valarray < double> > weights;

	// constant value
	const double eps = 1.0e-16;
	const double d0 = 2.0 / 3.0;		//constant values for UNIFORM grid 2/3
	const double d1 = 1 - d0;

	//1D case
	if (nDim == 1) {
		weights.resize(4);
		for (auto& r : weights) r.resize(size);

		for (int i = 0; i < size; i++) {
			// compute smothness coefficients
			double betta0 = 0.5 * (points[1].x - points[0].x) * (values[1][i] - value[i]) / (points[1].x - point.x);
			betta0 *= betta0;
			betta0 += eps;
			betta0 *= betta0;

			double betta1 = 0.5 * (points[1].x - points[0].x) * (values[0][i] - value[i]) / (points[0].x - point.x);
			betta1 *= betta1;
			betta1 += eps;
			betta1 *= betta1;

			// compute weights for two interpolations (for left and right)
			weights[0][i] = d0 / betta0;
			weights[1][i] = d1 / betta1;
			weights[2][i] = d1 / betta0;
			weights[3][i] = d0 / betta1;

			double a0 = weights[0][i] + weights[1][i];
			double a1 = weights[2][i] + weights[3][i];

			// norm our weights
			weights[0][i] /= a0;
			weights[1][i] /= a0;
			weights[2][i] /= a1;
			weights[3][i] /= a1;

			// compute gradients
			gradientsL[i] = Vector( (value[i] - values[0][i]) / (point.x - points[0].x), 0, 0);
			gradientsR[i] = Vector( (value[i] - values[1][i]) / (point.x - points[1].x), 0, 0);
		};
	};

	//2D case

	//3D case

	return std::move(WENO2PointsStencil());
};


#endif
