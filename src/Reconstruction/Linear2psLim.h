#ifndef TurboStructured_Reconstruction_Linear2psLim
#define TurboStructured_Reconstruction_Linear2psLim

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/Linear2PointsStencil.h>


// Linear reconstruction for 2 points stencil with a limiter
template<typename LimiterType> class Linear2psLim : public Linear2PointsStencil {
private:
	LimiterType lim;
public:

	// reconstruction at the point
	virtual std::valarray<double> SampleSolution(Vector const& point) override {
		std::valarray<double> res(nValues);
		for (int i = 0; i < nValues; i++) res[i] = vals[i] + lim * (grads[i] * point);
		return res;
	};

	//constructor
	Linear2psLim() { };
	Linear2psLim(std::valarray<double> _vals, int _nDims, int _nValues, std::vector<Vector> _grads) :
		Linear2PointsStencil(_vals, _nDims, _nValues, _grads) {};

	// Set limiter
	void SetLimiter(LimiterType _lim) {
		lim = _lim;
	};
};

template<>
Linear2psLim<double> ComputeReconstruction< Linear2psLim<double> >(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	size_t size = value.size();
	std::vector< Vector > gradients(size);
	double grad_l, grad_r, grad;		// spatial derivatives computed by two neighbour stencils

	// ENO procedure for X direction
	for (auto i = 0; i < size; i++) {
		grad_l = (value[i] - values[0][i]) / (point.x - points[0].x);
		grad_r = (values[1][i] - value[i]) / (points[1].x - point.x);
		if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
		else grad = grad_r;
		gradients[i].x = grad;
	};

	// 2D case
	if (nDim > 1) {
		for (auto i = 0; i < size; i++) {
			grad_l = (value[i] - values[2][i]) / (point.y - points[2].y);
			grad_r = (values[3][i] - value[i]) / (points[3].y - point.y);
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			gradients[i].y = grad;
		};
	};

	// 3D case
	if (nDim > 2) {
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[4][i]) / (point.z - points[4].z);
			grad_r = (values[5][i] - value[i]) / (points[5].z - point.z);
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			gradients[i].z = grad;
		};
	};

	// Create reconstruction with limiter
	auto res = Linear2psLim<double>(value, nDim, size, gradients);
	double lim = 0;	// just first order
	res.SetLimiter(lim);

	return std::move(res);
};


#endif
