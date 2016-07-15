#ifndef TurboStructured_Reconstruction_Linear2psLim
#define TurboStructured_Reconstruction_Linear2psLim

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include "utility/Limiters.h"
#include <Reconstruction/Linear2PointsStencil.h>


// Linear reconstruction for 2 points stencil with a limiter
template<typename limiterType> class Linear2psLim : public Linear2PointsStencil {
public:
	limiterType lim;

	// reconstruction at the point
	virtual std::valarray<double> SampleSolution(Vector const& point) override {
		std::valarray<double> res(nValues);
		for (int i = 0; i < nValues; i++) res[i] = vals[i] + lim[i] * (grads[i] * point);
		return res;
	};

	//constructor
	Linear2psLim() { };
	Linear2psLim(Linear2PointsStencil& base) : Linear2PointsStencil(base) { };
	Linear2psLim(std::valarray<double> _vals, int _nDims, int _nValues, std::vector<Vector> _grads) :
		Linear2PointsStencil(_vals, _nDims, _nValues, _grads) {};

};

template<>
Linear2psLim<limBarsJespersen> ComputeReconstruction< Linear2psLim<limBarsJespersen> >(std::vector<std::valarray<double> >& values, std::vector<CellInfo>& cells, std::valarray<double>& value, CellInfo& cell, int nDim) {
	size_t size = value.size();

	// First compute linear part
	auto linear_rec = ComputeReconstruction<Linear2PointsStencil>(values, cells, value, cell, nDim);

	// Call constructor of the reconstruction
	auto res = Linear2psLim<limBarsJespersen>(linear_rec);

	// Compute reconstructions at all face centers
	std::vector<valarray<double> > projs;

	// 1D case
	projs.push_back(linear_rec.SampleSolution({ -0.5 * cell.hx, 0, 0 }) - linear_rec.vals);
	projs.push_back(linear_rec.SampleSolution({ 0.5 * cell.hx, 0, 0 }) - linear_rec.vals);

	// 2D case
	if (nDim > 1) {
		projs.push_back(linear_rec.SampleSolution({ 0, -0.5 * cell.hy, 0 }) - linear_rec.vals);
		projs.push_back(linear_rec.SampleSolution({ 0, 0.5 * cell.hy, 0 }) - linear_rec.vals);
	};

	// 3D case
	if (nDim > 2) {
		projs.push_back(linear_rec.SampleSolution({ 0, 0, -0.5 * cell.hz }) - linear_rec.vals);
		projs.push_back(linear_rec.SampleSolution({ 0, 0, 0.5 * cell.hz }) - linear_rec.vals);
	};

	// Compute limiter values 
	res.lim.ComputeLimiterValues(values, value, projs, nDim);

	return std::move(res);
};

template<>
Linear2psLim<limVenkatar> ComputeReconstruction< Linear2psLim<limVenkatar> >(std::vector<std::valarray<double> >& values, std::vector<CellInfo>& cells, std::valarray<double>& value, CellInfo& cell, int nDim) {
	size_t size = value.size();

	// First compute linear part
	auto linear_rec = ComputeReconstruction<Linear2PointsStencil>(values, cells, value, cell, nDim);

	// Call constructor of the reconstruction
	auto res = Linear2psLim<limVenkatar>(linear_rec);

	// Compute reconstructions at all face centers
	std::vector<valarray<double> > projs;

	// 1D case
	projs.push_back(linear_rec.SampleSolution({ -0.5 * cell.hx, 0, 0 }) - linear_rec.vals);
	projs.push_back(linear_rec.SampleSolution({ 0.5 * cell.hx, 0, 0 }) - linear_rec.vals);

	// 2D case
	if (nDim > 1) {
		projs.push_back(linear_rec.SampleSolution({ 0, -0.5 * cell.hy, 0 }) - linear_rec.vals);
		projs.push_back(linear_rec.SampleSolution({ 0, 0.5 * cell.hy, 0 }) - linear_rec.vals);
	};

	// 3D case
	if (nDim > 2) {
		projs.push_back(linear_rec.SampleSolution({ 0, 0, -0.5 * cell.hz }) - linear_rec.vals);
		projs.push_back(linear_rec.SampleSolution({ 0, 0, 0.5 * cell.hz }) - linear_rec.vals);
	};

	// Compute limiter values
	res.lim.SendCellInfo(cell, nDim);
	res.lim.ComputeLimiterValues(values, value, projs, nDim);

	return std::move(res);
};

#endif
