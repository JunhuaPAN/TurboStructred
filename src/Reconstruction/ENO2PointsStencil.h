#ifndef TurboStructured_Reconstruction_ENO2PointsStencil
#define TurboStructured_Reconstruction_ENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/Linear2PointsStencil.h>


//Reconstruction by ENO interpolation procedure for 2 points stencil
class ENO2PointsStencil : public Linear2PointsStencil {
public:

	//constructor
	ENO2PointsStencil() { };
	ENO2PointsStencil(std::valarray<double> _vals, int _nDims, int _nValues, std::vector<Vector> _grads) :
		Linear2PointsStencil(_vals, _nDims, _nValues, _grads) {};
};

template<>
ENO2PointsStencil ComputeReconstruction<ENO2PointsStencil>(std::vector<std::valarray<double> >& values, std::vector<CellInfo>& cells, std::valarray<double>& value, CellInfo& cell, int nDim) {
	size_t size = value.size();
	std::vector< Vector > gradients(size);
	double grad_l, grad_r, grad;		// spatial derivatives computed by two neighbour stencils

	// ENO procedure for X direction
	for (auto i = 0; i < size; i++) {
		grad_l = (value[i] - values[0][i]) / (cell.x - cells[0].x);
		grad_r = (values[1][i] - value[i]) / (cells[1].x - cell.x);
		if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
		else grad = grad_r;
		gradients[i].x = grad;
	};

	//2D case
	if (nDim > 1) {
		for (auto i = 0; i < size; i++) {
			grad_l = (value[i] - values[2][i]) / (cell.y - cells[2].y);
			grad_r = (values[3][i] - value[i]) / (cells[3].y - cell.y);
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			gradients[i].y = grad;
		};
	};

	//3D case
	if (nDim > 2) {
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[4][i]) / (cell.z - cells[4].z);
			grad_r = (values[5][i] - value[i]) / (cells[5].z - cell.z);
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			gradients[i].z = grad;
		};
	};
	
	return std::move(ENO2PointsStencil(value, nDim, size, gradients));
};


#endif
