#ifndef TurboStructured_Reconstruction_ENO2PointsStencil
#define TurboStructured_Reconstruction_ENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/IReconstruction.h>


//Reconstruction by ENO interpolation procedure for 2 points stencil
class ENO2PointsStencil : public IReconstruction {
public:

	//constructor
	ENO2PointsStencil() { };
	ENO2PointsStencil(std::valarray<double>& _xR, std::valarray<double>& _xL, std::valarray<double>& _yR, std::valarray<double>& _yL, std::valarray<double>& _zR, std::valarray<double>& _zL, int _nDims, int _nValues) :
	IReconstruction(_xR, _xL, _yR, _yL, _zR, _zL, _nDims, _nValues) {};

};

template<>
ENO2PointsStencil ComputeReconstruction<ENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	//compute Gradients
  size_t size = value.size();
	std::valarray< Vector > gradients(size);
	double grad_l, grad_r, grad;		// spatial derivatives computed by two neighbour stencils

	//1D case
	if (nDim == 1) {
		for (auto i = 0; i < size; i++) {
			grad_l = (value[i] - values[0][i]) / (point.x - points[0].x);
			grad_r = (values[1][i] - value[i]) / (points[1].x - point.x);
			if (std::abs(grad_l) < std::abs(grad_r)) gradients[i] = Vector(grad_l, 0, 0);
			else gradients[i] = Vector(grad_r, 0, 0);
		};
	};

	//2D case
	if (nDim == 2) {
		for (auto i = 0; i < size; i++) {
			grad_l = (value[i] - values[2][i]) / (point.y - points[2].y);
			grad_r = (values[3][i] - value[i]) / (points[3].y - point.y);
			if (std::abs(grad_l) < std::abs(grad_r)) gradients[i] = gradients[i] + Vector(0, grad_l, 0);
			else gradients[i] = gradients[i] + Vector(0, grad_r, 0);
		};
	};

	//3D case
	if (nDim == 3) {
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[4][i]) / (point.z - points[4].z);
			grad_r = (values[5][i] - value[i]) / (points[5].z - point.z);
			if (std::abs(grad_l) < std::abs(grad_r)) gradients[i] = gradients[i] + Vector(0, 0, grad_l);
			else gradients[i] = gradients[i] + Vector(0, 0, grad_r);
		};
	};
	
	//ENO2PointsStencil res = ENO2PointsStencil(value, point, gradients);
	return std::move(ENO2PointsStencil());
};


#endif
