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
	ENO2PointsStencil(CellReconstruction& _recons, int _nDims, int _nValues) :
	IReconstruction(_recons, _nDims, _nValues) {};

};

template<>
ENO2PointsStencil ComputeReconstruction<ENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim, double gamma) {
	// compute Gradients
	size_t size = value.size();
	std::valarray< Vector > gradients(size);
	double grad_l, grad_r, grad;		// spatial derivatives computed by two neighbour stencils
	CellReconstruction recons;

	//1D case
	recons.xR.resize(size);
	recons.xL.resize(size);
	double deltaL = point.x - points[0].x;
	double deltaR = points[1].x - point.x;
	double delta_x = deltaL * deltaR / (deltaL + deltaR);

	// ENO procedure
	for (auto i = 0; i < size; i++) {
		grad_l = (value[i] - values[0][i]) / deltaL;
		grad_r = (values[1][i] - value[i]) / deltaR;
		if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
		else grad = grad_r;
		recons.xR[i] = value[i] + grad * delta_x;
		recons.xL[i] = value[i] - grad * delta_x;
	};

	//2D case
	if (nDim > 1) {
		recons.yR.resize(size);
		recons.yL.resize(size);
		double deltaL = point.y - points[2].y;
		double deltaR = points[3].y - point.y;
		double delta_y = deltaL * deltaR / (deltaL + deltaR);

		for (auto i = 0; i < size; i++) {
			grad_l = (value[i] - values[2][i]) / deltaL;
			grad_r = (values[3][i] - value[i]) / deltaR;
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			recons.yR[i] = value[i] + grad * delta_y;
			recons.yL[i] = value[i] - grad * delta_y;
		};
	};

	//3D case
	if (nDim > 2) {
		recons.zR.resize(size);
		recons.zL.resize(size);
		double deltaL = point.z - points[4].z;
		double deltaR = points[5].z - point.z;
		double delta_z = deltaL * deltaR / (deltaL + deltaR);
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[4][i]) / deltaL;
			grad_r = (values[5][i] - value[i]) / deltaR;
			if (std::abs(grad_l) < std::abs(grad_r)) grad = grad_l;
			else grad = grad_r;
			recons.zR[i] = value[i] + grad * delta_z;
			recons.zL[i] = value[i] - grad * delta_z;
		};
	};
	
	return std::move(ENO2PointsStencil(recons, nDim, size));
};


#endif
