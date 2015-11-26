#ifndef TurboStructured_Reconstruction_ENOCharactVars
#define TurboStructured_Reconstruction_ENOCharactVars

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/IReconstruction.h>
#include <Reconstruction/ENO2PointsStencil.h>


//Reconstruction by ENO interpolation procedure for 2 points stencil
class ENO2CharactVars : public IReconstruction {
public:

	//constructor
	ENO2CharactVars() { };
	ENO2CharactVars(CellReconstruction& _recons, int _nDims, int _nValues) :
		IReconstruction(_recons, _nDims, _nValues) {};

};


template<>
ENO2CharactVars ComputeReconstruction<ENO2CharactVars>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim, double gamma) {

	// Compute Gradients
	size_t size = value.size();
	CellReconstruction recons;
	RoeLineariser RL(size, gamma);

	// initialize the coordinates vectors (trivial for X direction)
	std::vector<Vector> newPoints{ points[0], points[1] };
	Vector newPoint = point;

	// initialize subset of values vector
	std::vector< std::valarray< double > > newValues{ values[0], values[1] };

	// lambda implies ENO 1D procedure for both faces in the direction dir
	// return first for left reconstruction and second for right one
	auto Reconst1D = [&](Direction dir) -> std::pair<std::valarray<double>, std::valarray<double> > {
		std::pair<std::valarray<double>, std::valarray<double> > res;

		// compute eigenvectors matrix in left face
		RL.SetDirection(dir);
		RL.PrepareTransitionMatrix(newValues[0], value);
		auto CharValues = RL.GroupTransition(values);		// convert all Ui in stencill cells
		auto CharValue = RL.SingleTransition(value);		// convert values in first cell

		// compute left reconstruction
		res.first = RL.InverseTransition(ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, newPoint, 1, gamma).recons.xL);

		// compute eigenvectors matrix in left face
		RL.PrepareTransitionMatrix(value, newValues[1]);
		CharValues = RL.GroupTransition(values);			// convert all Ui in stencill cells
		CharValue = RL.SingleTransition(value);				// convert values in first cell

		// compute right reconstruction
		res.second = RL.InverseTransition(ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, newPoint, 1, gamma).recons.xR);

		return res;
	};
	
	// 1D case
	auto rec = Reconst1D(Direction::XDirection);
	recons.xL = rec.first;
	recons.xR = rec.second;
		
	// 2D
	if (nDim > 1) {
		// update coordinates
		newPoints[0] = { points[2].y, 0, 0 };
		newPoints[1] = { points[3].y, 0, 0 };
		newPoint = { point.y, 0 , 0};

		// update values
		newValues = { values[2], values[3] };

		// compute reconstructions
		rec = Reconst1D(Direction::YDirection);
		recons.yL = rec.first;
		recons.yR = rec.second;
	};

	// 3D
	if (nDim > 2) {
		// update coordinates
		newPoints[0] = { points[2].z, 0, 0 };
		newPoints[1] = { points[3].z, 0, 0 };
		newPoint = { point.z, 0 , 0 };

		// update values
		newValues = { values[4], values[5] };

		// compute reconstructions
		rec = Reconst1D(Direction::ZDirection);
		recons.zL = rec.first;
		recons.zR = rec.second;
	};

	return std::move(ENO2CharactVars(recons, nDim, size));
};


#endif

