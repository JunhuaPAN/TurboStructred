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
	// points in stencil
	const int Npoints = 2;

	// Compute Gradients
	size_t size = value.size();
	CellReconstruction recons;
	RoeLineariser RL(size, gamma);

	// hange the coordinates (trivial for X direction)
	std::vector<Vector> newPoints = points;
	Vector newPoint = point;

	// initialize structures for characterictic variables
	std::vector<std::valarray<double> > CharValues;
	std::valarray<double> CharValue;

	// compute left eigenvectors matrix in left face
	RL.PrepareTransitionMatrix(values[Npoints - 2], value);
	CharValues = RL.GroupTransition(values);		// convert all Ui in stencill cells
	CharValue = RL.SingleTransition(value);		// convert values in first cell

	// compute left reconstruction
	recons.xL = RL.InverseTransition( ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, point, 1, gamma).recons.xL);

	// compute left eigenvectors matrix in left face
	RL.PrepareTransitionMatrix(value, values[Npoints - 1]);
	CharValues = RL.GroupTransition(values);		// convert all Ui in stencill cells
	CharValue = RL.SingleTransition(value);		// convert values in first cell

	// compute right reconstruction
	recons.xR = RL.InverseTransition(ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, point, 1, gamma).recons.xR);
	
	// 2D

	// 3D

	return std::move(ENO2CharactVars(recons, nDim, size));
};


#endif

