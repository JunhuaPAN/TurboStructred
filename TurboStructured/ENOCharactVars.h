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



std::vector<std::valarray<double> > GroupTransition(std::vector<std::valarray<double> >& values) {
	// characteristic variables in stencil
	std::vector<std::valarray<double> > CharValues;
	for (auto r : values) {
		std::valarray<double> CharVal(size);
		for (int i = 0; i < size; i++) CharVal[i] = (Rinv[i] * r).sum();
		CharValues.push_back(CharVal);
	};
	return std::move(CharValues);
};
std::valarray<double> SingleTransition(std::valarray<double>& value) {
	// convert values in our cell
	std::valarray<double> CharValue(size);
	for (int i = 0; i < size; i++) CharValue[i] = (Rinv[i] * value).sum();
	return std::move(CharValue);
};
std::valarray<double> InverseTransition(std::valarray<double>& value) {
	// convert values in our cell
	std::valarray<double> CharValue(size);
	for (int i = 0; i < size; i++) CharValue[i] = (R[i] * value).sum();
	return std::move(CharValue);
};
void PrepareTransitionMatrix(std::valarray<double>& UL, std::valarray<double>& UR) {
	ComputeRoeAverageState(UL, UR);
	ComputeLeftEigenVectors();
	ComputeRightEigenVectors();
};


template<>
ENO2CharactVars ComputeReconstruction<ENO2CharactVars>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim, double gamma) {
	// points in stencil
	const int Npoints = 2;

	// Compute Gradients
	size_t size = value.size();
	CellReconstruction recons;
	// RoeLineralisator RL(size);

	// change the coordinates (trivial for X direction)
	std::vector<Vector> newPoints = points;
	Vector newPoint = point;

	// initialize structures for characterictic variables
	std::vector<std::valarray<double> > CharValues;
	std::valarray<double> CharValue;

	// compute left eigenvectors matrix in left face
	ComputeRoeAverageState(values[Npoints - 1], value);
	ComputeLeftEigenVectors();
	ComputeRightEigenVectors();
	CharValues = GroupTransition(values);		// convert all Ui in stencill cells
	CharValue = SingleTransition(value);		// convert values in first cell

	// compute left reconstruction
	recons.xL = InverseTransition( ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, point, 1).recons.xL);

	// compute left eigenvectors matrix in left face
	ComputeRoeAverageState(values[Npoints - 1], value);
	ComputeLeftEigenVectors();
	ComputeRightEigenVectors();
	std::vector<std::valarray<double> > CharValues = GroupTransition(values);		// convert all Ui in stencill cells
	std::valarray<double> CharValue = SingleTransition(value);						// convert values in first cell

	// compute right reconstruction
	recons.xR = InverseTransition(ComputeReconstruction<ENO2PointsStencil>(CharValues, newPoints, CharValue, point, 1).recons.xR);







	return std::move(ENO2CharactVars(recons, nDim, size));
};


#endif

