#ifndef TurboStructured_RiemannSolvers_RiemannSolver
#define TurboStructured_RiemannSolvers_RiemannSolver

#include <vector>
#include <valarray>
#include "utility\Vector.h"

//Solution information structure
class RiemannProblemSolutionResult {
public:
	std::vector<double> Fluxes;
	double MaxEigenvalue; //Maximal local eigenvalue of Jacobian
	Vector Velocity; //Interface velocity estimate
	double Pressure; //Interface pressure estimate
};

//Base class for all riemann solvers
class RiemannSolver {
public:
	//Get required width of dummy cell layer
	virtual int GetDummyCellLayerSize() = 0;

	//Compute flux given stencil values by reference
	virtual RiemannProblemSolutionResult ComputeFlux(std::vector<std::valarray<double> > &values, Vector faceNormal) = 0;

};

#endif