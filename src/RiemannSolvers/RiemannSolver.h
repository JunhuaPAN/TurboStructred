#ifndef TurboStructured_RiemannSolvers_RiemannSolver
#define TurboStructured_RiemannSolvers_RiemannSolver

#include <vector>
#include <valarray>
#include "utility/Vector.h"

// All type of Riemann's Problem solvers
enum class RPSolver {
	GodunovSolver,
	RoePikeSolver,
	NoSolver
};

//Solution information structure
class RiemannProblemSolutionResult {
public:
	std::vector<double> Fluxes;
	double MaxEigenvalue; //Maximal local eigenvalue of Jacobian
	Vector Velocity; //Interface velocity estimate
	double Pressure; //Interface pressure estimate
};

//Base class for all Riemann solvers
class RiemannSolver {
public:

	//	Compute flux given stencil values by reference
	virtual RiemannProblemSolutionResult ComputeFlux(std::valarray<double> &valueL, std::valarray<double> &valueR, Vector faceNormal) = 0;

};

#endif
