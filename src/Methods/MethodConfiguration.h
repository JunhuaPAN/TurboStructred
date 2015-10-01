#ifndef TurboStructured_Methods_MethodConfiguration
#define TurboStructured_Methods_MethodConfiguration

#include "RiemannSolvers/RiemannSolver.h"
#include "Reconstruction/IReconstruction.h"

//Method configuration class
struct MethodConfiguration {

	// Method Settings
	RPSolver RiemannProblemSolver { RPSolver::NoSolver };
	Reconstruction ReconstructionType { Reconstruction::PiecewiseConstant };

	// Parapemeters values
	int RungeKuttaOrder { 0 };
	double CFL { 0.25 };
	double Eps { 0.0 };
	double OperatingPresure { 0 };
};

#endif
