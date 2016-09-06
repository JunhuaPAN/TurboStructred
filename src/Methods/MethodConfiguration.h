#ifndef TurboStructured_Methods_MethodConfiguration
#define TurboStructured_Methods_MethodConfiguration

// Enumerate all tupes of reconstruction here
enum class Reconstruction {
	PiecewiseConstant,
	Linear2PointsStencil,
	ENO2PointsStencil,
	Linear2psLim
};

// All type of Riemann's Problem solvers
enum class RPSolver {
	GodunovSolver,
	RoePikeSolver,
	NoSolver
};

//Method configuration class
struct MethodConfiguration {

	// Method Settings
	RPSolver RiemannProblemSolver { RPSolver::NoSolver };
	Reconstruction ReconstructionType { Reconstruction::PiecewiseConstant };

	// Parapemeters values
	int RungeKuttaOrder { 0 };
	double CFL { 0.25 };
	double Eps { 0.0 };
};

#endif
