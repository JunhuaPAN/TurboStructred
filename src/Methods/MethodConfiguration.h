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

/// Enumerate all available types of limiters
enum class LimiterType {
	BarsJespersen,		///< from Blazek CFD principles and apps 
	Venkatakrishnan		///< from Blazek
};

//Method configuration class
struct MethodConfiguration {

	// Method Settings
	RPSolver RiemannProblemSolver { RPSolver::NoSolver };
	Reconstruction ReconstructionType { Reconstruction::PiecewiseConstant };
	LimiterType GeneralLimitter { LimiterType::BarsJespersen };

	// Parapemeters values
	int RungeKuttaOrder { 0 };
	double CFL { 0.25 };
	double Eps { 0.0 };
};

#endif
