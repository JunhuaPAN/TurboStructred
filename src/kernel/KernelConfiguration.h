#ifndef TurboStructured_kernel_KernelConfiguration
#define TurboStructured_kernel_KernelConfiguration

#include "Methods\MethodConfiguration.h"
#include "Vector.h"

enum class BoundaryConditionType {
	Wall,
	SymmetryX,
	SymmetryY,
	SymmetryZ,
	MovingWall,
	General,
	Natural
};

class BoundaryConditionConfiguration {
public:
	BoundaryConditionType BCType;
	Vector Velocity;
	double Gamma;
};

//Class that manages all configurable parameters
class KernelConfiguration {
public:
	//Grid sizes and configuration
	int nDims; //Number of dimensions
	int nX; //Number of cells in x dimension
	int nY; //Number of cells in y dimension
	int nZ; //Number of cells in z dimension	
	bool isPeriodicX; //X periodicity
	bool isPeriodicY; //Y periodicity
	bool isPeriodicZ; //Z periodicity
	double LX; //X size
	double LY; //Y size
	double LZ; //Z size

	//Gas model parameters
	double Gamma;
	double Viscosity;
	double ThermalConductivity;

	//model configuration
	bool IsViscousFlow;
	bool IsExternalForceRequared;
	bool IsUnifromAccelerationRequared;

	//Solution method parameters
	enum class Method {
		HybridFVM,
		ExplicitRungeKuttaFVM,
		HybridGeneralEOSOnePhase,
		HybridBarotropicEOSOnePhase,
		HybridBarotropicEOSTwoPhase
	} SolutionMethod;

	MethodConfiguration methodConfiguration;

	//Run parameters
	double MaxTime;
	double MaxIteration;
	double SaveSolutionSnapshotTime;	
	int SaveSolutionSnapshotIterations;
	int ResidualOutputIterations;
	bool DebugOutputEnabled;
	
	//Boundary conditions configuration
	BoundaryConditionConfiguration xLeftBoundary;
	BoundaryConditionConfiguration xRightBoundary;
	BoundaryConditionConfiguration yLeftBoundary;
	BoundaryConditionConfiguration yRightBoundary;
	BoundaryConditionConfiguration zLeftBoundary;
	BoundaryConditionConfiguration zRightBoundary;

	//External potential forces
	Vector Sigma;			//dP/dn
	Vector UniformAcceleration;		// g acceleration

	//Computation continuation regime
	std::string fileToLoad;
	bool ContinueComputation;

	//Default values
	KernelConfiguration() {
		Viscosity = 0;
		ThermalConductivity = 0;

		IsViscousFlow = false;
		IsExternalForceRequared = false;
		IsUnifromAccelerationRequared = false;
		DebugOutputEnabled = false;
		ContinueComputation = false;
	};
};

#endif