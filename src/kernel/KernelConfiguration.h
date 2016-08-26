#ifndef TurboStructured_kernel_KernelConfiguration
#define TurboStructured_kernel_KernelConfiguration

#include <map>
#include "Methods/MethodConfiguration.h"
#include "Vector.h"
#include "BoundaryConditions/BoundaryConditionConfiguration.h"

// Grid configuration part
struct BlockNode {
	double pos{ 0.0 };				// Coordinate of node
	double q_com{ 1.0 };	// Compression coefficient (q>1); Expansion (q<1)
	int N_cells{ -1 };			// Number of cells in block

	// constructors
	BlockNode(double _pos) : pos(_pos) {};
	BlockNode() {};
};

//Class that manages all configurable parameters
class KernelConfiguration {
public:
	//Grid configuration
	int nDims{ 1 }; //Number of dimensions
	int nX{ 1 }; //Number of cells in x dimension
	int nY{ 1 }; //Number of cells in y dimension
	int nZ{ 1 }; //Number of cells in z dimension
	bool isUniformAlongY{ true };
	bool isUniformAlongZ{ true };
	bool isPeriodicX{ true }; //X periodicity
	bool isPeriodicY{ true }; //Y periodicity
	bool isPeriodicZ{ true }; //Z periodicity
	double LX{ 1.0 }; //X size
	double LY{ 1.0 }; //Y size
	double LZ{ 1.0 }; //Z size
	// Set uniform block distribution
	std::vector< BlockNode > CompressionX{ std::vector< BlockNode >(1) };
	std::vector< BlockNode > CompressionY{ std::vector< BlockNode >(1) };
	std::vector< BlockNode > CompressionZ{ std::vector< BlockNode >(1) };

	//Gas model parameters
	double Gamma{ 1.4 };
	double Viscosity{ 0.0 };
	double ThermalConductivity{ 0.0 };
	double MolarMass{ 0.0 };
	double UniversalGasConstant{ 8.3144598 };

	//model configuration
	bool IsViscousFlow{ false };
	bool IsExternalForceRequared{ false };
	bool IsUnifromAccelerationRequared{ false };

	//Solution method parameters
	int DummyLayerSize{ 1 };
	enum class Method {
		ExplicitRungeKuttaFVM,
	} SolutionMethod { Method::ExplicitRungeKuttaFVM  };

	MethodConfiguration methodConfiguration;

	//Run parameters
	int MaxIteration{ 1000 };
	double MaxTime{ 1.0 };
	double SaveSolutionTime{ 0 };
	double SaveSliceTime{ 0 };
	int SaveSolutionIterations{ 0 };
	int SaveSliceIterations{ 0 };
	int ResidualOutputIterations{ 0 };
	bool DebugOutputEnabled{ false };
	
	// Boundary conditions configuration
	std::map<int, BoundaryConditionConfiguration> MyConditions;

	// By default all sides have just one BC
	DomainBCinfo xLeftBoundary;		
	DomainBCinfo xRightBoundary;
	DomainBCinfo yLeftBoundary;
	DomainBCinfo yRightBoundary;
	DomainBCinfo zLeftBoundary;
	DomainBCinfo zRightBoundary;

	//External potential forces
	Vector Sigma{ Vector(0,0,0) };			            //!< Pressure gradient
	Vector UniformAcceleration{ Vector(0,0,0) };		//!< Uniform external gravity field

	//Computation continuation regime
	std::string fileToLoad{ "" };
	bool ContinueComputation{ false };	
};

#endif
