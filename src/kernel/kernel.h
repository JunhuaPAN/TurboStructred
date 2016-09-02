#ifndef TurboStructured_kernel_kernel
#define TurboStructured_kernel_kernel

#include <array>
#include <vector>
#include <chrono>
#include <thread>
#include <functional>
#include <sstream>
#include <cassert>
#include <fstream>
#include <valarray>
#include <memory>

#include "KernelConfiguration.h"
#include "grid.h"
#include "GasProperties.h"
#include "ParallelManager.h"
#include "utility/Vector.h"
#include "utility/Timer.h"
#include "utility/Slices.h"
#include "utility/NumericQuadrature.h"
#include "utility/Stencil.h"
#include "utility/VariablesTransition.h"
#include "RiemannSolvers/RiemannSolversList.h"
#include "BoundaryConditions/BCGeneral.h"
#include "Sensors/Sensors.h"


// Step info
class StepInfo {
public:
	double Time;
	double TimeStep;	
	int Iteration;
	std::vector<double> Residual;
	double NextSolutionSnapshotTime;
	double NextSliceSnapshotTime;
};


// Calculation kernel
class Kernel {
public:
	virtual void IterationStep() = 0;		// main function
	virtual void InitializeMethod(KernelConfiguration& config) = 0;	// initialization of parameters of the Method

	// Parallel run information
	std::unique_ptr<ParallelManager> pManager;
	int rank; //Rank of process

	// Structured grid and dimensions number
	int nDims;
	Grid grid;

	// Current step information
	StepInfo stepInfo;

	// Gas model information
	int nVariables;				// number of conservative variables
	GasProperties gas_prop;		// all parameters of the gas
	ValuesTransition compute;	// functions for computation of primituve values

	bool isViscousFlow;
	bool isGradientRequired;
	bool isExternalAccelaration;
	bool isExternalForce;

	// Sensors operating
	bool isSensorEnable{ false };
	int SaveSensorRecordIterations{ 1 };

	// External forces
	Vector Sigma; //Potential force like presure gradient
	Vector UniformAcceleration;	//external uniform acceleration

	// Boundary conditions TO DO refactor
	std::map<int, std::unique_ptr<BoundaryConditions::BCGeneral> > bConditions;

	BCinfo xLeftBC;
	BCinfo xRightBC;
	BCinfo yLeftBC;
	BCinfo yLeftBCspecial;
	BCinfo yRightBC;
	BCinfo zLeftBC;
	BCinfo zRightBC;

	// Solution data
	std::valarray<double> values;
	std::valarray<double> residual;
	std::vector<Slice> slices;

	// Get cell values
	inline double* getCellValues(int i, int j, int k) {
		int sBegin = grid.getSerialIndexLocal(i, j, k) * nVariables;
		return &values[sBegin];
	};

	// Calculation parameters	
	int MaxIteration;
	double MaxTime;
	double SaveSolutionTime;
	double SaveSliceTime;
	double SaveBinarySolTime;
	int SaveSolutionIterations;
	int SaveSliceIterations;
	int SaveBinarySolIterations;

	int ResidualOutputIterations{ };
	bool ContinueComputation{ false };
	bool DebugOutputEnabled{ false };

	//Array of sensors in use
	std::vector<std::unique_ptr<Sensor>> Sensors;

	//Constructor
	Kernel(int* argc, char **argv[]) : pManager(new ParallelManager(argc, argv)) {
		nVariables = 5;	//default value
		rank = pManager->getRank();
	};

	// Set initial conditions using exact values in cell centers
	void SetInitialConditions(std::function<std::vector<double>(const Vector& r)> initF) {
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					//Obtain cell data
					double x = grid.CoordinateX[i];
					double y = grid.CoordinateY[j];
					double z = grid.CoordinateZ[k];
					int sBegin = grid.getSerialIndexLocal(i, j, k) * nVariables;
					std::vector<double> U = initF(Vector(x, y, z));
					for (int varn = 0; varn < nVariables; varn++) values[sBegin + varn] = U[varn];
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Initial conditions written.\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	// Set initial conditions using cell averaged values
	void SetInitialConditions(std::function<std::vector<double>(const Vector& r)> initF, NumericQuadrature& Int) {
		CellInfo cell;
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					//Obtain cell data
					cell.x = grid.CoordinateX[i];
					cell.y = grid.CoordinateY[j];
					cell.z = grid.CoordinateZ[k];
					cell.hx = grid.hx[i];
					cell.hy = grid.hy[j];
					cell.hz = grid.hz[k];
					double cellV = cell.hx * cell.hy * cell.hz;

					int sBegin = grid.getSerialIndexLocal(i, j, k) * nVariables;
					values[slice(sBegin, nVariables, 1)] = Int.IntegrateGroap(cell, initF, nVariables) / cellV;
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Initial conditions written.\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	// Initialize kernel  ( TO DO READ THE CONFIG FILE)
	virtual void Init(KernelConfiguration& config) {
		// Initialize MPI		
		nDims = config.nDims;

		// Output settings 
		DebugOutputEnabled = config.DebugOutputEnabled;

		// Initialize global grid
		grid.InitGlobal(config);

		// Initialize cartezian topology
		pManager->InitCartesianTopology(grid);

		// Compute local grid sizes and initialize local grid
		grid.nlocalX = grid.nX / pManager->dimsCart[0];
		grid.nlocalY = grid.nY / pManager->dimsCart[1];
		grid.nlocalZ = grid.nZ / pManager->dimsCart[2];
		grid.iMin = pManager->rankCart[0] * grid.nlocalX + grid.dummyCellLayersX;
		grid.iMax = (pManager->rankCart[0] + 1) * grid.nlocalX + grid.dummyCellLayersX - 1;
		grid.jMin = pManager->rankCart[1] * grid.nlocalY + grid.dummyCellLayersY;
		grid.jMax = (pManager->rankCart[1] + 1) * grid.nlocalY + grid.dummyCellLayersY - 1;
		grid.kMin = pManager->rankCart[2] * grid.nlocalZ + grid.dummyCellLayersZ;
		grid.kMax = (pManager->rankCart[2] + 1) * grid.nlocalZ + grid.dummyCellLayersZ - 1;
		grid.InitLocal(config);		// end of grid initialization

		//Sync
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->_rankCart << ", iMin = " << grid.iMin << ", iMax = " << grid.iMax << "\n";
			std::cout << "rank = " << pManager->_rankCart << ", jMin = " << grid.jMin << ", jMax = " << grid.jMax << "\n";
			std::cout << "rank = " << pManager->_rankCart << ", kMin = " << grid.kMin << ", kMax = " << grid.kMax << "\n";
		};
		pManager->Barrier();

		//Initialize gas model parameters and Riemann solver
		gas_prop.gamma = config.Gamma;
		gas_prop.thermalConductivity = config.ThermalConductivity;
		gas_prop.molarMass = config.MolarMass;
		gas_prop.universalGasConstant = config.UniversalGasConstant;
		gas_prop.viscosity = config.Viscosity;
		compute.BindGasProperties(gas_prop);

		// Set all flags as configuration requires
		if (config.IsViscousFlow == true) {
			isViscousFlow = true;
			isGradientRequired = true;
		} else {
			isViscousFlow = false;
			isGradientRequired = false;
		};
		if (config.IsExternalForceRequared == true) {
			isExternalForce = true;
		} else {
			isExternalForce = false;
		};
		if (config.IsUnifromAccelerationRequared == true) {
			isExternalAccelaration = true;
		} else {
			isExternalAccelaration = false;
		};
		ContinueComputation = config.ContinueComputation;

		//Allocate data structures
		values.resize(nVariables * grid.nlocalXAll * grid.nlocalYAll * grid.nlocalZAll);
		residual.resize(nVariables * grid.nlocalXAll * grid.nlocalYAll * grid.nlocalZAll);

		//Initialize calculation parameters
		MaxTime = config.MaxTime;
		MaxIteration = config.MaxIteration;
		SaveSolutionTime = config.SaveSolutionTime;
		SaveSliceTime = config.SaveSliceTime;
		SaveSolutionIterations = config.SaveSolutionIterations;
		SaveSliceIterations = config.SaveSliceIterations;
		ResidualOutputIterations = config.ResidualOutputIterations;

		// Initialize step information
		if (ContinueComputation == false) {
			stepInfo.Time = 0;
			stepInfo.Iteration = 0;
			stepInfo.NextSolutionSnapshotTime = SaveSolutionTime;
			stepInfo.NextSliceSnapshotTime = SaveSliceTime;
		};

		// Initialize boundary conditions
		InitBoundaryConditions(config);

		// External forces
		Sigma = config.Sigma;
		UniformAcceleration = config.UniformAcceleration;

		// Initialize method by method configuration
		InitializeMethod(config);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	// Initialize boundary conditions
	virtual void InitBoundaryConditions(KernelConfiguration& config) {
		for (auto r : config.MyConditions) {
			bConditions[r.first] = BoundaryConditions::CreateBC(r.second);
			bConditions[r.first]->BindGasProperties(gas_prop);
			bConditions[r.first]->loadConfiguration(r.second);
		};

		// Check if periodic condition is used then initialize 
		if (!grid.IsPeriodicX) {
			if (config.xLeftBoundary.isComplex) xLeftBC = BCinfo(config.xLeftBoundary.getMarker);		// for complex case transmit the function getMarker
			else xLeftBC = BCinfo(config.xLeftBoundary.GetMarker());		// for trivial case set marker value

			// Do the same for right side
			if (config.xRightBoundary.isComplex) xRightBC = BCinfo(config.xRightBoundary.getMarker);
			else xRightBC = BCinfo(config.xRightBoundary.GetMarker());
		};

		// Repeat process for other directions
		if ((!grid.IsPeriodicY) && (nDims > 1)) {
			if (config.yLeftBoundary.isComplex) yLeftBC = BCinfo(config.yLeftBoundary.getMarker);
			else yLeftBC = BCinfo(config.yLeftBoundary.GetMarker());

			if (config.yRightBoundary.isComplex) yRightBC = BCinfo(config.yRightBoundary.getMarker);
			else yRightBC = BCinfo(config.yRightBoundary.GetMarker());
		};
		if ((!grid.IsPeriodicZ) && (nDims > 2)) {
			if (config.zLeftBoundary.isComplex) zLeftBC = BCinfo(config.zLeftBoundary.getMarker);
			else zLeftBC = BCinfo(config.zLeftBoundary.GetMarker());

			if (config.zRightBoundary.isComplex) zRightBC = BCinfo(config.zRightBoundary.getMarker);
			else zRightBC = BCinfo(config.zRightBoundary.GetMarker());
		};
	};

	// Update solution
	void UpdateSolution(double dt) {
		stepInfo.Residual.resize(nVariables, 0);		

		for (auto nv = 0; nv < nVariables; nv++) {
			std::slice val_slice(nv, grid.nCellsLocalAll, nVariables);
			std::valarray<double> c = dt * ((std::valarray<double>)residual[val_slice] / grid.volumes);
			values[val_slice] += c;
			stepInfo.Residual[nv] = ( grid.volumes * abs((std::valarray<double>)residual[val_slice]) ).sum();
		};

		//Aggregate
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = pManager->Sum(stepInfo.Residual[nv]);
	};

	// Exchange values between processors
	void ExchangeValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//Get parallel run info
		MPI_Comm comm = pManager->getComm();

		//Buffers
		std::vector<double> bufferToRecv;
		std::vector<double> bufferToSend;
		int rankLeft{ -1 }; // if equals -1 no left neighbour
		int rankRight{ -1 }; // if equals -1 no right neighbour

		//X direction exchange		

		//Allocate buffers
		int layerSize = grid.nlocalY * grid.nlocalZ;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);

		//Determine neighbours' ranks
		int rankL = pManager->GetRankByCartesianIndexShift(-1, 0, 0);
		int rankR = pManager->GetRankByCartesianIndexShift(+1, 0, 0);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		for (int layer = 1; layer <= grid.dummyCellLayersX; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int iSend = grid.iMin + layer - 1; // layer index to send
			int iRecv = grid.iMax + layer; // layer index to recv

									  // Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (j = grid.jMin; j <= grid.jMax; j++) {
					for (k = grid.kMin; k <= grid.kMax; k++) {
						int idxBuffer = (j - grid.jMin) + (k - grid.kMin)* grid.nlocalY; //Exclude x index
						double* U = getCellValues(iSend, j, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = U[nv];
					};
				};
			};

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; iSend = " << iSend <<
					"; iRecv = " << iRecv <<
					std::endl << std::flush;
			};

			//Determine recive number
			if (rankR != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankL, rankR, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (j = grid.jMin; j <= grid.jMax; j++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					int idxBuffer = (j - grid.jMin) + (k - grid.kMin)* grid.nlocalY; //Exclude x index
					double* U = getCellValues(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) U[nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			iSend = grid.iMax - layer + 1; // layer index to send
			iRecv = grid.iMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; iSend = " << iSend <<
					"; iRecv = " << iRecv <<
					std::endl << std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (j = grid.jMin; j <= grid.jMax; j++) {
					for (k = grid.kMin; k <= grid.kMax; k++) {
						int idxBuffer = (j - grid.jMin) + (k - grid.kMin)* grid.nlocalY; //Exclude x index
						double *U = getCellValues(iSend, j, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = U[nv];
					};
				};
			};

			//Determine recive number
			if (rankL != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankR, rankL, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (j = grid.jMin; j <= grid.jMax; j++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					int idxBuffer = (j - grid.jMin) + (k - grid.kMin)* grid.nlocalY; //Exclude x index
					double *U = getCellValues(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) U[nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange values X-direction executed\n";
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction exchange		

		//Allocate buffers
		layerSize = grid.nlocalX * grid.nlocalZ;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		for (int layer = 1; layer <= grid.dummyCellLayersY; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int jSend = grid.jMin + layer - 1; // layer index to send
			int jRecv = grid.jMax + layer; // layer index to recv

			// Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = grid.iMin; i <= grid.iMax; i++) {
					for (k = grid.kMin; k <= grid.kMax; k++) {
						int idxBuffer = (i - grid.iMin) + (k - grid.kMin) * grid.nlocalX; //Exclude y index
						int idxValues = grid.getSerialIndexLocal(i, jSend, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = values[idxValues * nVariables + nv];
					};
				};
			};

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; jSend = " << jSend <<
					"; jRecv = " << jRecv <<
					std::endl << std::flush;
			};

			//Determine recive number
			if (rankR != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankL, rankR, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					int idxBuffer = (i - grid.iMin) + (k - grid.kMin) * grid.nlocalX; //Exclude y index
					int idxValues = grid.getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			jSend = grid.jMax - layer + 1; // layer index to send
			jRecv = grid.jMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; jSend = " << jSend <<
					"; jRecv = " << jRecv <<
					std::endl << std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = grid.iMin; i <= grid.iMax; i++) {
					for (k = grid.kMin; k <= grid.kMax; k++) {
						int idxBuffer = (i - grid.iMin) + (k - grid.kMin) * grid.nlocalX; //Exclude y index
						int idxValues = grid.getSerialIndexLocal(i, jSend, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = values[idxValues * nVariables + nv];
					};
				};
			};

			//Determine recive number
			if (rankL != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankR, rankL, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					int idxBuffer = (i - grid.iMin) + (k - grid.kMin) * grid.nlocalX; //Exclude y index
					int idxValues = grid.getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange values Y-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction exchange		

		//Allocate buffers
		layerSize = grid.nlocalX * grid.nlocalY;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, 0, -1);
		rankR = pManager->GetRankByCartesianIndexShift(0, 0, +1);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};


		for (int layer = 1; layer <= grid.dummyCellLayersZ; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int kSend = grid.kMin + layer - 1; // layer index to send
			int kRecv = grid.kMax + layer; // layer index to recv

									  // Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = grid.iMin; i <= grid.iMax; i++) {
					for (j = grid.jMin; j <= grid.jMax; j++) {
						int idxBuffer = (i - grid.iMin) + (j - grid.jMin) * grid.nlocalX; //Exclude z index
						int idxValues = grid.getSerialIndexLocal(i, j, kSend);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = values[idxValues * nVariables + nv];
					};
				};
			};

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; kSend = " << kSend <<
					"; kRecv = " << kRecv <<
					std::endl << std::flush;
			};

			//Determine recive number
			if (rankR != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankL, rankR, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (j = grid.jMin; j <= grid.jMax; j++) {
					int idxBuffer = (i - grid.iMin) + (j - grid.jMin) * grid.nlocalX; //Exclude z index
					int idxValues = grid.getSerialIndexLocal(i, j, kRecv);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			kSend = grid.kMax - layer + 1; // layer index to send
			kRecv = grid.kMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; kSend = " << kSend <<
					"; kRecv = " << kRecv <<
					std::endl << std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = grid.iMin; i <= grid.iMax; i++) {
					for (j = grid.jMin; j <= grid.jMax; j++) {
						int idxBuffer = (i - grid.iMin) + (j - grid.jMin) * grid.nlocalX; //Exclude z index
						int idxValues = grid.getSerialIndexLocal(i, j, kSend);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = values[idxValues * nVariables + nv];
					};
				};
			};

			//Determine recive number
			if (rankL != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankR, rankL, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (j = grid.jMin; j <= grid.jMax; j++) {
					int idxBuffer = (i - grid.iMin) + (j - grid.jMin) * grid.nlocalX; //Exclude z index
					int idxValues = grid.getSerialIndexLocal(i, j, kRecv);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange values Z-direction executed\n";
			std::cout.flush();
		};

		//Sync
		pManager->Barrier();

		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank << ", Exchange values finished." <<
				std::endl << std::flush;
		};

	}; // function

	// Compute dummy cell values as result of boundary conditions and interprocessor exchange communication
	void ComputeDummyCellValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank << ", Dummy cell processing started." <<
				std::endl << std::flush;
		};

		//Current face and cell information
		Vector faceNormalL;
		Vector faceNormalR;
		Vector faceCenter;
		Vector cellCenter;

		//X direction		
		faceNormalL = Vector(-1.0, 0.0, 0.0);
		faceNormalR = Vector(1.0, 0.0, 0.0);
		if (!grid.IsPeriodicX) {
			for (j = grid.jMin; j <= grid.jMax; j++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					//Inner cell
					cellCenter.y = grid.CoordinateY[j];
					cellCenter.z = grid.CoordinateZ[k];

					//And face
					faceCenter.y = grid.CoordinateY[j];
					faceCenter.z = grid.CoordinateZ[k];

					if (DebugOutputEnabled) {
						std::cout << "rank = " << rank <<
							"; f.y = " << faceCenter.y <<
							"; f.z = " << faceCenter.z <<
							std::endl << std::flush;
					};

					for (int layer = 1; layer <= grid.dummyCellLayersX; layer++) {
						if (pManager->rankCart[0] == 0) {
							//Left border
							i = grid.iMin - layer; // layer index
							int iIn = grid.iMin + layer - 1; // opposite index
							cellCenter.x = grid.CoordinateX[iIn];
							faceCenter.x = (grid.CoordinateX[iIn] + grid.CoordinateX[i]) / 2.0;		// TO DO WRONG for nonuniform grid
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(iIn, j, k);

							//Apply left boundary conditions
							auto bcMarker = xLeftBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[0] == pManager->dimsCart[0] - 1) {
							//Right border
							i = grid.iMax + layer; // layer index
							int iIn = grid.iMax - layer + 1; // opposite index
							cellCenter.x = grid.CoordinateX[iIn];
							faceCenter.x = (grid.CoordinateX[iIn] + grid.CoordinateX[i]) / 2.0;		// TO DO WRONG for nonuniform grid
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(iIn, j, k);

							//Apply right boundary conditions						
							auto bcMarker = xRightBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};
					};
				};
			};
		}; //X direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Dummy values X-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction		
		faceNormalL = Vector(0.0, -1.0, 0.0);
		faceNormalR = Vector(0.0, 1.0, 0.0);
		if (!grid.IsPeriodicY) {
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (k = grid.kMin; k <= grid.kMax; k++) {
					//Inner cell
					cellCenter.x = grid.CoordinateX[i];
					cellCenter.z = grid.CoordinateZ[k];

					//And face
					faceCenter.x = grid.CoordinateX[i];
					faceCenter.z = grid.CoordinateZ[k];

					//Debug info
					if (DebugOutputEnabled) {
						std::cout << "rank = " << rank <<
							"; f.x = " << faceCenter.x <<
							"; f.z = " << faceCenter.z <<
							std::endl << std::flush;
					};

					for (int layer = 1; layer <= grid.dummyCellLayersY; layer++) {
						if (pManager->rankCart[1] == 0) {
							//Left border
							j = grid.jMin - layer; // layer index
							int jIn = grid.jMin + layer - 1; // opposite index
							cellCenter.y = grid.CoordinateY[jIn];
							faceCenter.y = (grid.CoordinateY[jIn] + grid.CoordinateY[j]) / 2.0;	//for uniform dummy cells		TO DO MODIFY
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(i, jIn, k);

							//Apply left boundary conditions						
							auto bcMarker = yLeftBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[1] == pManager->dimsCart[1] - 1) {
							//Right border
							j = grid.jMax + layer; // layer index
							int jIn = grid.jMax - layer + 1; // opposite index
							cellCenter.y = grid.CoordinateY[jIn];
							faceCenter.y = (grid.CoordinateY[jIn] + grid.CoordinateY[j]) / 2.0;
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(i, jIn, k);

							//Apply right boundary conditions						
							auto bcMarker = yRightBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};
					};
				};
			};
		}; //Y direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Dummy values Y-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction		
		faceNormalL = Vector(0.0, 0.0, -1.0);
		faceNormalR = Vector(0.0, 0.0, 1.0);
		if (!grid.IsPeriodicZ) {
			for (i = grid.iMin; i <= grid.iMax; i++) {
				for (j = grid.jMin; j <= grid.jMax; j++) {
					//Inner cell
					cellCenter.x = grid.CoordinateX[i];
					cellCenter.y = grid.CoordinateY[j];

					//And face
					faceCenter.x = grid.CoordinateX[i];
					faceCenter.y = grid.CoordinateY[j];

					//Debug info
					if (DebugOutputEnabled) {
						std::cout << "rank = " << rank <<
							"; f.x = " << faceCenter.x <<
							"; f.y = " << faceCenter.y <<
							std::endl << std::flush;
					};

					for (int layer = 1; layer <= grid.dummyCellLayersY; layer++) {

						if (pManager->rankCart[2] == 0) {
							//Left border
							k = grid.kMin - layer; // layer index
							int kIn = grid.kMin + layer - 1; // opposite index
							cellCenter.z = grid.CoordinateZ[kIn];
							faceCenter.z = (grid.CoordinateZ[kIn] + grid.CoordinateZ[k]) / 2.0;
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(i, j, kIn);

							//Apply left boundary conditions						
							auto bcMarker = zLeftBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[2] == pManager->dimsCart[2] - 1) {
							//Right border
							k = grid.kMax + layer; // layer index
							int kIn = grid.kMax - layer + 1; // opposite index
							cellCenter.z = grid.CoordinateZ[kIn];
							faceCenter.z = (grid.CoordinateZ[kIn] + grid.CoordinateZ[k]) / 2.0;
							int idx = grid.getSerialIndexLocal(i, j, k);
							int idxIn = grid.getSerialIndexLocal(i, j, kIn);

							//Apply right boundary conditions						
							auto bcMarker = zRightBC.getMarker(faceCenter);
							auto dValues = bConditions[bcMarker]->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};
					};
				};
			};
		}; //Z direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Dummy values Z-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	// Compute gradient of value in ijk cell
	Vector ComputeGradientNonUniformStructured(int i, int j, int k, std::function<double(double*)> func) {
		Vector grad(0, 0, 0);
		double dux = func(getCellValues(i + 1, j, k)) - func(getCellValues(i - 1, j, k));
		double dudx = dux / (grid.CoordinateX[i + 1] - grid.CoordinateX[i - 1]);
		grad.x = dudx;

		// 2D case
		if (nDims > 1) {
			double duy = func(getCellValues(i, j + 1, k)) - func(getCellValues(i, j - 1, k));
			double dudy = duy / (grid.CoordinateY[j + 1] - grid.CoordinateY[j - 1]);
			grad.y = dudy;
		};

		// 3D case
		if (nDims > 2) {
			double duz = func(getCellValues(i, j, k + 1)) - func(getCellValues(i, j, k - 1));
			double dudz = duz / (grid.CoordinateZ[k + 1] - grid.CoordinateZ[k - 1]);
			grad.z = dudz;
		};
		return grad;
	};

	// Compute source terms
	virtual void ProcessExternalForces() {
		//right handside part - mass foces
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					int idx = grid.getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = grid.hx[i] * grid.hy[j] * grid.hz[k];

					double *V = getCellValues(i, j, k);
					double u = V[1] / V[0];
					double v = V[2] / V[0];
					double w = V[3] / V[0];
					Vector velocity = Vector(u, v, w); //Cell velocity

													   //Compute total residual
					residual[idx * nVariables];										//ro
					residual[idx * nVariables + 1] += volume * Sigma.x;				//rou
					residual[idx * nVariables + 2] += volume * Sigma.y;				//rov
					residual[idx * nVariables + 3] += volume * Sigma.z;				//row
					residual[idx * nVariables + 4] += Sigma * velocity * volume;	//roE
				};
			};
		}; //end right hand side part
	};
	virtual void ProcessExternalAcceleration() {
		//right handside part - mass foces
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					int idx = grid.getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = grid.hx[i] * grid.hy[j] * grid.hz[k];
					double *V = getCellValues(i, j, k);
					double u = V[1] / V[0];
					double v = V[2] / V[0];
					double w = V[3] / V[0];
					Vector velocity = Vector(u, v, w); //Cell velocity

					//Mass forces represents body accelerations (per unit mass)
					double ro = V[0];
					residual[idx * nVariables + 1] += ro * volume * UniformAcceleration.x;				//rou
					residual[idx * nVariables + 2] += ro * volume * UniformAcceleration.y;				//rov
					residual[idx * nVariables + 3] += ro * volume * UniformAcceleration.z;				//row
					residual[idx * nVariables + 4] += ro * UniformAcceleration * velocity * volume;	//roE
				};
			};
		};
	};

	// Save and load solution functions
	void SaveSolution(std::string fname) {

		int rank = pManager->getRank();

		// write down the total time
		std::ofstream ofs;
		if (pManager->IsMaster()) {
			ofs.open(fname, std::ios::binary | std::ios::out);
			ofs.write(reinterpret_cast<char*>(&stepInfo.Time), sizeof stepInfo.Time);
			ofs.close();
		};

		// Wait for previous process to finish writing
		if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

		// Reopen file for writing and write rank
		ofs.open(fname, std::ios::binary | std::ios_base::app);
		ofs.write(reinterpret_cast<char*>(&rank), sizeof rank);

		// Write our solution
		for (int k = grid.kMin; k <= grid.kMax; k++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int i = grid.iMin; i <= grid.iMax; i++) {
					// write down all cell values
					auto val = getCellValues(i, j, k);
					ofs.write(reinterpret_cast<char*>(val), nVariables * sizeof(double) );
				};
			};
		};
		ofs.close();

		// Signal to next process to begin writing
		if (!pManager->IsLastNode()) {
			pManager->Signal(rank + 1);
		};
		//Syncronize
		pManager->Barrier();
		return;
	};
	void LoadSolution(std::string fname) {

		int rank = pManager->getRank();
		int Nproc = pManager->getProcessorNumber();

		// read down the total time
		std::ifstream ifs;
		ifs.open(fname, std::ios::binary | std::ios::in);
		ifs.read(reinterpret_cast<char*>(&stepInfo.Time), sizeof stepInfo.Time);

		// start to read file
		for (int np = 0; np < Nproc; np++) {

			// read the rank of process that write the part of solution
			int read_rank;
			ifs.read(reinterpret_cast<char*>(&read_rank), sizeof rank);

			// if rank the same as current - read the data
			if (rank == read_rank) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					for (int j = grid.jMin; j <= grid.jMax; j++) {
						for (int i = grid.iMin; i <= grid.iMax; i++) {
							// read all cell values
							auto val = getCellValues(i, j, k);
							ifs.read(reinterpret_cast<char*>(val), nVariables * sizeof(double));
						};
					};
				};
				// finish of data reading
				break;	
			} // if ranks are differrent
			else {
				double read_var;
				for (int i = 0; i < grid.nCellsLocal * nVariables; i++) ifs.read(reinterpret_cast<char*>(&read_var), sizeof read_var);
			};
		};

		ifs.close();

		//Syncronize
		pManager->Barrier();
		if (pManager->IsLastNode()) std::cout <<"rank " << rank << ": solution is loaded" << std::endl;
		return;
	};

	// Tecplot format saver
	virtual void SaveSolutionToTecplot(std::string fname) {
		//Tecplot version    
		int rank = pManager->getRank();

		//	1D tecplot style
		if(nDims == 1) {

			// Open the file
			if (pManager->IsMaster()) {
				//Create file
				std::ofstream ofs(fname);
				ofs << std::scientific;

				// Header				
				ofs << "VARIABLES = ";
				ofs << "\"" << "X" << "\" ";
				ofs << "\"" << "ro" << "\" ";
				ofs << "\"" << "u" << "\" ";
				ofs << "\"" << "v" << "\" ";
				ofs << "\"" << "w" << "\" ";
				ofs << "\"" << "P" << "\" ";
				ofs << "\"" << "e" << "\" ";
				ofs << std::endl;
				ofs.close();
			};
			pManager->Barrier();

			//Wait for previous process to finish writing
			if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

			//Reopen file for writing
			std::ofstream ofs(fname, std::ios_base::app);
			ofs << std::scientific;

			//Solution
			for (int i = grid.iMin; i <= grid.iMax; i++) {
				//Obtain cell data
				double x = grid.CoordinateX[i];
				double* U = getCellValues(i, grid.jMin, grid.kMin);
				double ro = U[0];
				double u = U[1] / ro;
				double v = U[2] / ro;
				double w = U[3] / ro;
				double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
				double P = (gas_prop.gamma - 1.0) * ro * e;

				//Write to file
				ofs << x << " ";
				ofs << ro << " ";
				ofs << u << " ";
				ofs << v << " ";
				ofs << w << " ";
				ofs << P << " ";
				ofs << e << " ";
				ofs << std::endl;
			};	//	end cycle
			ofs.close();

			// Signal to next process to begin writing
			if (!pManager->IsLastNode()) {
				pManager->Signal(rank + 1);
			};

			//Syncronize
			pManager->Barrier();

			return;
		};	//end if

		//	2D tecplot style
		if (nDims == 2) {
			// Open the file
			if (pManager->IsMaster()) {
				//Create file
				std::ofstream ofs(fname);
				ofs << std::scientific;

				// Header				
				ofs << "VARIABLES = ";
				ofs << R"("X" )";
				ofs << R"("Y" )";
				ofs << R"("Rho" )";
				ofs << R"("u" )";
				ofs << R"("v" )";
				ofs << R"("w" )";
				ofs << R"("P" )";
				ofs << R"("e")";
				ofs << R"("T")";
				ofs << std::endl;
				ofs << R"(ZONE T=")" << stepInfo.Time << R"(")";
				ofs << std::endl;
				ofs << "N=" << (grid.nX + 1) * (grid.nY + 1) << ", E=" << grid.nX * grid.nY << ", F=FEBLOCK, ET=QUADRILATERAL";
				ofs << std::endl;
				ofs << "VARLOCATION = (NODAL, NODAL";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED)";
				ofs << std::endl;

				// Write coordinates
				for (int j = grid.dummyCellLayersY; j <= grid.dummyCellLayersY + grid.nY; j++) {
					for (int i = grid.dummyCellLayersX; i <= grid.dummyCellLayersX + grid.nX; i++) {
						ofs << grid.CoordinateX[i] - 0.5 * grid.hx[i];
						ofs << std::endl;
					};
				};
				for (int j = grid.dummyCellLayersY; j <= grid.dummyCellLayersY + grid.nY; j++) {
					for (int i = grid.dummyCellLayersX; i <= grid.dummyCellLayersX + grid.nX; i++) {
						ofs << grid.CoordinateY[j] - 0.5 * grid.hy[j];
						ofs << std::endl;
					};
				};
				ofs.close();
			};
			pManager->Barrier();

			// Compute variables
			auto idx_first = grid.getSerialIndexLocal(grid.iMin, grid.jMin, 0);
			std::valarray<size_t> s{ size_t(grid.nlocalY), size_t(grid.nlocalX) };
			std::valarray<size_t> str{ size_t(grid.nlocalXAll * nVariables), size_t(nVariables) };

			// Density
			std::valarray<double> rho = values[std::gslice(idx_first * nVariables, s, str)];

			// X Velocity
			std::valarray<double> u = values[std::gslice(idx_first * nVariables + 1, s, str)];
			u /= rho;

			// Y Velocity
			std::valarray<double> v = values[std::gslice(idx_first * nVariables + 2, s, str)];
			v /= rho;

			// Z Velocity
			std::valarray<double> w = values[std::gslice(idx_first * nVariables + 3, s, str)];
			w /= rho;

			// Internal Energy and Pressure
			std::valarray<double> e = values[std::gslice(idx_first * nVariables + 4, s, str)];
			e /= rho;
			e = e - 0.5 * (u * u + v * v + w * w);
			std::valarray<double> P = (gas_prop.gamma - 1.0) * (rho * e);

			// Temperature
			double cond = (gas_prop.gamma - 1.0) * gas_prop.molarMass / gas_prop.universalGasConstant;
			std::valarray<double> T = cond * e;

			// write all data
			pManager->WriteData(rho, fname);
			pManager->WriteData(u, fname);
			pManager->WriteData(v, fname);
			pManager->WriteData(w, fname);
			pManager->WriteData(P, fname);
			pManager->WriteData(e, fname);
			pManager->WriteData(T, fname);

			if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

			// Reopen file for writing
			std::ofstream ofs(fname, std::ios_base::app);

			// write connectivity list
			int bl, br, ul, ur;		// indexes for left/right top/bottom nodes
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				bl = (grid.nX + 1) * (j - grid.dummyCellLayersY) + (grid.iMin -  grid.dummyCellLayersX) + 1;		// indexes starts from 1 in Tecplot format
				for (int i = grid.iMin; i <= grid.iMax; i++) {
					br = bl + 1;
					ul = bl + (grid.nX + 1);
					ur = ul + 1;
					ofs << bl << ' ';
					ofs << ul << ' ';
					ofs << ur << ' ';
					ofs << br;
					ofs << std::endl;
					bl++;
				};
			};
			ofs.close();

			if(!pManager->IsLastNode()) pManager->Signal(rank + 1);
			pManager->Barrier();
		};	//end 2D

		//	3D tecplot style
		if (nDims == 3) {
			// Open the file
			if (pManager->IsMaster()) {
				//Create file
				std::ofstream ofs(fname);
				ofs << std::scientific;

				// Header				
				ofs << "VARIABLES = ";
				ofs << R"("X" )";
				ofs << R"("Y" )";
				ofs << R"("Z" )";
				ofs << R"("Rho" )";
				ofs << R"("u" )";
				ofs << R"("v" )";
				ofs << R"("w" )";
				ofs << R"("P" )";
				ofs << R"("e")";
				ofs << R"("T")";
				ofs << std::endl;
				ofs << R"(ZONE T=")" << stepInfo.Time << R"(")";
				ofs << std::endl;
				ofs << "N=" << (grid.nX + 1) * (grid.nY + 1) * (grid.nZ + 1) << ", E=" << grid.nX * grid.nY * grid.nZ << ", F=FEBLOCK, ET=BRICK";
				ofs << std::endl;
				ofs << "VARLOCATION = (NODAL, NODAL, NODAL";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED)";
				ofs << std::endl;

				// Write coordinates
				for (int k = grid.dummyCellLayersZ; k <= grid.dummyCellLayersZ + grid.nZ; k++) {
					for (int j = grid.dummyCellLayersY; j <= grid.dummyCellLayersY + grid.nY; j++) {
						for (int i = grid.dummyCellLayersX; i <= grid.dummyCellLayersX + grid.nX; i++) {
							ofs << grid.CoordinateX[i] - 0.5 * grid.hx[i];
							ofs << std::endl;
						};
					};
				};
				for (int k = grid.dummyCellLayersZ; k <= grid.dummyCellLayersZ + grid.nZ; k++) {
					for (int j = grid.dummyCellLayersY; j <= grid.dummyCellLayersY + grid.nY; j++) {
						for (int i = grid.dummyCellLayersX; i <= grid.dummyCellLayersX + grid.nX; i++) {
							ofs << grid.CoordinateY[j] - 0.5 * grid.hy[j];
							ofs << std::endl;
						};
					};
				};
				for (int k = grid.dummyCellLayersZ; k <= grid.dummyCellLayersZ + grid.nZ; k++) {
					for (int j = grid.dummyCellLayersY; j <= grid.dummyCellLayersY + grid.nY; j++) {
						for (int i = grid.dummyCellLayersX; i <= grid.dummyCellLayersX + grid.nX; i++) {
							ofs << grid.CoordinateZ[k] - 0.5 * grid.hz[k];
							ofs << std::endl;
						};
					};
				};
				ofs.close();
			};
			pManager->Barrier();

			// Compute variables
			auto idx_first = grid.getSerialIndexLocal(grid.iMin, grid.jMin, grid.kMin);
			std::valarray<size_t> s{ size_t(grid.nlocalZ), size_t(grid.nlocalY), size_t(grid.nlocalX) };
			std::valarray<size_t> str{ size_t(grid.nlocalYAll * grid.nlocalXAll * nVariables), size_t(grid.nlocalXAll * nVariables), size_t(nVariables) };

			// Density
			std::valarray<double> rho = values[std::gslice(idx_first * nVariables, s, str)];

			// X Velocity
			std::valarray<double> u = values[std::gslice(idx_first * nVariables + 1, s, str)];		// ro * u

			// Y Velocity
			std::valarray<double> v = values[std::gslice(idx_first * nVariables + 2, s, str)];		// ro * v
			
			// Z Velocity
			std::valarray<double> w = values[std::gslice(idx_first * nVariables + 3, s, str)];		// ro * w

			// Mean kinetic energy
			std::valarray<double> k = 0.5 * (u * u + v * v + w * w) / rho;
			u /= rho;		// u
			v /= rho;		// v
			w /= rho;		// w
			
			// Internal Energy and Pressure
			std::valarray<double> e = values[std::gslice(idx_first * nVariables + 4, s, str)];		// ro * E
			e = e - k;		// ro * e
			std::valarray<double> P = (gas_prop.gamma - 1.0) * (e);
			e /= rho;		// e


			// Temperature
			double cond = (gas_prop.gamma - 1.0) * gas_prop.molarMass / gas_prop.universalGasConstant;
			std::valarray<double> T = cond * e;

			// Write all data
			pManager->WriteData(rho, fname);
			pManager->WriteData(u, fname);
			pManager->WriteData(v, fname);
			pManager->WriteData(w, fname);
			pManager->WriteData(P, fname);
			pManager->WriteData(e, fname);
			pManager->WriteData(T, fname);

			if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

			// Reopen file for writing
			std::ofstream ofs(fname, std::ios_base::app);

			// write connectivity list
			int fbl, fbr, ful, fur, bbl, bbr, bul, bur;		// indexes for front/back, top/bottom, left/right nodes
			for (int k = grid.kMin; k <= grid.kMax; k++) {
				for (int j = grid.jMin; j <= grid.jMax; j++) {
					fbl = (grid.nX + 1) * (grid.nY + 1) * (k - grid.dummyCellLayersZ) + (grid.nX + 1) * (j - grid.dummyCellLayersY) + (grid.iMin - grid.dummyCellLayersX) + 1;		// indexes starts from 1 in Tecplot format
					for (int i = grid.iMin; i <= grid.iMax; i++) {
						fbr = fbl + 1;
						ful = fbl + (grid.nX + 1);
						fur = ful + 1;
						bbl = fbl + (grid.nX + 1) * (grid.nY + 1);
						bbr = bbl + 1;
						bul = bbl + (grid.nX + 1);
						bur = bul + 1;
						ofs << bbr << ' ';
						ofs << bur << ' ';
						ofs << bul << ' ';
						ofs << bbl << ' ';
						ofs << fbr << ' ';
						ofs << fur << ' ';
						ofs << ful << ' ';
						ofs << fbl;
						ofs << std::endl;
						fbl++;
					};
				};
			};
			ofs.close();

			if (!pManager->IsLastNode()) pManager->Signal(rank + 1);
			pManager->Barrier();
		};	//end 2D

		return;
	};

	// Save 2D slice
	void Save2DSliceToTecplot(std::string fname, int I, int J, int K) {

		std::ofstream ofs;
		ofs << std::scientific;
		int rank = pManager->getRank();

		//Open file
		if (pManager->IsMaster()) {
			//std::cout<<"File created"<<std::endl<<std::flush;
			ofs.open(fname, std::ios_base::out);
			ofs << "VARIABLES = ";
			int dims = 0;
			if (I == -1) {
				ofs << "\"" << "X" << "\" ";
				dims++;
			};
			if (J == -1) {
				ofs << "\"" << "Y" << "\" ";
				dims++;
			};
			if (K == -1) {
				ofs << "\"" << "Z" << "\" ";
				dims++;
			};

			ofs << "\"" << "ro" << "\" ";
			ofs << "\"" << "u" << "\" ";
			ofs << "\"" << "v" << "\" ";
			ofs << "\"" << "w" << "\" ";
			ofs << "\"" << "P" << "\" ";
			ofs << "\"" << "e" << "\" ";
			ofs << std::endl;

			ofs << "ZONE T=\"1\"";	//zone name
			ofs << std::endl;

			if (I == -1) ofs << "I=" << grid.nX;
			else ofs << "I=" << 1;
			if (J == -1) ofs << "J=" << grid.nY;
			else ofs << "J=" << 1;
			if (K == -1) ofs << "K=" << grid.nZ;
			else ofs << "K=" << 1;
			ofs << "F=POINT\n";
			ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
		}
		else {
			//Wait for previous process to finish writing
			pManager->Wait(rank - 1);

			//Open file for modification
			//std::cout<<"File opened"<<std::endl<<std::flush;
			ofs.open(fname, std::ios_base::app);
		};

		//Solution output
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			if ((I != -1) && (i != I)) continue;
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				if ((J != -1) && (j != J)) continue;
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					if ((K != -1) && (k != K)) continue;
					//Obtain cell data
					double x = grid.CoordinateX[i];
					double y = grid.CoordinateY[j];
					double z = grid.CoordinateZ[k];
					double* U = getCellValues(i, j, k);
					double ro = U[0];
					double u = U[1] / ro;
					double v = U[2] / ro;
					double w = U[3] / ro;
					double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
					double P = (gas_prop.gamma - 1.0) * ro * e;

					//Write to file
					if (I == -1) ofs << x << " ";
					if (J == -1) ofs << y << " ";
					if (K == -1) ofs << z << " ";
					ofs << ro << " ";
					ofs << u << " ";
					ofs << v << " ";
					ofs << w << " ";
					ofs << P << " ";
					ofs << e << " ";
					ofs << std::endl;
				};
			};
		}; // Solution output

		   //Close file and signal to next process to begin writing
		ofs.close();
		//std::cout<<"File closed"<<std::endl<<std::flush;
		if (rank != pManager->getProcessorNumber() - 1) {
			pManager->Signal(rank + 1);
		};

		//Syncronize
		pManager->Barrier();
	};

	// Save 1D Slice
	void SaveSliceToTecplot(std::string fname, Slice &s) {
		int I = s.i;
		int J = s.j;
		int K = s.k;

		std::ofstream ofs;
		ofs << std::scientific;
		int rank = pManager->getRank();

		//Open file
		if (pManager->IsMaster()) {
			ofs.open(fname, std::ios_base::out);
			ofs << "VARIABLES = ";
			int dims = 0;
			if (I == -1) {
				ofs << "\"" << "X" << "\" ";
				dims++;
			};
			if (J == -1) {
				ofs << "\"" << "Y" << "\" ";
				dims++;
			};
			if (K == -1) {
				ofs << "\"" << "Z" << "\" ";
				dims++;
			};

			ofs << "\"" << "ro" << "\" ";
			ofs << "\"" << "u" << "\" ";
			ofs << "\"" << "v" << "\" ";
			ofs << "\"" << "w" << "\" ";
			ofs << "\"" << "P" << "\" ";
			ofs << "\"" << "e" << "\" ";
			ofs << "\"" << "T" << "\" ";

			if (I != -1) {
				ofs << "\"" << "x_cnst" << "\" ";
				dims++;
			};
			if (J != -1) {
				ofs << "\"" << "y_cnst" << "\" ";
				dims++;
			};
			if (K != -1) {
				ofs << "\"" << "z_cnst" << "\" ";
				dims++;
			};
			ofs << std::endl;

			ofs << "ZONE T=\"1\"";	//zone name
			ofs << std::endl;

			if (I == -1) ofs << "I=" << grid.nX;
			else ofs << "I=" << 1;
			if (J == -1) ofs << "J=" << grid.nY;
			else ofs << "J=" << 1;
			if (K == -1) ofs << "K=" << grid.nZ;
			else ofs << "K=" << 1;
			ofs << "F=POINT\n";
			ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
		}
		else {
			//Wait for previous process to finish writing
			pManager->Wait(rank - 1);

			ofs.open(fname, std::ios_base::app);
		};

		//Solution output
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			if ((I != -1) && (i != I)) continue;
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				if ((J != -1) && (j != J)) continue;
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					if ((K != -1) && (k != K)) continue;
					//Obtain cell data
					double x = grid.CoordinateX[i];
					double y = grid.CoordinateY[j];
					double z = grid.CoordinateZ[k];
					double* U = getCellValues(i, j, k);
					double ro = U[0];
					double u = U[1] / ro;
					double v = U[2] / ro;
					double w = U[3] / ro;
					double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
					double P = (gas_prop.gamma - 1.0) * ro * e;
					double T = (gas_prop.gamma - 1.0) * gas_prop.molarMass * e / gas_prop.universalGasConstant;

					//Write to file
					if (I == -1) ofs << x << " ";
					if (J == -1) ofs << y << " ";
					if (K == -1) ofs << z << " ";
					ofs << ro << " ";
					ofs << u << " ";
					ofs << v << " ";
					ofs << w << " ";
					ofs << P << " ";
					ofs << e << " ";
					ofs << T << " ";
					if (I != -1) ofs << x << " ";
					if (J != -1) ofs << y << " ";
					if (K != -1) ofs << z << " ";
					ofs << std::endl;
				};
			};
		}; // Solution output

		   //Close file and signal to next process to begin writing
		ofs.close();
		//std::cout<<"File closed"<<std::endl<<std::flush;
		if (rank != pManager->getProcessorNumber() - 1) {
			pManager->Signal(rank + 1);
		};

		//Syncronize
		pManager->Barrier();
	};

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int SolutionTimePrecision = 0;
		int SliceTimePrecision = 0;
		if (SaveSolutionTime > 0) {
			SolutionTimePrecision = static_cast<int>(1 - std::floor(std::log10(SaveSolutionTime)));
		};
		if (SaveSliceTime > 0) {
			SliceTimePrecision = static_cast<int>(1 - std::floor(std::log10(SaveSliceTime)));
		};

		//Open history file
		std::ofstream ofs;
		if (pManager->IsMaster()) {
			ofs.open("history.dat", std::ios_base::out);
			//Header
			ofs << "VARIABLES = ";
			ofs << "\"" << "time" << "\" ";
			ofs << "\"" << "roR" << "\" ";
			ofs << "\"" << "rouR" << "\" ";
			if (nDims > 1) ofs << "\"" << "rovR" << "\" ";
			if (nDims > 2) ofs << "\"" << "rowR" << "\" ";
			if (nVariables > 4) ofs << "\"" << "roeR" << "\" ";
			ofs << std::endl;
		};

		//Start timers
		Timer workTimer;
		workTimer.Start();
		workTimer.Resume();

		Timer iterTimer;
		iterTimer.Start();

		//Calc loop		
		if (pManager->IsMaster()) {
			std::cout << "Calculation started!" << std::endl << std::flush;
		};
		pManager->Barrier();

		//Processor info output
		rank = pManager->getRank();
		if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);
		std::cout << "Info for node rank : " << rank << "" << std::endl <<
			"rankX = " << pManager->rankCart[0] <<
			", iMin = " << grid.iMin <<
			", iMax = " << grid.iMax
			<< std::endl;
		if (nDims > 1) {
			std::cout <<
				"rankY = " << pManager->rankCart[1] <<
				", jMin = " << grid.jMin <<
				", jMax = " << grid.jMax
				<< std::endl;
		};
		if (nDims > 2) {
			std::cout <<
				"rankZ = " << pManager->rankCart[2] <<
				", kMin = " << grid.kMin <<
				", kMax = " << grid.kMax
				<< std::endl;
		};
		if (!pManager->IsLastNode()) pManager->Signal(rank + 1);
		pManager->Barrier();

		while(true) {
			//Calculate one time step
			IterationStep();

			//Output step information						
			if (pManager->IsMaster() && (ResidualOutputIterations != 0) && (stepInfo.Iteration % ResidualOutputIterations) == 0) {
				iterTimer.Stop();
				double iterationsPackTimeAverage = iterTimer.ElapsedTimeMilliseconds();
				iterationsPackTimeAverage /= ResidualOutputIterations;
				iterTimer.Start();
				std::cout << "Iteration = " << stepInfo.Iteration <<
					"; Total time = " << stepInfo.Time <<
					"; Time step = " << stepInfo.TimeStep <<
					"; RMSrou = " << stepInfo.Residual[1] <<
					"; avgIterTime = " << iterationsPackTimeAverage << " ms" <<
					std::endl;
			};

			//Sensors recording
			if ((isSensorEnable == true) && (stepInfo.Iteration % SaveSensorRecordIterations == 0))
			{
				for (auto &r : Sensors) {
					r->UpdateIteration(stepInfo.Iteration);
					r->UpdateTimer(stepInfo.Time);
					r->Process(values);
				};
			};

			//Solution and slice snapshots every few iterations
			if ((SaveSolutionIterations != 0) && (stepInfo.Iteration % SaveSolutionIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				//snapshotFileName.str(std::string());
				snapshotFileName << "solution, It = " << stepInfo.Iteration << ".dat";
				SaveSolutionToTecplot(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout << "Solution has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
				};
			};
			if ((SaveSliceIterations != 0) && (stepInfo.Iteration % SaveSliceIterations) == 0) {
				//Save snapshots
				for (auto i = 0; i < slices.size(); i++) {
					std::stringstream snapshotFileName;
					//snapshotFileName.str(std::string());
					snapshotFileName << "slice" << i << ", It =" << stepInfo.Iteration << ".dat";
					SaveSliceToTecplot(snapshotFileName.str(), slices[i]);

					if (pManager->IsMaster()) {
						std::cout << "Slice has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
					};
				};
			};
			
			//Every fixed time interval
			if ((SaveSolutionTime > 0) && (stepInfo.NextSolutionSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				//snapshotFileName.str(std::string());
				snapshotFileName << std::fixed;
				snapshotFileName.precision(SolutionTimePrecision);
				snapshotFileName << "solution, t = " << stepInfo.Time << ".dat";
				SaveSolutionToTecplot(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout << "Solution has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
				};

				//Adjust next snapshot time
				stepInfo.NextSolutionSnapshotTime += SaveSolutionTime;
			};
			if ((SaveSliceTime > 0) && (stepInfo.NextSliceSnapshotTime == stepInfo.Time)) {
				//Save snapshots
				for (auto i = 0; i < slices.size(); i++) {
					std::stringstream snapshotFileName;
					//snapshotFileName.str(std::string());
					snapshotFileName << "slice" << i;
					snapshotFileName << std::fixed;
					snapshotFileName.precision(SliceTimePrecision);
					snapshotFileName << ",T=" << stepInfo.Time << ".dat";
					SaveSliceToTecplot(snapshotFileName.str(), slices[i]);

					if (pManager->IsMaster()) {
						std::cout << "Slice has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
					};
				};

				//Adjust next snapshot time
				stepInfo.NextSliceSnapshotTime += SaveSliceTime;
			};

			//Save convergence history		
			if (pManager->IsMaster() && (stepInfo.Iteration % ResidualOutputIterations)) {
				ofs << stepInfo.Time << " ";
				ofs << stepInfo.Residual[0] << " ";
				ofs << stepInfo.Residual[1] << " ";
				if (nDims > 1) ofs << stepInfo.Residual[2] << " ";
				if (nDims > 2) ofs << stepInfo.Residual[3] << " ";
				if (nVariables > 4) ofs << stepInfo.Residual[4] << " ";
				ofs << std::endl;
			};

			//Convergence criteria
			if (stepInfo.Iteration == MaxIteration) {
				if (pManager->IsMaster()) {
					std::cout << "Maximal number of iterations reached.\n";
				};
				break;
			};

			if (stepInfo.Time >= MaxTime) {
				if (pManager->IsMaster()) {
					std::cout << "Maximal time reached.\n";
				};
				break;
			};

			//Synchronize
			pManager->Barrier();
		};

		//Synchronize		
		pManager->Barrier();

		//Calculate max idle time
		double idleTime = pManager->getIdleTime<std::milli>(); //ms
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", local idle time = " << idleTime / 1e3 << " seconds.\n";
			std::cout.flush();
		};

		double maxIdleTime = pManager->Max(idleTime);

		if (pManager->IsMaster()) {
			//Close history file
			ofs.flush();
			ofs.close();

			std::cout << "Calculation finished!\n";
			//Stop timer
			workTimer.Stop();
			std::cout << "Total work time = " << workTimer.ElapsedTimeMilliseconds() / 1e3 << " seconds.\n";
			std::cout << "Maximum idle time = " << maxIdleTime / 1e3 << " seconds.\n";
		};

	};

	//Finalize kernel
	void Finalize() {
		//Finalize MPI
		pManager->Finalize();
	};
	
};



#endif
