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

#include "KernelConfiguration.h"
#include "grid.h"
#include "ParallelManager.h"
#include "utility/Vector.h"
#include "utility/Matrix.h"
#include "utility/Timer.h"
#include "RiemannSolvers/RiemannSolversList.h"
#include "BoundaryConditions/BCGeneral.h"
#include "Sensors/Sensors.h"

// #include "cgnslib.h"

// Step info
class StepInfo {
public:
	double Time;
	double TimeStep;	
	int Iteration;
	std::vector<double> Residual;
	double NextSnapshotTime;
};


// Calculation kernel
class Kernel {
public:
	virtual void IterationStep() = 0; // main function

	// Parallel run information
	std::unique_ptr<ParallelManager> pManager;
	int rank; //Rank of process

	// Structured grid and dimensions number
	int nDims;
	Grid g;

	// Current step information
	StepInfo stepInfo;

	// Gas model information
	int nVariables; // number of conservative variables
	double gamma;
	double thermalConductivity;
	double viscosity;
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

	// Boundary conditions
	std::unique_ptr<BoundaryConditions::BCGeneral> xLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> xRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zRightBC;

	// Solution data
	std::valarray<double> values;
	std::valarray<double> residual;

	// Get serial index for cell
	inline int getSerialIndexGlobal(int i, int j, int k) {
		int sI = (k * g.nXAll * g.nYAll + j * g.nXAll + i);
		return sI;
	};

	// Get local index for cell
	inline int getSerialIndexLocal(int i, int j, int k) {
		int sI = (k - g.kMin + g.dummyCellLayersZ) * g.nlocalXAll * g.nlocalYAll + (j - g.jMin + g.dummyCellLayersY) * g.nlocalXAll + (i - g.iMin + g.dummyCellLayersX);
		return sI;
	};

	// Get cell values
	inline double* getCellValues(int i, int j, int k) {
		int sBegin = getSerialIndexLocal(i, j, k) * nVariables;
		return &values[sBegin];
	};

	// Transition beatween variables (trivial by default)
	inline std::valarray<double> ForvardVariablesTransition(double* vals) {
		std::valarray<double> res(nVariables);
		double rho = vals[0];
		double u = vals[1] / rho;
		double v = vals[2] / rho;
		double w = vals[3] / rho;
		double rho_e = vals[4] - 0.5 * rho * (u * u + v * v + w * w);

		// Pressure
		res[0] = rho_e * (gamma - 1.0);

		// Velosities
		res[1] = u;
		res[2] = v;
		res[3] = w;

		// Internal energy
		res[4] = rho_e / rho;

		return std::move(res);
	};

	// Inverse transition
	inline std::valarray<double> InverseVariablesTransition(std::valarray<double> &vals) {
		std::valarray<double> res(nVariables);
		double rho_e = vals[0] / (gamma - 1.0);
		double rho = rho_e / vals[4];
		double rho_u = rho * vals[1];
		double rho_v = rho * vals[2];
		double rho_w = rho * vals[3];
		double rho_E = rho_e + 0.5 * (rho_u * rho_u + rho_v * rho_v + rho_w * rho_w) / rho;

		//fill conservative result array
		res[0] = rho;
		res[1] = rho_u;
		res[2] = rho_v;
		res[3] = rho_w;
		res[4] = rho_E;

		return std::move(res);
	};

	// Calculation parameters	
	int MaxIteration;
	double MaxTime;
	double SaveSolutionSnapshotTime;
	int SaveSolutionSnapshotIterations;
	int ResidualOutputIterations;
	bool ContinueComputation;
	bool DebugOutputEnabled;

	//Array of sensors in use
	std::vector<std::unique_ptr<Sensor>> Sensors;

	//Constructor
	Kernel(int* argc, char **argv[]) : pManager(new ParallelManager(argc, argv)) {
		nVariables = 5;	//default value
		rank = pManager->getRank();
	};

	//Set initial conditions
	void SetInitialConditions(std::function<std::vector<double>(const Vector& r)> initF) {
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					//Obtain cell data
					double x = g.CoordinateX[i];
					double y = g.CoordinateY[j];
					double z = g.CoordinateZ[k];
					int sBegin = getSerialIndexLocal(i, j, k) * nVariables;
					std::vector<double> U = initF(Vector(x, y, z));
					for (int varn = 0; varn < nVariables; varn++) values[sBegin + varn] = U[varn];

        //  std::cout << "rank = " << rank << ", x = " << x << ", i = " << i << std::endl << std::flush; //MPI DebugMessage
        //  std::cout << "rank = " << rank << ", sBegin = " << sBegin << ", ro = " << values[sBegin] << ", rou = " << values[sBegin + 1] << std::endl << std::flush; //MPI DebugMessage
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

	//Initialize kernel
	virtual void Init(KernelConfiguration& config) {
		// Initialize MPI		
		nDims = config.nDims;

		// Output settings 
		DebugOutputEnabled = config.DebugOutputEnabled;

		// Initialize global grid
		g.InitGlobal(config);

		// Initialize cartezian topology
		pManager->InitCartesianTopology(g);

		// Compute local gris sizes
		g.nlocalX = g.nX / pManager->dimsCart[0];
		g.nlocalY = g.nY / pManager->dimsCart[1];
		g.nlocalZ = g.nZ / pManager->dimsCart[2];
		g.iMin = pManager->rankCart[0] * g.nlocalX + g.dummyCellLayersX;
		g.iMax = (pManager->rankCart[0] + 1) * g.nlocalX + g.dummyCellLayersX - 1;
		g.jMin = pManager->rankCart[1] * g.nlocalY + g.dummyCellLayersY;
		g.jMax = (pManager->rankCart[1] + 1) * g.nlocalY + g.dummyCellLayersY - 1;
		g.kMin = pManager->rankCart[2] * g.nlocalZ + g.dummyCellLayersZ;
		g.kMax = (pManager->rankCart[2] + 1) * g.nlocalZ + g.dummyCellLayersZ - 1;
		g.InitLocal(config);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->_rankCart << ", iMin = " << g.iMin << ", iMax = " << g.iMax << "\n";
			std::cout << "rank = " << pManager->_rankCart << ", jMin = " << g.jMin << ", jMax = " << g.jMax << "\n";
			std::cout << "rank = " << pManager->_rankCart << ", kMin = " << g.kMin << ", kMax = " << g.kMax << "\n";
		};
		
		//Sync
		pManager->Barrier();

		//Initialize gas model parameters and Riemann solver
		gamma = config.Gamma;
		viscosity = config.Viscosity;
		thermalConductivity = config.ThermalConductivity;
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
		values.resize(nVariables * g.nlocalXAll * g.nlocalYAll * g.nlocalZAll);
		residual.resize(nVariables * g.nlocalXAll * g.nlocalYAll * g.nlocalZAll);

		//Initialize calculation parameters
		MaxTime = config.MaxTime;
		MaxIteration = config.MaxIteration;
		SaveSolutionSnapshotTime = config.SaveSolutionSnapshotTime;
		SaveSolutionSnapshotIterations = config.SaveSolutionSnapshotIterations;
		ResidualOutputIterations = config.ResidualOutputIterations;
		stepInfo.Time = 0;
		stepInfo.Iteration = 0;
		stepInfo.NextSnapshotTime = stepInfo.Time + SaveSolutionSnapshotTime;

		//Initialize boundary conditions
		InitBoundaryConditions(config);

		//External forces
		Sigma = config.Sigma;
		UniformAcceleration = config.UniformAcceleration;

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Initialize boundary conditions
	virtual void InitBoundaryConditions(KernelConfiguration& config) {
		if (!g.IsPeriodicX) {
			xLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!g.IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};
		if ((!g.IsPeriodicZ) && (nDims > 2)) {
			zLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			zRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			zLeftBC->loadConfiguration(config.zLeftBoundary);
			zRightBC->loadConfiguration(config.zRightBoundary);
		};
	};

	//Update solution
	void UpdateSolution(double dt) {
		stepInfo.Residual.resize(nVariables);
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = 0;

		for (int i = g.iMin; i <= g.iMax; i++)
		{
			for (int j = g.jMin; j <= g.jMax; j++)
			{
				for (int k = g.kMin; k <= g.kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume{ g.hx[i] * g.hy[j] * g.hz[k] };
					//Update cell values
					for (int nv = 0; nv < nVariables; nv++) {
						values[idx * nVariables + nv] += residual[idx * nVariables + nv] * dt / volume;
						//Compute total residual
						stepInfo.Residual[nv] += abs(residual[idx * nVariables + nv]);			// L1 norm
					};
				};
			};
		};

		//Aggregate
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = pManager->Sum(stepInfo.Residual[nv]);
	};

	//Update Residual
	void UpdateResiduals() {
		stepInfo.Residual.resize(nVariables);
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = 0;

		for (int i = g.iMin; i <= g.iMax; i++)
		{
			for (int j = g.jMin; j <= g.jMax; j++)
			{
				for (int k = g.kMin; k <= g.kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = g.hx[i] * g.hy[j] * g.hz[k];
					//Update cell values
					for (int nv = 0; nv < nVariables; nv++) {
						//Compute total residual
						stepInfo.Residual[nv] += abs(residual[idx * nVariables + nv]);
					};
				};
			};
		};

		//Aggregate
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = pManager->Sum(stepInfo.Residual[nv]);
	};

	//Exchange values between processors
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
		int layerSize = g.nlocalY * g.nlocalZ;
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

		for (int layer = 1; layer <= g.dummyCellLayersX; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int iSend = g.iMin + layer - 1; // layer index to send
			int iRecv = g.iMax + layer; // layer index to recv

									  // Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (j = g.jMin; j <= g.jMax; j++) {
					for (k = g.kMin; k <= g.kMax; k++) {
						int idxBuffer = (j - g.jMin) + (k - g.kMin)* g.nlocalY; //Exclude x index
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
			for (j = g.jMin; j <= g.jMax; j++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					int idxBuffer = (j - g.jMin) + (k - g.kMin)* g.nlocalY; //Exclude x index
					double* U = getCellValues(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) U[nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			iSend = g.iMax - layer + 1; // layer index to send
			iRecv = g.iMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; iSend = " << iSend <<
					"; iRecv = " << iRecv <<
					std::endl << std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (j = g.jMin; j <= g.jMax; j++) {
					for (k = g.kMin; k <= g.kMax; k++) {
						int idxBuffer = (j - g.jMin) + (k - g.kMin)* g.nlocalY; //Exclude x index
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
			for (j = g.jMin; j <= g.jMax; j++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					int idxBuffer = (j - g.jMin) + (k - g.kMin)* g.nlocalY; //Exclude x index
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
		layerSize = g.nlocalX * g.nlocalZ;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);
		/*std::cout<<"rank = "<<rank<<
		"; rankL = "<<rankL<<
		"; rankR = "<<rankR<<
		std::endl<<std::flush;*/

		for (int layer = 1; layer <= g.dummyCellLayersY; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int jSend = g.jMin + layer - 1; // layer index to send
			int jRecv = g.jMax + layer; // layer index to recv

									  // Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = g.iMin; i <= g.iMax; i++) {
					for (k = g.kMin; k <= g.kMax; k++) {
						int idxBuffer = (i - g.iMin) + (k - g.kMin) * g.nlocalX; //Exclude y index
						int idxValues = getSerialIndexLocal(i, jSend, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = values[idxValues * nVariables + nv];
					};
				};
			};

			//Determine recive number
			if (rankR != -1) {
				nRecv = layerSize * nVariables;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankL, rankR, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write to recieving layer of dummy cells
			for (i = g.iMin; i <= g.iMax; i++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					int idxBuffer = (i - g.iMin) + (k - g.kMin) * g.nlocalX; //Exclude y index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			jSend = g.jMax - layer + 1; // layer index to send
			jRecv = g.jMin - layer; // layer index to recv

								  // Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = g.iMin; i <= g.iMax; i++) {
					for (k = g.kMin; k <= g.kMax; k++) {
						int idxBuffer = (i - g.iMin) + (k - g.kMin) * g.nlocalX; //Exclude y index
						int idxValues = getSerialIndexLocal(i, jSend, k);
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
			for (i = g.iMin; i <= g.iMax; i++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					int idxBuffer = (i - g.iMin) + (k - g.kMin) * g.nlocalX; //Exclude y index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
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
		layerSize = g.nlocalX * g.nlocalY;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, 0, -1);
		rankR = pManager->GetRankByCartesianIndexShift(0, 0, +1);
		//rankL = pManager->GetRankByCartesianIndex(0, 0, -1);
		//rankR = pManager->GetRankByCartesianIndex(0, 0, 1);

		if (DebugOutputEnabled) {
			std::cout << "Exchange Z-direction started" << std::endl << std::flush;
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};


		for (int layer = 1; layer <= g.dummyCellLayersZ; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int kSend = g.kMin + layer - 1; // layer index to send
			int kRecv = g.kMax + layer; // layer index to recv

									  // Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = g.iMin; i <= g.iMax; i++) {
					for (j = g.jMin; j <= g.jMax; j++) {
						int idxBuffer = (i - g.iMin) + (j - g.jMin) * g.nlocalX; //Exclude z index
						int idxValues = getSerialIndexLocal(i, j, kSend);
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
			for (i = g.iMin; i <= g.iMax; i++) {
				for (j = g.jMin; j <= g.jMax; j++) {
					int idxBuffer = (i - g.iMin) + (j - g.jMin) * g.nlocalX; //Exclude z index
					int idxValues = getSerialIndexLocal(i, j, kRecv);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			kSend = g.kMax - layer + 1; // layer index to send
			kRecv = g.kMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout << "rank = " << rank <<
					"; kSend = " << kSend <<
					"; kRecv = " << kRecv <<
					std::endl << std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = g.iMin; i <= g.iMax; i++) {
					for (j = g.jMin; j <= g.jMax; j++) {
						int idxBuffer = (i - g.iMin) + (j - g.jMin) * g.nlocalX; //Exclude z index
						int idxValues = getSerialIndexLocal(i, j, kSend);
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
			for (i = g.iMin; i <= g.iMax; i++) {
				for (j = g.jMin; j <= g.jMax; j++) {
					int idxBuffer = (i - g.iMin) + (j - g.jMin) * g.nlocalX; //Exclude z index
					int idxValues = getSerialIndexLocal(i, j, kRecv);
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
			std::cout << "rank = " << rank << "Exchange values finished." <<
				std::endl << std::flush;
		};

	}; // function

	//Compute dummy cell values as result of boundary conditions and interprocessor exchange communication
	void ComputeDummyCellValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank << "Dummy cell processing started." <<
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
		if (!g.IsPeriodicX) {
			for (j = g.jMin; j <= g.jMax; j++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					//Inner cell
					cellCenter.y = g.CoordinateY[j];
					cellCenter.z = g.CoordinateZ[k];

					//And face
					faceCenter.y = g.CoordinateY[j];
					faceCenter.z = g.CoordinateZ[k];

					if (DebugOutputEnabled) {
						std::cout << "rank = " << rank <<
							"; f.y = " << faceCenter.y <<
							"; f.z = " << faceCenter.z <<
							std::endl << std::flush;
					};

					for (int layer = 1; layer <= g.dummyCellLayersX; layer++) {
						if (pManager->rankCart[0] == 0) {
							//Left border
							i = g.iMin - layer; // layer index
							int iIn = g.iMin + layer - 1; // opposite index
							cellCenter.x = g.CoordinateX[iIn];
							faceCenter.x = (g.CoordinateX[iIn] + g.CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);

							//Apply left boundary conditions						
							std::valarray<double> dValues = xLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[0] == pManager->dimsCart[0] - 1) {
							//Right border
							i = g.iMax + layer; // layer index
							int iIn = g.iMax - layer + 1; // opposite index
							cellCenter.x = g.CoordinateX[iIn];
							faceCenter.x = (g.CoordinateX[iIn] + g.CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);

							//Apply right boundary conditions						
							std::valarray<double> dValues = xRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
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
		if (!g.IsPeriodicY) {
			for (i = g.iMin; i <= g.iMax; i++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					//Inner cell
					cellCenter.x = g.CoordinateX[i];
					cellCenter.z = g.CoordinateZ[k];

					//And face
					faceCenter.x = g.CoordinateX[i];
					faceCenter.z = g.CoordinateZ[k];

					for (int layer = 1; layer <= g.dummyCellLayersY; layer++) {
						if (pManager->rankCart[1] == 0) {
							//Left border
							j = g.jMin - layer; // layer index
							int jIn = g.jMin + layer - 1; // opposite index
							cellCenter.y = g.CoordinateY[jIn];
							faceCenter.y = (g.CoordinateY[jIn] + g.CoordinateY[j]) / 2.0;	//for uniform dummy cells		TO DO MODIFY
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);

							//Apply left boundary conditions						
							std::valarray<double> dValues = yLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[1] == pManager->dimsCart[1] - 1) {
							//Right border
							j = g.jMax + layer; // layer index
							int jIn = g.jMax - layer + 1; // opposite index
							cellCenter.y = g.CoordinateY[jIn];
							faceCenter.y = (g.CoordinateY[jIn] + g.CoordinateY[j]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);

							//Apply right boundary conditions						
							std::valarray<double> dValues = yRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
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
		if (!g.IsPeriodicZ) {
			for (i = g.iMin; i <= g.iMax; i++) {
				for (j = g.jMin; j <= g.jMax; j++) {
					//Inner cell
					cellCenter.x = g.CoordinateX[i];
					cellCenter.y = g.CoordinateY[j];

					//And face
					faceCenter.x = g.CoordinateX[i];
					faceCenter.y = g.CoordinateY[j];

					for (int layer = 1; layer <= g.dummyCellLayersY; layer++) {

						if (pManager->rankCart[2] == 0) {
							//Left border
							k = g.kMin - layer; // layer index
							int kIn = g.kMin + layer - 1; // opposite index
							cellCenter.z = g.CoordinateZ[kIn];
							faceCenter.z = (g.CoordinateZ[kIn] + g.CoordinateZ[k]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, j, kIn);

							//Apply left boundary conditions						
							std::valarray<double> dValues = zLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[2] == pManager->dimsCart[2] - 1) {
							//Right border
							k = g.kMax + layer; // layer index
							int kIn = g.kMax - layer + 1; // opposite index
							cellCenter.z = g.CoordinateZ[kIn];
							faceCenter.z = (g.CoordinateZ[kIn] + g.CoordinateZ[k]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, j, kIn);

							//Apply right boundary conditions						
							std::valarray<double> dValues = zRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
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

	//Compute source terms
	virtual void ProcessExternalForces() {
		//right handside part - mass foces
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = g.hx[i] * g.hy[j] * g.hz[k];

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
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = g.hx[i] * g.hy[j] * g.hz[k];
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

	virtual void SaveSolution(std::string fname) {
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
			for (int i = g.iMin; i <= g.iMax; i++) {
				//Obtain cell data
				double x = g.CoordinateX[i];
				double* U = getCellValues(i, g.jMin, g.kMin);
				double ro = U[0];
				double u = U[1] / ro;
				double v = U[2] / ro;
				double w = U[3] / ro;
				double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
				double P = (gamma - 1.0) * ro * e;

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

		//2D tecplot style
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
				ofs << std::endl;
				ofs << R"(ZONE T=")" << stepInfo.Time << R"(")";
				ofs << std::endl;
				ofs << "N=" << (g.nX + 1) * (g.nY + 1) << ", E=" << g.nX * g.nY << ", F=FEBLOCK, ET=QUADRILATERAL";
				ofs << std::endl;
				ofs << "VARLOCATION = (NODAL, NODAL";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED)";
				ofs << std::endl;

				// Write coordinates
				for (int j = g.dummyCellLayersY; j <= g.dummyCellLayersY + g.nY; j++) {
					for (int i = g.dummyCellLayersX; i <= g.dummyCellLayersX + g.nX; i++) {
						ofs << g.CoordinateX[i] - 0.5 * g.hx[i];
						ofs << std::endl;
					};
				};
				for (int j = g.dummyCellLayersY; j <= g.dummyCellLayersY + g.nY; j++) {
					for (int i = g.dummyCellLayersX; i <= g.dummyCellLayersX + g.nX; i++) {
						ofs << g.CoordinateY[j] - 0.5 * g.hy[j];
						ofs << std::endl;
					};
				};
				ofs.close();
			};
			pManager->Barrier();

			// Compute variables
			auto idx_first = getSerialIndexLocal(g.iMin, g.jMin, 0);
			std::valarray<size_t> s{ size_t(g.nlocalY), size_t(g.nlocalX) };
			std::valarray<size_t> str{ size_t(g.nlocalXAll * nVariables), size_t(nVariables) };

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
			std::valarray<double> P = (gamma - 1.0) * (rho * e);

			// write all data
			pManager->WriteData(rho, fname);
			pManager->WriteData(u, fname);
			pManager->WriteData(v, fname);
			pManager->WriteData(w, fname);
			pManager->WriteData(P, fname);
			pManager->WriteData(e, fname);

			if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

			// Reopen file for writing
			std::ofstream ofs(fname, std::ios_base::app);

			// write connectivity list
			int bl, br, ul, ur;		// indexes for left/right top/bottom nodes
			for (int j = g.jMin; j <= g.jMax; j++) {
				bl = (g.nX + 1) * (j - g.dummyCellLayersY) + (g.iMin -  g.dummyCellLayersX) + 1;		// indexes starts from 1 in Tecplot format
				for (int i = g.iMin; i <= g.iMax; i++) {
					br = bl + 1;
					ul = bl + (g.nX + 1);
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

		//3D tecplot style
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
				ofs << std::endl;
				ofs << R"(ZONE T=")" << stepInfo.Time << R"(")";
				ofs << std::endl;
				ofs << "N=" << (g.nX + 1) * (g.nY + 1) * (g.nZ + 1) << ", E=" << g.nX * g.nY * g.nZ << ", F=FEBLOCK, ET=BRICK";
				ofs << std::endl;
				ofs << "VARLOCATION = (NODAL, NODAL, NODAL";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED";
				ofs << ", CELLCENTERED)";
				ofs << std::endl;

				// Write coordinates
				for (int k = g.dummyCellLayersZ; k <= g.dummyCellLayersZ + g.nZ; k++) {
					for (int j = g.dummyCellLayersY; j <= g.dummyCellLayersY + g.nY; j++) {
						for (int i = g.dummyCellLayersX; i <= g.dummyCellLayersX + g.nX; i++) {
							ofs << g.CoordinateX[i] - 0.5 * g.hx[i];
							ofs << std::endl;
						};
					};
				};
				for (int k = g.dummyCellLayersZ; k <= g.dummyCellLayersZ + g.nZ; k++) {
					for (int j = g.dummyCellLayersY; j <= g.dummyCellLayersY + g.nY; j++) {
						for (int i = g.dummyCellLayersX; i <= g.dummyCellLayersX + g.nX; i++) {
							ofs << g.CoordinateY[j] - 0.5 * g.hy[j];
							ofs << std::endl;
						};
					};
				};
				for (int k = g.dummyCellLayersZ; k <= g.dummyCellLayersZ + g.nZ; k++) {
					for (int j = g.dummyCellLayersY; j <= g.dummyCellLayersY + g.nY; j++) {
						for (int i = g.dummyCellLayersX; i <= g.dummyCellLayersX + g.nX; i++) {
							ofs << g.CoordinateZ[k] - 0.5 * g.hz[k];
							ofs << std::endl;
						};
					};
				};
				ofs.close();
			};
			pManager->Barrier();

			// Compute variables
			auto idx_first = getSerialIndexLocal(g.iMin, g.jMin, g.kMin);
			std::valarray<size_t> s{ size_t(g.nlocalZ), size_t(g.nlocalY), size_t(g.nlocalX) };
			std::valarray<size_t> str{ size_t(g.nlocalYAll * g.nlocalXAll * nVariables), size_t(g.nlocalXAll * nVariables), size_t(nVariables) };

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
			std::valarray<double> P = (gamma - 1.0) * (rho * e);

			// write all data
			pManager->WriteData(rho, fname);
			pManager->WriteData(u, fname);
			pManager->WriteData(v, fname);
			pManager->WriteData(w, fname);
			pManager->WriteData(P, fname);
			pManager->WriteData(e, fname);

			if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);

			// Reopen file for writing
			std::ofstream ofs(fname, std::ios_base::app);

			// write connectivity list
			int fbl, fbr, ful, fur, bbl, bbr, bul, bur;		// indexes for front/back, top/bottom, left/right nodes
			for (int k = g.kMin; k <= g.kMax; k++) {
				for (int j = g.jMin; j <= g.jMax; j++) {
					fbl = (g.nX + 1) * (g.nY + 1) * (k - g.dummyCellLayersZ) + (g.nX + 1) * (j - g.dummyCellLayersY) + (g.iMin - g.dummyCellLayersX) + 1;		// indexes starts from 1 in Tecplot format
					for (int i = g.iMin; i <= g.iMax; i++) {
						fbr = fbl + 1;
						ful = fbl + (g.nX + 1);
						fur = ful + 1;
						bbl = fbl + (g.nX + 1) * (g.nY + 1);
						bbr = bbl + 1;
						bul = bbl + (g.nX + 1);
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

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = static_cast<int>(1 - std::floor(std::log10(SaveSolutionSnapshotTime)));
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
			", iMin = " << g.iMin <<
			", iMax = " << g.iMax
			<< std::endl;
		if (nDims > 1) {
			std::cout <<
				"rankY = " << pManager->rankCart[1] <<
				", jMin = " << g.jMin <<
				", jMax = " << g.jMax
				<< std::endl;
		};
		if (nDims > 2) {
			std::cout <<
				"rankZ = " << pManager->rankCart[2] <<
				", kMin = " << g.kMin <<
				", kMax = " << g.kMax
				<< std::endl;
		};
		if (!pManager->IsLastNode()) pManager->Signal(rank + 1);
		pManager->Barrier();

		for (int iteration = 0; iteration <= MaxIteration; iteration++) {
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
			if ((isSensorEnable == true) && (iteration % SaveSensorRecordIterations == 0))
			{
				for (auto &r : Sensors) {
					r->UpdateIteration(iteration);
					r->UpdateTimer(stepInfo.Time);
					r->Process(values);
				};
			};

			//Solution snapshots
			//Every few iterations
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName << "dataI" << stepInfo.Iteration << ".dat";
				SaveSolution(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout << "Solution has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
				};
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (stepInfo.NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName << std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);
				snapshotFileName << "dataT" << stepInfo.Time << ".dat";
				SaveSolution(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout << "Solution has been written to file \"" << snapshotFileName.str() << "\"" << std::endl;
				};

				//Adjust next snapshot time
				stepInfo.NextSnapshotTime += SaveSolutionSnapshotTime;
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
