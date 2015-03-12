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

#include "KernelConfiguration.h"
#include "ParallelManager.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "RiemannSolvers\RoeSolverPerfectGasEOS.h"
#include "BoundaryConditions\BCGeneral.h"

#include "cgnslib.h"

//Step info
class StepInfo {
public:
	double Time;
	double TimeStep;	
	int Iteration;
	std::vector<double> Residual;
	double NextSnapshotTime;
};

//Calculation kernel
class Kernel {
public:
	//Interface part
	virtual int GetDummyCellLayerSize() = 0; // request as much dummy cell layers as required
	virtual void IterationStep() = 0; // main function

	//Parallel run information
	std::unique_ptr<ParallelManager> pManager;
	int rank; //Rank of process
	int iMin; 
	int iMax;
	int jMin;
	int jMax;
	int kMin;
	int kMax;

	//Current step information
	StepInfo stepInfo;

	//Grid information	
	int nDims; //Number of dimensions
	int nCellsLocalAll; //Total number of local cells including dummy
	int nCellsLocal;
	int nX; //Number of cells in x dimension
	int nY; //Number of cells in y dimension
	int nZ; //Number of cells in z dimension
	int nXAll; //Number of cells in x dimension including dummy layers
	int nYAll; //Number of cells in y dimension including dummy layers
	int nZAll; //Number of cells in z dimension including dummy layers

	//Same for local cells
	int nlocalX;
	int nlocalY;
	int nlocalZ;
	int nlocalXAll;
	int nlocalYAll;
	int nlocalZAll;
	double hx, hy, hz;	//cell sizes
	bool IsPeriodicX; //X periodicity
	bool IsPeriodicY; //Y periodicity
	bool IsPeriodicZ; //Z periodicity
	std::vector<double> CoordinateX; //Cell center coordinates
	std::vector<double> CoordinateY; //Cell center coordinates
	std::vector<double> CoordinateZ; //Cell center coordinates

	//int dummyCellLayers; //number of dummy cell layers
	int dummyCellLayersX;
	int dummyCellLayersY;
	int dummyCellLayersZ;

	//Gas model information
	int nVariables; // number of conservative variables
	double gamma;
	double thermalConductivity;
	double viscosity;
	bool isViscousFlow;
	bool isGradientRequired;

	//External forces
	Vector Sigma; //Potential force	


	//Boundary conditions
	std::unique_ptr<BoundaryConditions::BCGeneral> xLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> xRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zRightBC;

	//Solution data
	std::vector<double> values;	
	std::vector<double> residual;

	//Get serial index for cell
	inline int getSerialIndexGlobal(int i, int j, int k) {
		int sI = (k * nXAll * nYAll + j * nXAll + i);
		return sI;
	};

	inline int getSerialIndexLocal(int i, int j, int k) {
		int sI =  (k - kMin + dummyCellLayersZ) * nlocalXAll * nlocalYAll + (j - jMin + dummyCellLayersY) * nlocalXAll + (i - iMin + dummyCellLayersX);
		return sI;
	};
	
	//Get cell values
	inline double* getCellValues(int i, int j, int k) {
		int sBegin = getSerialIndexLocal(i, j, k) * nVariables;		
		return &values[sBegin];
	};

	//Calculation parameters	
	double MaxIteration;
	double MaxTime;
	double SaveSolutionSnapshotTime;
	int SaveSolutionSnapshotIterations;
	int ResidualOutputIterations;

	bool DebugOutputEnabled; //

	//Constructor
	Kernel(int* argc, char **argv[]) : pManager(new ParallelManager(argc, argv)) { };

	//Set initial conditions
	void SetInitialConditions(std::function<std::vector<double>(Vector r)> initF) {		
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					int sBegin = getSerialIndexLocal(i, j, k) * nVariables;
					std::vector<double> U = initF(Vector(x,y,z));
					for (int varn = 0; varn < nVariables; varn++) values[sBegin + varn] = U[varn];					
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Initial conditions written.\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Initialize kernel
	virtual void Init(KernelConfiguration& config) {
		//Initialize MPI		
		nDims = config.nDims;

		//Initialize local grid
		int dummyCellLayers = GetDummyCellLayerSize(); //number of dummy cell layers				
		nX = config.nX;
		IsPeriodicX = config.isPeriodicX;
		dummyCellLayersX = dummyCellLayers;
		nY = 1;
		IsPeriodicY = true;
		dummyCellLayersY = 0;
		nZ = 1;
		IsPeriodicZ = true;
		dummyCellLayersZ = 0;
		if (nDims > 1) {
			nY = config.nY;
			IsPeriodicY = config.isPeriodicY;
			dummyCellLayersY = dummyCellLayers;
		};
		if (nDims > 2) {
			nZ = config.nZ;
			IsPeriodicZ = config.isPeriodicZ;
			dummyCellLayersZ = dummyCellLayers;
		};

		//Initialize cartezian topology
		StructuredGridInfo gridInfo;
		gridInfo.nDims = nDims;
		gridInfo.nX = nX;
		gridInfo.nY = nY;
		gridInfo.nZ = nZ;
		gridInfo.IsPeriodicX = IsPeriodicX;
		gridInfo.IsPeriodicY = IsPeriodicY;
		gridInfo.IsPeriodicZ = IsPeriodicZ;
		pManager->InitCartesianTopology(gridInfo);
		
		nlocalX = nX / pManager->dimsCart[0];
		nlocalY = nY / pManager->dimsCart[1];
		nlocalZ = nZ / pManager->dimsCart[2];	
		iMin = pManager->rankCart[0] * nlocalX + dummyCellLayersX;
		iMax = (pManager->rankCart[0]+1) * nlocalX + dummyCellLayersX - 1;
		jMin = pManager->rankCart[1] * nlocalY + dummyCellLayersY;
		jMax = (pManager->rankCart[1]+1) * nlocalY + dummyCellLayersY - 1;
		kMin = pManager->rankCart[2] * nlocalZ + dummyCellLayersZ;
		kMax = (pManager->rankCart[2]+1) * nlocalZ + dummyCellLayersZ - 1;		
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->_rankCart<<", iMin = "<<iMin<<", iMax = "<<iMax<<"\n";
			std::cout<<"rank = "<<pManager->_rankCart<<", jMin = "<<jMin<<", jMax = "<<jMax<<"\n";
			std::cout<<"rank = "<<pManager->_rankCart<<", kMin = "<<kMin<<", kMax = "<<kMax<<"\n";
		};

		//Sync
		pManager->Barrier();

		double Lx = config.LX;
		double Ly = config.LY;
		double Lz = config.LZ;	
		if (nDims < 3) Lz = 0;
		if (nDims < 2) Ly = 0;

		//Generate cells (uniform grid)
		nXAll = nX + 2 * dummyCellLayersX;
		nYAll = nY + 2 * dummyCellLayersY;
		nZAll = nZ + 2 * dummyCellLayersZ;
		nlocalXAll = nlocalX + 2 * dummyCellLayersX;
		nlocalYAll = nlocalY + 2 * dummyCellLayersY;
		nlocalZAll = nlocalZ + 2 * dummyCellLayersZ;
		nCellsLocal = nlocalX * nlocalY * nlocalZ;
		nCellsLocalAll = nlocalXAll * nlocalYAll * nlocalZAll;
		CoordinateX.resize(nXAll);
		CoordinateY.resize(nYAll);
		CoordinateZ.resize(nZAll);
		double dx = Lx / nX;
		hx = dx;
		double xMin = 0;
		xMin = xMin - (dummyCellLayersX * dx) + 0.5 * dx;
		double xMax = Lx;
		xMax = xMax + (dummyCellLayersX * dx) - 0.5 * dx;
		//for (int i = iMin - dummyCellLayersX; i < iMax + dummyCellLayersX; i++) {
		for (int i = 0; i < nXAll; i++) {
			double x = xMin + (xMax - xMin) * 1.0 * i / (nXAll - 1);
			CoordinateX[i] = x;
		};
		CoordinateX[iMax + dummyCellLayersX] = xMax;

		double dy = Ly / nY;
		hy = dy;
		double yMin = 0;
		yMin = yMin - (dummyCellLayersY * dy) + 0.5 * dy;
		double yMax = Ly;
		yMax = yMax + (dummyCellLayersY * dy) - 0.5 * dy;
		//for (int j = jMin - dummyCellLayersY; j < jMax + dummyCellLayersY; j++) {
		for (int j = 0; j < nY + dummyCellLayersY; j++) {
			double y = yMin + (yMax - yMin) * 1.0 * j / (nYAll - 1);					
			CoordinateY[j] = y;
		};
		CoordinateY[jMax + dummyCellLayersY] = yMax;

		double dz = Lz / nZ;
		double zMin = 0;
		zMin = zMin - (dummyCellLayersZ * dz) + 0.5 * dz;
		double zMax = Lz;
		zMax = zMax + (dummyCellLayersZ * dz) - 0.5 * dz;
		//for (int k = kMin - dummyCellLayersZ; k < kMax + dummyCellLayersZ; k++) {
		for (int k = 0; k < kMax + dummyCellLayersZ; k++) {
			double z = zMin + (zMax - zMin) * 1.0 * k / (nZAll - 1);
			CoordinateZ[k] = z;
		};
		CoordinateZ[kMax + dummyCellLayersZ] = zMax;

		//Cell linear sizes
		hx = hy = hz = 1.0;
		hx = dx;
		if (nDims > 1) hy = dy;
		if (nDims > 2) hz = dz;

		//Initialize gas model parameters and riemann solver
		gamma = config.Gamma;
		nVariables = config.nVariables;	
		viscosity = config.Viscosity;
		thermalConductivity = config.ThermalConductivity;
		if(config.IsViscousFlow == true) {
			isViscousFlow = true;
			isGradientRequired = true;
		} else {
			isViscousFlow = false;
			isGradientRequired = false;
		};

		//Allocate data structures
		values.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);	
		residual.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);

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
		if (!gridInfo.IsPeriodicX) {
			xLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!gridInfo.IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};
		if ((!gridInfo.IsPeriodicZ) && (nDims > 2)) {
			zLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			zRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			zLeftBC->loadConfiguration(config.zLeftBoundary);
			zRightBC->loadConfiguration(config.zRightBoundary);
		};

		//

		//External forces
		Sigma = Vector(config.Sigma, 0, 0);

		//Output settings
		DebugOutputEnabled = config.DebugOutputEnabled;
		
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Update solution
	void UpdateSolution(double dt) {
		//Compute cell volume
		double volume = hx * hy * hz;
		stepInfo.Residual.resize(nVariables);
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = 0;


		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);

					//Update cell values
					for(int nv = 0; nv < nVariables; nv++) values[idx * nVariables + nv] += residual[idx * nVariables + nv] * dt / volume;

					//Compute total residual
					stepInfo.Residual[0] += abs(residual[idx * nVariables]);			//ro
					stepInfo.Residual[1] += abs(residual[idx * nVariables + 1]);		//rou
					stepInfo.Residual[2] += abs(residual[idx * nVariables + 2]);		//rov
					stepInfo.Residual[3] += abs(residual[idx * nVariables + 3]);		//row
					stepInfo.Residual[4] += abs(residual[idx * nVariables + 4]);		//roE
				};
			};
		};

		//Aggregate
		stepInfo.Residual[0] = pManager->Sum(stepInfo.Residual[0]);
		stepInfo.Residual[1] = pManager->Sum(stepInfo.Residual[1]);
		stepInfo.Residual[2] = pManager->Sum(stepInfo.Residual[2]);
		stepInfo.Residual[3] = pManager->Sum(stepInfo.Residual[3]);
		stepInfo.Residual[4] = pManager->Sum(stepInfo.Residual[4]);
	};

	//Exchange values between processormos
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
		int rankLeft; // if equals -1 no left neighbour
		int rankRight; // if equals -1 no right neighbour
		
		//X direction exchange		

		//Allocate buffers
		int layerSize = nlocalY * nlocalZ;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);				

		//Determine neighbours' ranks
		int rankL = pManager->GetRankByCartesianIndexShift(-1, 0, 0);
		int rankR = pManager->GetRankByCartesianIndexShift(+1, 0, 0);

		for (int layer = 1; layer <= dummyCellLayersX; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int iSend = iMin + layer - 1; // layer index to send
			int iRecv = iMax + layer; // layer index to recv

			// Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (j = jMin; j <= jMax; j++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
						int idxValues = getSerialIndexLocal(iSend, j, k);
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
			for (j = jMin; j <= jMax; j++) {
				for (k = kMin; k <= kMax; k++) {
					int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
					int idxValues = getSerialIndexLocal(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			iSend = iMax - layer + 1; // layer index to send
			iRecv = iMin - layer; // layer index to recv

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (j = jMin; j <= jMax; j++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
						int idxValues = getSerialIndexLocal(iSend, j, k);
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
			for (j = jMin; j <= jMax; j++) {
				for (k = kMin; k <= kMax; k++) {
					int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
					int idxValues = getSerialIndexLocal(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values X-direction executed\n";
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction exchange		

		//Allocate buffers
		layerSize = nlocalX * nlocalZ;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);				

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);

		for (int layer = 1; layer <= dummyCellLayersY; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int jSend = jMin + layer - 1; // layer index to send
			int jRecv = jMax + layer; // layer index to recv

			// Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = iMin; i <= iMax; i++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (i - iMin) + (k-kMin) * nX; //Exclude y index
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
			for (i = iMin; i <= iMax; i++) {
				for (k = kMin; k <= kMax; k++) {
					int idxBuffer = (i - iMin) + (k-kMin) * nX; //Exclude y index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			jSend = jMax - layer + 1; // layer index to send
			jRecv = jMin - layer; // layer index to recv

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = iMin; i <= iMax; i++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (i - iMin) + (k-kMin) * nX; //Exclude y index
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
			for (i = iMin; i <= iMax; i++) {
				for (k = kMin; k <= kMax; k++) {
					int idxBuffer = (i - iMin) + (k-kMin) * nX; //Exclude y index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values Y-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction exchange		

		//Allocate buffers
		layerSize = nlocalX * nlocalY;
		bufferToSend.resize(layerSize * nVariables);
		bufferToRecv.resize(layerSize * nVariables);				

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, 0, -1);
		rankR = pManager->GetRankByCartesianIndexShift(0, 0, +1);

		for (int layer = 1; layer <= dummyCellLayersZ; layer++) {
			// Minus direction exchange
			int nSend = 0; //
			int nRecv = 0; //
			int kSend = kMin + layer - 1; // layer index to send
			int kRecv = kMax + layer; // layer index to recv

			// Prepare values to send
			if (rankL != -1) {
				nSend = layerSize * nVariables;
				for (i = iMin; i <= iMax; i++) {
					for (j = jMin; j <= jMax; j++) {
						int idxBuffer = (i - iMin) + (j - jMin) * nX; //Exclude z index
						int idxValues = getSerialIndexLocal(i, j, kSend);
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
			for (i = iMin; i <= iMax; i++) {
				for (j = jMin; j <= jMax; j++) {
					int idxBuffer = (i - iMin) + (j - jMin) * nX; //Exclude z index
					int idxValues = getSerialIndexLocal(i, j, kRecv);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			kSend = kMax - layer + 1; // layer index to send
			kRecv = kMin - layer; // layer index to recv

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = iMin; i <= iMax; i++) {
					for (j = jMin; j <= jMax; j++) {
						int idxBuffer = (i - iMin) + (j - jMin) * nX; //Exclude z index
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
			for (i = iMin; i <= iMax; i++) {
				for (j = jMin; j <= jMax; j++) {
					int idxBuffer = (i - iMin) + (j - jMin) * nX; //Exclude z index
					int idxValues = getSerialIndexLocal(i, j, kRecv);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values Z-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

	}; // function

	//Compute dummy cell values as result of boundary conditions and interprocessor exchange communication
	void ComputeDummyCellValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//Interprocessor exchange
		ExchangeValues();

		//Current face and cell information
		Vector faceNormalL;
		Vector faceNormalR;
		Vector faceCenter;
		Vector cellCenter;
		
		//X direction		
		faceNormalL = Vector(-1.0, 0.0, 0.0);
		faceNormalR = Vector(1.0, 0.0, 0.0);
		if (!IsPeriodicX) {
			for (j = jMin; j <= jMax; j++) {
				for (k = kMin; k <= kMax; k++) {				
					//Inner cell
					cellCenter.y = CoordinateY[j];
					cellCenter.z = CoordinateZ[k];

					//And face
					faceCenter.y = CoordinateY[j];
					faceCenter.z = CoordinateZ[k];

					for (int layer = 1; layer <= dummyCellLayersX; layer++) {		
						if (pManager->rankCart[0] == 0) {
							//Left border
							i = iMin - layer; // layer index
							int iIn = iMin + layer - 1; // opposite index
							cellCenter.x = CoordinateX[iIn];
							faceCenter.x = (CoordinateX[iIn] + CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);
					
							//Apply left boundary conditions						
							std::vector<double> dValues = xLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[0] == pManager->dimsCart[0]-1) {
							//Right border
							i = iMax + layer; // layer index
							int iIn = iMax - layer + 1; // opposite index
							cellCenter.x = CoordinateX[iIn];
							faceCenter.x = (CoordinateX[iIn] + CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);
					
							//Apply right boundary conditions						
							std::vector<double> dValues = xRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};		
						};
					};
				};
			};
		}; //X direction

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Dummy values X-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction		
		faceNormalL = Vector(0.0, -1.0, 0.0);
		faceNormalR = Vector(0.0, 1.0, 0.0);
		for (i = iMin; i <= iMax; i++) {
			for (k = kMin; k <= kMax; k++) {				
				//Inner cell
				cellCenter.x = CoordinateX[i];
				cellCenter.z = CoordinateZ[k];

				//And face
				faceCenter.x = CoordinateX[i];
				faceCenter.z = CoordinateZ[k];

				for (int layer = 1; layer <= dummyCellLayersY; layer++) {
					if (!IsPeriodicY) {
						if (pManager->rankCart[1] == pManager->dimsCart[1]-1) {
							//Left border
							j = jMin - layer; // layer index
							int jIn = jMin + layer - 1; // opposite index
							cellCenter.y = CoordinateY[jIn];
							faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);
					
							//Apply left boundary conditions						
							std::vector<double> dValues = yLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[1] == pManager->dimsCart[1]-1) {
							//Right border
							j = jMax + layer; // layer index
							int jIn = jMax - layer + 1; // opposite index
							cellCenter.y = CoordinateY[jIn];
							faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);
					
							//Apply right boundary conditions						
							std::vector<double> dValues = yRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};
					};
				};
			};
		}; //Y direction

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Dummy values Y-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Save solution CGNS
	void SaveSolutionPCGNS(std::string fname) {
		//int comm_size, comm_rank;
		//int tot_nnodes, tot_nelems, nnodes, nelems;
		//int F, B, Z, E, S, Fs, A, Cx, Cy, Cz;
		//int i, j, k, n, nn, ne;
		//float *x, *y, *z, *d;
		//cgsize_t sizes[3], *e, start, end, ncells;
 
		///* open the file and create base and zone */
		//sizes[0] = tot_nnodes;
		//sizes[1] = tot_nelems;
		//sizes[2] = 0;

		///* the default here is to use MPI_COMM_WORLD,
		//   but this allows assigning of another communicator
		//cgp_mpi_comm(MPI_COMM_WORLD); */

		//if (cg_open(fname.c_str(), CG_MODE_WRITE, &F) ||
		//	cg_base_write(F, "Base", nDims, nDims, &B) ||
		//	cg_zone_write(F, B, "Zone", sizes, Unstructured, &Z))
		//	cg_error_exit();

		///* print info */
		//if (comm_rank == 0) {
		//	printf("writing %d coordinates and %d hex elements to %s\n",
		//		tot_nnodes, tot_nelems, outfile);
		//}

		///* create data nodes for coordinates */
		//if (cgp_coord_write(F, B, Z, RealSingle, "CoordinateX", &Cx) ||
		//	cgp_coord_write(F, B, Z, RealSingle, "CoordinateY", &Cy) ||
		//	cgp_coord_write(F, B, Z, RealSingle, "CoordinateZ", &Cz))
		//	cgp_error_exit();
 
		///* number of nodes and range this process will write */
		//nnodes = (tot_nnodes + comm_size - 1) / comm_size;
		//start  = nnodes * comm_rank + 1;
		//end    = nnodes * (comm_rank + 1);
		//if (end > tot_nnodes) end = tot_nnodes;
  //  
		///* create the coordinate data for this process */
		//x = (float *)malloc(nnodes * sizeof(float));
		//y = (float *)malloc(nnodes * sizeof(float));
		//z = (float *)malloc(nnodes * sizeof(float));
		//nn = 0;
		//for (n = 1, k = 0; k < NODES_PER_SIDE; k++) {
		//	for (j = 0; j < NODES_PER_SIDE; j++) {
		//		for (i = 0; i < NODES_PER_SIDE; i++, n++) {
		//			if (n >= start && n <= end) {
		//				x[nn] = (float)i;
		//				y[nn] = (float)j;
		//				z[nn] = (float)k;
		//				nn++;
		//			}
		//		}
		//	}
		//}

		///* write the coordinate data in parallel to the queue */
		//if (cgp_queue_set(1) ||
		//	cgp_coord_write_data(F, B, Z, Cx, &start, &end, x) ||
		//	cgp_coord_write_data(F, B, Z, Cy, &start, &end, y) ||
		//	cgp_coord_write_data(F, B, Z, Cz, &start, &end, z))
		//	cgp_error_exit();
  //  
		///* write out the queued coordinate data */
		//if (cgp_queue_flush()) cgp_error_exit();
		//cgp_queue_set(0);

		///* create data node for elements */
		//if (cgp_section_write(F, B, Z, "Hex", HEXA_8, 1, tot_nelems, 0, &E))
		//	cgp_error_exit();
 
		///* number of elements and range this process will write */
		//nelems = (tot_nelems + comm_size - 1) / comm_size;
		//start  = nelems * comm_rank + 1;
		//end    = nelems * (comm_rank + 1);
		//if (end > tot_nelems) end = tot_nelems;
  //  
		///* create the hex element data for this process */
		//e = (cgsize_t *)malloc(8 * nelems * sizeof(cgsize_t));
		//nn = 0;
		//for (n = 1, k = 1; k < NODES_PER_SIDE; k++) {
		//	for (j = 1; j < NODES_PER_SIDE; j++) {
		//		for (i = 1; i < NODES_PER_SIDE; i++, n++) {
		//			if (n >= start && n <= end) {
		//				ne = i + NODES_PER_SIDE*((j-1)+NODES_PER_SIDE*(k-1));
		//				e[nn++] = ne;
		//				e[nn++] = ne + 1;
		//				e[nn++] = ne + 1 + NODES_PER_SIDE;
		//				e[nn++] = ne + NODES_PER_SIDE;
		//				ne += NODES_PER_SIDE * NODES_PER_SIDE;
		//				e[nn++] = ne;
		//				e[nn++] = ne + 1;
		//				e[nn++] = ne + 1 + NODES_PER_SIDE;
		//				e[nn++] = ne + NODES_PER_SIDE;
		//			}
		//		}
		//	}
		//}

		///* write the element connectivity in parallel */
		//if (cgp_elements_write_data(F, B, Z, E, start, end, e))
		//	cgp_error_exit();

		///* create a centered solution */
		//if (cg_sol_write(F, B, Z, "Solution", CellCenter, &S) ||
		//	cgp_field_write(F, B, Z, S, RealSingle, "CellIndex", &Fs))
		//	cgp_error_exit();

		///* create the field data for this process */
		//d = (float *)malloc(nelems * sizeof(float));
		//nn = 0;
		//for (n = 1; n <= tot_nelems; n++) {
		//	if (n >= start && n <= end) {
		//		d[nn] = (float)n;
		//		nn++;
		//	}
		//}

		///* write the solution field data in parallel */
		//if (cgp_field_write_data(F, B, Z, S, Fs, &start, &end, d))
		//	cgp_error_exit();

		///* create user data under the zone and duplicate solution data */
		//ncells = tot_nelems;
		//if (cg_goto(F, B, "Zone_t", 1, NULL) ||
		//	cg_user_data_write("User Data") ||
		//	cg_gorel(F, "User Data", 0, NULL) ||
		//	cg_array_write("CellIndex", RealSingle, 1, &ncells, &A))
		//	cg_error_exit();

		///* write the array data in parallel */
		//if (cg_array_write_data(A, &start, &end, d))
		//	cg_error_exit();

		///* close the file and terminate MPI */
		//cg_close(F);   
		//return;

	};

	void SaveSolutionStructuredCGNS(std::string fname) {
		int flag = 0;
		int rank = pManager->getRank();
		int np = pManager->getProcessorNumber();
		MPI_Status status;

		//Function to determine value position in buffer
		std::vector<int> index(3, 0);
		std::vector<double> outBuffer(nCellsLocal);
		auto getIndex = [] (std::vector<int> index, int nDims, int dims[]) {
			int ind = index[0];
			if (nDims > 1) ind += (dims[0]) * index[1];
			if (nDims > 2) ind += (dims[1]) * (dims[0]) * index[2];
		/*	int ind = index[0];
			if (nDims > 1) ind = ind * dims[1] + index[1];
			if (nDims > 2) ind = ind * dims[2] + index[2];*/
			return ind;
		};

		//Indexes inside CGNS database
		int F, B, Z, S, C, Fs;
 
		//Create file
		if (pManager->IsMasterCart()) {
			//Delete file if exists
			if (remove(fname.c_str()))
				std::cout << "Failed to delete '" << fname << "': " << strerror(errno) << '\n';
			//else
				//std::cout << '\'' << fname << "' successfully deleted.\n";

			/* Open the file and create base and zone */
			std::vector<cgsize_t> sizes;

			//Number of vertexes in each direction
			int nVertexes = 1;
			nVertexes *= nX + 1;
			sizes.push_back(nX + 1);
			if (nDims > 1) {
				nVertexes *= nY + 1;
				sizes.push_back(nY + 1);
			};
			if (nDims > 2) {
				nVertexes *= nZ + 1;
				sizes.push_back(nZ + 1);
			};

			//Number of cells in each direction
			sizes.push_back(nX);
			if (nDims > 1) sizes.push_back(nY);
			if (nDims > 2) sizes.push_back(nZ);


			int cellDim = nDims;
			int physDim = nDims;
			if (cg_open(fname.c_str(), CG_MODE_WRITE, &F) ||
				cg_base_write(F, "Base", cellDim, physDim, &B) ||
				cg_zone_write(F, B, "Zone", &sizes.front(), Structured, &Z))
				cg_error_exit();

			//Write solution node
			if (cg_sol_write(F, B, Z, "Solution", CellCenter, &S))
				cg_error_exit();

			//Close file
			cg_close(F);
		};

		//Wait for file structure to be written before writing parts of solution on each process
		pManager->Barrier();

		std::cout<<"rank = "<<rank<<", file has been created.\n";
		std::cout.flush();		

		if (!pManager->IsMaster()) {
			//And have to go to same zone node
			B = 1;
			Z = 1;
			S = 1;
		};

		/* create data nodes for coordinates */

		

		int dimsC[3];
		dimsC[0] = nlocalX;
		dimsC[1] = nlocalY;
		dimsC[2] = nlocalZ;
		auto getCellIndex = std::bind(getIndex, std::placeholders::_1, nDims, dimsC);

		//Open file
		if (cg_open(fname.c_str(), CG_MODE_MODIFY, &F)) 
			cg_error_exit();


		/*} catch (std::exception exc) {
			std::cout<<"rank = "<<rank<<", exception msg : "<<exc.what()<<"\n";
			std::cout.flush();
		};*/
		//};

		
		std::cout<<"rank = "<<rank<<", np = "<<np<<"\n";
		std::cout.flush();		

		int currentRank = 0;
		for (currentRank = 0; currentRank < 1; currentRank++) {

			if (currentRank != rank) {
		
	//			pManager->Barrier();
		//		continue;
			};

			//Begin of critical section
			//Create vertex coordinate arrays
			std::vector<cgsize_t> startIndex(0);
			startIndex.push_back(iMin - dummyCellLayersX + 1);
			startIndex.push_back((nDims > 1) ? jMin - dummyCellLayersY + 1 : 1);
			startIndex.push_back((nDims > 2) ? kMin - dummyCellLayersZ + 1 : 1);
			std::vector<cgsize_t> endIndex(0);
			endIndex.push_back(iMax - dummyCellLayersX + 1);
			endIndex.push_back((nDims > 1) ? jMax - dummyCellLayersY + 1 : 1);
			endIndex.push_back((nDims > 2) ? kMax - dummyCellLayersZ + 1 : 1);
			int nVertexesLocal = 1;
			int dimsV[3];
			for (int dn = 0; dn < nDims; dn++) {
				if (pManager->rankCart[dn] == pManager->dimsCart[dn] - 1) endIndex[dn]++;
				dimsV[dn] = endIndex[dn] - startIndex[dn] + 1;
				nVertexesLocal *= dimsV[dn];
			};

			auto getVertexIndex = std::bind(getIndex, std::placeholders::_1, nDims, dimsV);

			std::cout<<"rank = "<<rank<<", buffer size = "<<outBuffer.size()<<"\n";
			std::cout<<"iMin = "<<startIndex[0]<<", iMax = "<<endIndex[0]<<"\n";
			std::cout<<"jMin = "<<startIndex[1]<<", jMax = "<<endIndex[1]<<"\n";
			std::cout<<"kMin = "<<startIndex[2]<<", kMax = "<<endIndex[2]<<"\n";
			std::cout.flush();

			//Write x-coords
			std::vector<double> outBufferX;
			outBufferX.resize(nVertexesLocal);
			for (index[0] = 0; index[0] <= (endIndex[0] - startIndex[0]); index[0]++) {
				for (index[1] = 0; index[1] <= (endIndex[1] - startIndex[1]); index[1]++) {
					for (index[2] = 0; index[2] <= (endIndex[2] - startIndex[2]); index[2]++) {
						int i = index[0] + iMin;
						int ind = getVertexIndex(index);
						outBufferX[ind] = (CoordinateX[i - 1] + CoordinateX[i]) / 2.0;
					};
				};
			};
			if (cg_coord_partial_write(F, B, Z, RealDouble, "CoordinateX", &startIndex.front(), &endIndex.front(), &outBufferX.front(), &C))
				cg_error_exit();

			pManager->Barrier();

			if (nDims > 1) {
				std::vector<double> outBufferY;
				outBufferY.resize(nVertexesLocal);

				//Write y-coords
				for (index[0] = 0; index[0] <= (endIndex[0] - startIndex[0]); index[0]++) {
					for (index[1] = 0; index[1] <= (endIndex[1] - startIndex[1]); index[1]++) {
						for (index[2] = 0; index[2] <= (endIndex[2] - startIndex[2]); index[2]++) {
							int j = index[1] + jMin;
							int ind = getVertexIndex(index);
							outBufferY[ind] = (CoordinateY[j - 1] + CoordinateY[j]) / 2.0;
						};
					};
				};
				if (cg_coord_partial_write(F, B, Z, RealDouble, "CoordinateY", &startIndex.front(), &endIndex.front(), &outBufferY.front(), &C))
					cg_error_exit();
			};

			if (nDims > 2) {
				std::vector<double> outBufferZ;
				outBufferZ.resize(nVertexesLocal);

				//Write z-coords
				for (index[0] = 0; index[0] <= (endIndex[0] - startIndex[0]); index[0]++) {
					for (index[1] = 0; index[1] <= (endIndex[1] - startIndex[1]); index[1]++) {
						for (index[2] = 0; index[2] <= (endIndex[2] - startIndex[2]); index[2]++) {
							int k = index[2] + kMin;
							int ind = getVertexIndex(index);
							outBufferZ[ind] = (CoordinateZ[k - 1] + CoordinateZ[k]) / 2.0;
						};
					};
				};
				if (cg_coord_partial_write(F, B, Z, RealDouble, "CoordinateZ", &startIndex.front(), &endIndex.front(), &outBufferZ.front(), &C))
					cg_error_exit();
			};


			std::cout<<"rank = "<<rank<<", coordinates written.\n";
			std::cout.flush();	


			std::cout<<"rank = "<<rank<<", Solution writing begin\n";
			std::cout.flush();

			//Write solution data part for each process
			startIndex.clear();
			startIndex.push_back(iMin - dummyCellLayersX + 1);
			startIndex.push_back((nDims > 1) ? jMin - dummyCellLayersY + 1 : 1);
			startIndex.push_back((nDims > 2) ? kMin - dummyCellLayersZ + 1 : 1);
			endIndex.clear();
			endIndex.push_back(iMax - dummyCellLayersX + 1);
			endIndex.push_back((nDims > 1) ? jMax - dummyCellLayersY + 1 : 1);
			endIndex.push_back((nDims > 2) ? kMax - dummyCellLayersZ + 1 : 1);

			std::cout<<"rank = "<<rank<<", buffer size = "<<outBuffer.size()<<"\n";
			std::cout<<"iMin = "<<startIndex[0]<<", iMax = "<<endIndex[0]<<"\n";
			std::cout<<"jMin = "<<startIndex[1]<<", jMax = "<<endIndex[1]<<"\n";
			std::cout<<"kMin = "<<startIndex[2]<<", kMax = "<<endIndex[2]<<"\n";
			std::cout.flush();
		
			//Function that returns density
			auto getDensity = [] (double *U) {
				return U[0];
			};


			//Writing arbitrary function applied to solution for each cell into out buffer
			//auto fillBufferWith = [startIndex, endIndex, getIndex] (std::vector<double> outBuffer) {
			for (index[0] = 0; index[0] <= (endIndex[0] - startIndex[0]); index[0]++) {
				for (index[1] = 0; index[1] <= (endIndex[1] - startIndex[1]); index[1]++) {
					for (index[2] = 0; index[2] <= (endIndex[2] - startIndex[2]); index[2]++) {
						int bIndex = getCellIndex(index);
						int i = index[0] + iMin;
						int j = index[1] + jMin;
						int k = index[2] + kMin;
						double *U = getCellValues(i, j, k);
						double val = getDensity(U);
						//std::cout<<"rank = "<<rank<<", bIndex : "<<bIndex<<"\n";
						//std::cout.flush();
						outBuffer[bIndex] = val;
					};
				};
			};

			//Density
			if (cg_field_partial_write(F, B, Z, S, RealDouble, "Density", &startIndex.front(), &endIndex.front(), &outBuffer.front(), &Fs))
				cg_error_exit(); 
		
			std::cout<<"rank = "<<rank<<", Density field written\n";
			std::cout.flush();

			/* close the file */
			cg_close(F); 

			//End of critical section
			pManager->Barrier();
		};


		//Final sync
		pManager->Barrier();

		//Read coordinates
		/*if (rank == 0) {
			if (cg_open(fname.c_str(), CG_MODE_READ, &F)) 
				cg_error_exit();

			std::vector<cgsize_t> startIndex(0);
			startIndex.push_back(iMin - dummyCellLayersX + 1);
			startIndex.push_back((nDims > 1) ? jMin - dummyCellLayersY + 1 : 1);
			startIndex.push_back((nDims > 2) ? kMin - dummyCellLayersZ + 1 : 1);
			std::vector<cgsize_t> endIndex(0);
			endIndex.push_back(iMax - dummyCellLayersX + 1);
			endIndex.push_back((nDims > 1) ? jMax - dummyCellLayersY + 1 : 1);
			endIndex.push_back((nDims > 2) ? kMax - dummyCellLayersZ + 1 : 1);
			int nVertexesLocal = 1;
			int dimsV[3];
			for (int dn = 0; dn < nDims; dn++) {
				if (pManager->rankCart[dn] == pManager->dimsCart[dn] - 1) endIndex[dn]++;
				dimsV[dn] = endIndex[dn] - startIndex[dn] + 1;
				nVertexesLocal *= dimsV[dn];
			};
			outBuffer.resize(nVertexesLocal);

			cg_coord_read(F, B, Z, "CoordinateX", RealDouble, &startIndex.front(), &endIndex.front(), &outBuffer.front());

			std::ofstream ofs("CoordX.txt");
			for (int i = 0; i<outBuffer.size(); i++) ofs<<outBuffer[i]<<std::endl;
			ofs.close();


			cg_coord_read(F, B, Z, "CoordinateY", RealDouble, &startIndex.front(), &endIndex.front(), &outBuffer.front());

			std::ofstream ofs2("CoordY.txt");
			for (int i = 0; i<outBuffer.size(); i++) ofs2<<outBuffer[i]<<std::endl;
			ofs2.close();

			cg_close(F);
			exit(0);
		};*/

		return;
	};

	//Save solution
	void SaveSolution(std::string fname) {
		//Tecpol version
		std::ofstream ofs;

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->_rankCart<<", Solution save begin\n";
			std::cout.flush();
		};
		pManager->Barrier();
		
		//Header
		if (pManager->IsMasterCart()) {
			ofs.open(fname, std::ios_base::out);
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"X"<<"\" ";
			if (nDims > 1) ofs<<"\""<<"Y"<<"\" ";
			if (nDims > 2) ofs<<"\""<<"Z"<<"\" ";
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			if (nDims > 1) ofs<<"\""<<"v"<<"\" ";
			if (nDims > 2) ofs<<"\""<<"w"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"e"<<"\" ";
			ofs<<std::endl;

			std::string zoneTitle = "Time = ";
			ofs<<"ZONE T = \""<<zoneTitle<<stepInfo.Time<<"\", ";

			ofs<<"I = "<<nX<<", ";
			if (nDims > 1) ofs<<"J = "<<nY<<", ";
			if (nDims > 2) ofs<<"K = "<<nZ<<", ";
			ofs<<"DATAPACKING=POINT"<<std::endl;
		};

		//Obtain indexes for each process
		std::vector<int> iMinAll(pManager->getProcessorNumber());
		pManager->GatherCounts(iMin, iMinAll);
		std::vector<int> iMaxAll(pManager->getProcessorNumber());
		pManager->GatherCounts(iMax, iMaxAll);
		std::vector<int> jMinAll(pManager->getProcessorNumber(), 0);
		std::vector<int> jMaxAll(pManager->getProcessorNumber(), 0);
		if (nDims > 1) {
			pManager->GatherCounts(jMin, jMinAll);
			pManager->GatherCounts(jMax, jMaxAll);
		};
		std::vector<int> kMinAll(pManager->getProcessorNumber(), 0);		
		std::vector<int> kMaxAll(pManager->getProcessorNumber(), 0);		
		if (nDims > 2) {
			pManager->GatherCounts(kMin, kMinAll);
			pManager->GatherCounts(kMax, kMaxAll);
		};

		//Send solution to master from slave nodes
		int masterRank = 0;
		if (!pManager->IsMasterCart()) {	
			MPI_Send(&values.front(), values.size(), MPI_LONG_DOUBLE, masterRank, 0, pManager->_commCart);
		};		
	

		//On master output each part
		if (pManager->IsMasterCart()) {
			std::cout<<"rank = "<<pManager->_rankCart<<", Solution writing begin\n";
			std::cout.flush();			

			//Output solution
			int masterRank = 0;
			std::vector<double> recvBuff;			
			for (int rI = 0; rI < pManager->dimsCart[0]; rI++) {
				for (int rJ = 0; rJ < pManager->dimsCart[1]; rJ++) {
					for (int rK = 0; rK < pManager->dimsCart[2]; rK++) {
						int rank = pManager->GetRankByCartesianIndex(rI, rJ, rK);
						std::cout<<"rank = "<<rank<<", Part of solution writing begin\n";
						std::cout.flush();

						//Change view to current process
						iMin = iMinAll[rank];
						iMax = iMaxAll[rank];
						if (nDims > 1) {
							jMin = jMinAll[rank];
							jMax = jMaxAll[rank];
						};
						if (nDims > 2) {
							kMin = kMinAll[rank];
							kMax = kMaxAll[rank];
						};

						//Determine size of buffer
						int size = (iMax - iMin + 2 * dummyCellLayersX + 1);
						if (nDims > 1) size *= (jMax - jMin + 2 * dummyCellLayersY + 1);
						if (nDims > 2) size *= (kMax - kMin + 2 * dummyCellLayersZ + 1);

						//Get values into buffer
						if (rank == masterRank) {
							//Do not recieve
							recvBuff = std::vector<double>(values.begin(), values.end());								
						} else {
							//Recieve
							recvBuff.resize(size * nVariables);								
							MPI_Status status;
							MPI_Recv(&recvBuff.front(), recvBuff.size(), MPI_LONG_DOUBLE, rank, 0, pManager->_commCart, &status);
						};
						
						//Solution output
						for (int i = iMin; i <= iMax; i++) {
							for (int j = jMin; j <= jMax; j++) {
								for (int k = kMin; k <= kMax; k++) {
									//Obtain cell data
									double x = CoordinateX[i];
									double y = CoordinateY[j];
									double z = CoordinateZ[k];
									int sidx = getSerialIndexLocal(i, j, k) * nVariables;
									double* U = &recvBuff[sidx];
									double ro = U[0];
									double u = U[1] / ro;
									double v = U[2] / ro;
									double w = U[3] / ro;
									double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
									double P = (gamma - 1.0) * ro * e;

									//Write to file
									ofs<<x<<" ";
									if (nDims > 1) ofs<<y<<" ";
									if (nDims > 2) ofs<<z<<" ";
									ofs<<ro<<" ";
									ofs<<u<<" ";
									if (nDims > 1) ofs<<v<<" ";
									if (nDims > 2) ofs<<w<<" ";
									ofs<<P<<" ";
									ofs<<e<<" ";
									ofs<<std::endl;
								};
							};
						}; // Solution output

						std::cout<<"rank = "<<rank<<", Part of solution written\n";
						std::cout.flush();

					}; // for rK
				}; // for rJ
			}; // for rI

			//Restore indexes
			iMin = iMinAll[masterRank];
			iMax = iMaxAll[masterRank];
			if (nDims > 1) {
				jMin = jMinAll[masterRank];
				jMax = jMaxAll[masterRank];
			};
			if (nDims > 2) {			
				kMin = kMinAll[masterRank];
				kMax = kMaxAll[masterRank];
			};

		};// if	

		//Close file
		if (pManager->IsMasterCart()) ofs.close();

		std::cout<<"rank = "<<pManager->getRank()<<", Solution written\n";
		std::cout.flush();
		//Sync
		pManager->Barrier();					

		return;

		//2D tecplot style
		if(nDims > 1) {			
			if (pManager->IsMaster()) {
				ofs.open(fname, std::ios_base::out);
				//Header				
				ofs<<"VARIABLES = ";
				ofs<<"\""<<"X"<<"\" ";
				ofs<<"\""<<"Y"<<"\" ";
				ofs<<"\""<<"Z"<<"\" ";
				ofs<<"\""<<"ro"<<"\" ";
				ofs<<"\""<<"u"<<"\" ";
				ofs<<"\""<<"v"<<"\" ";
				ofs<<"\""<<"w"<<"\" ";
				ofs<<"\""<<"P"<<"\" ";
				ofs<<"\""<<"e"<<"\" ";
				ofs<<std::endl;
			
				ofs << "ZONE T=\"1\"\nI=" << nX <<"J=" << nY << "K=" << nZ << "F=POINT\n";
				ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
			};

			////Solution
			//for (int k = kMin; k <= kMax; k++) {
			//	for (int j = jMin; j <= jMax; j++) {
			//		for (int i = iMin; i <= iMax; i++) {
			//			//Obtain cell data
			//			double x = CoordinateX[i];
			//			double y = CoordinateY[j];
			//			double z = CoordinateZ[k];
			//			double* U = getCellValues(i,j,k);
			//			double ro = U[0];
			//			double u = U[1] / ro;
			//			double v = U[2] / ro;
			//			double w = U[3] / ro;
			//			double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
			//			double P = (gamma - 1.0) * ro * e;

			//			//Write to file
			//			ofs<<x<<" ";
			//			ofs<<y<<" ";
			//			ofs<<z<<" ";
			//			ofs<<ro<<" ";
			//			ofs<<u<<" ";
			//			ofs<<v<<" ";
			//			ofs<<w<<" ";
			//			ofs<<P<<" ";
			//			ofs<<e<<" ";
			//			ofs<<std::endl;
			//		};
			//	};
			//};
		};

		if (pManager->IsMaster()) ofs.close();
		
		std::cout<<"rank = "<<pManager->getRank()<<", Solution written\n";
		std::cout.flush();
		//Sync
		pManager->Barrier();
	};

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = 1 - std::floor(std::log10(SaveSolutionSnapshotTime));
		};

		//Open history file
		std::ofstream ofs;
		if (pManager->IsMaster()) {
			ofs.open("history.dat", std::ios_base::out);
			//Header
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"time"<<"\" ";			
			ofs<<"\""<<"roR"<<"\" ";
			ofs<<"\""<<"rouR"<<"\" ";
			if (nDims > 1) ofs<<"\""<<"rovR"<<"\" ";
			if (nDims > 2) ofs<<"\""<<"rowR"<<"\" ";			
			ofs<<"\""<<"roeR"<<"\" ";
			ofs<<"\""<<"uMax"<<"\" ";
			ofs<<std::endl;			
		};

		//Start timer
		Timer workTimer;
		workTimer.Start();
		workTimer.Resume();

		//Calc loop		
		if (pManager->IsMaster()) {
			std::cout<<"Calculation started!\n";
		};

		for (int iteration = 0; iteration <= MaxIteration; iteration++) {
			//Exchange and boundary conditions
			ComputeDummyCellValues();

			//Calculate one time step
			IterationStep();

			//Output step information						
			if (pManager->IsMaster() && (ResidualOutputIterations != 0) && (stepInfo.Iteration % ResidualOutputIterations) == 0) {
				std::cout<<"Iteration = "<<stepInfo.Iteration<<"; Total time = "<< stepInfo.Time << "; Time step = " <<stepInfo.TimeStep << "; RMSrou = "<<stepInfo.Residual[1]<<"\n";			
			};			

			//Solution snapshots
			//Every few iterations
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".dat";						
				SaveSolution(snapshotFileName.str());
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (stepInfo.NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);								
				snapshotFileName<<"dataT"<<stepInfo.Time<<".dat";							
				SaveSolution(snapshotFileName.str());

				//Adjust next snapshot time
				stepInfo.NextSnapshotTime += SaveSolutionSnapshotTime;
			};

			//Find maximum x velocity	
			double uMax = 0;
			for (int i = iMin; i <= iMax; i++) {
				for (int j = jMin; j <= jMax; j++) {
					for (int k = kMin; k <= kMax; k++) {
						double *U = getCellValues(i, j, k);
						double u = U[1] / U[0];
						if (std::abs(u) > std::abs(uMax)) uMax = u;
					};
				};
			};

			double uMaxAbs = std::abs(uMax);
			uMaxAbs = pManager->Max(uMaxAbs);
			if (uMax > 0) {
				uMax = uMaxAbs;
			} else {
				uMax = -uMaxAbs;
			};
			
			//Save convergence history		
			if (pManager->IsMaster() && (stepInfo.Iteration % ResidualOutputIterations)) {
				ofs<<stepInfo.Time<<" ";			
				ofs<<stepInfo.Residual[0]<<" ";							
				ofs<<stepInfo.Residual[1]<<" ";			
				if (nDims > 1) ofs<<stepInfo.Residual[2]<<" ";			
				if (nDims > 2) ofs<<stepInfo.Residual[3]<<" ";			
				ofs<<stepInfo.Residual[4]<<" ";				
				ofs<<uMax<<" ";			
				ofs<<std::endl;	
			};

			//Convergence criteria
			if (stepInfo.Iteration == MaxIteration) {
				if (pManager->IsMaster()) {
					std::cout<<"Maximal number of iterations reached.\n";
				};
				break;
			};

			if (stepInfo.Time >= MaxTime) {
				if (pManager->IsMaster()) {
					std::cout<<"Maximal time reached.\n";
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
			std::cout<<"rank = "<<pManager->getRank()<<", local idle time = " << idleTime / 1e3 << " seconds.\n";
			std::cout.flush();
		};
		
		double maxIdleTime = pManager->Max(idleTime);

		if (pManager->IsMaster()) {
			//Close history file
			ofs.flush();
			ofs.close();

			std::cout<<"Calculation finished!\n";			
			//Stop timer
			workTimer.Stop();
			std::cout<<"Total work time = " << workTimer.ElapsedTimeMilliseconds() / 1e3 << " seconds.\n";			
			std::cout<<"Maximum idle time = " << maxIdleTime / 1e3 << " seconds.\n";			
		};
		
	};

	//Finalize kernel
	void Finilaze() {
		//Finalize MPI
		pManager->Finalize();
	};

};

#endif