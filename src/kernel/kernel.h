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
#include "ParallelManager.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "RiemannSolvers\RoeSolverPerfectGasEOS.h"
#include "BoundaryConditions\BCGeneral.h"
#include "Sensors\Sensors.h"

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
	//Interface part  TO DO REMOVE GetDummyCellLayerSize
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
	std::valarray<double> hx, hy, hz;	//cell sizes
	bool IsPeriodicX; //X periodicity
	bool IsPeriodicY; //Y periodicity
	bool IsPeriodicZ; //Z periodicity
	std::valarray<double> CoordinateX; //Cell center coordinates
	std::valarray<double> CoordinateY; //Cell center coordinates
	std::valarray<double> CoordinateZ; //Cell center coordinates

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
	bool isExternalAccelaration;
	bool isExternalForce;

	//sensors operating
	bool isSensorEnable;

	//External forces
	Vector Sigma; //Potential force like presure gradient
	Vector UniformAcceleration;	//external uniform acceleration

	//Boundary conditions
	std::unique_ptr<BoundaryConditions::BCGeneral> xLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> xRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> yRightBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zLeftBC;
	std::unique_ptr<BoundaryConditions::BCGeneral> zRightBC;

	//Solution data
	std::valarray<double> values;	
	std::valarray<double> residual;

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
	bool ContinueComputation;

	bool DebugOutputEnabled; //

	std::vector<std::shared_ptr<Sensor>> Sensors;

	//Constructor
	Kernel(int* argc, char **argv[]) : pManager(new ParallelManager(argc, argv)) {
			nVariables = 5;	//default value
	};

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
		//Output settings
		DebugOutputEnabled = config.DebugOutputEnabled;

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

		//nullify compression if we have uniform grid
		if(config.isUniformAlongX == false) {
			gridInfo.qx = config.qx;
			assert(gridInfo.IsPeriodicX == false);
		} else gridInfo.qx = 1;
		if(config.isUniformAlongY == false) {
			gridInfo.qy = config.qy;
			assert(gridInfo.IsPeriodicY == false);
		} else gridInfo.qy = 1;		
		if(config.isUniformAlongZ == false) {
			gridInfo.qz = config.qz;
			assert(gridInfo.IsPeriodicZ == false);
		} else gridInfo.qz = 1;

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
		hx.resize(nXAll);
		hy.resize(nYAll);
		hz.resize(nZAll);

		//fill cell centers positions and edges sizes
		double h_x = Lx / nX;			//uniform grid case
		if(gridInfo.qx != 1) h_x = 0.5 * Lx * (1.0 - gridInfo.qx) / (1.0 - pow(gridInfo.qx, 0.5 * nX));	//X step around the border
		double xl = - (dummyCellLayersX * h_x) + 0.5 * h_x;				//left cell (global) position
		double xr = Lx + dummyCellLayersX * h_x - 0.5 * h_x;			//right cell (global) position
		for(int i = 0; i < dummyCellLayersX; i++) {
			CoordinateX[i] = xl;
			CoordinateX[nXAll - 1 - i] = xr;
			hx[i] = h_x;
			hx[nXAll - 1 - i] = h_x;
			xl += h_x;
			xr -= h_x;
		};

		for (int i = iMin; i < 0.5 * (iMax + iMin + 1); i++) {
			CoordinateX[i] = xl;
			CoordinateX[nXAll - 1 - i] = xr;
			hx[i] = h_x;
			hx[nXAll - 1 - i] = h_x;
			xl += 0.5 * h_x * (1.0 + gridInfo.qx);
			xr -= 0.5 * h_x * (1.0 + gridInfo.qx);
			h_x *= gridInfo.qx;
		};

		double h_y = Ly / nY;			//uniform grid case
		if(gridInfo.qy != 1) h_y = 0.5 * Ly * (1.0 - gridInfo.qy) / (1.0 - pow(gridInfo.qy, 0.5 * nY));	//Y step around the border
		double yl = - (dummyCellLayersY * h_y) + 0.5 * h_y;			//	left cell (global) position
		double yr = Ly + dummyCellLayersY * h_y - 0.5 * h_y;		//	right cell (global) position
		for(int i = 0; i < dummyCellLayersY; i++) {
			CoordinateY[i] = yl;
			CoordinateY[nYAll - 1 - i] = yr;
			hy[i] = h_y;
			hy[nYAll - 1 - i] = h_y;
			yl += h_y;
			yr -= h_y;
		};

		for (int j = jMin; j < 0.5 * (jMax + jMin + 1); j++) {
			CoordinateY[j] = yl;
			CoordinateY[nYAll - 1 - j] = yr;
			hy[j] = h_y;
			hy[nYAll - 1 - j] = h_y;
			yl += 0.5 * h_y * (1.0 + gridInfo.qy);
			yr -= 0.5 * h_y * (1.0 + gridInfo.qy);
			h_y *= gridInfo.qy;
		};

		double h_z = Lz / nZ;			//uniform grid case
		if(gridInfo.qz != 1) h_z = 0.5 * Lz * (1.0 - gridInfo.qz) / (1.0 - pow(gridInfo.qz, 0.5 * nZ));	//Y step around the border
		double zl = - (dummyCellLayersZ * h_z) + 0.5 * h_z;				//left cell (global) position
		double zr = Lz + dummyCellLayersZ * h_z - 0.5 * h_z;			//right cell (global) position
		for(int i = 0; i < dummyCellLayersZ; i++) {
			CoordinateZ[i] = zl;
			CoordinateZ[nZAll - 1 - i] = zr;
			hz[i] = h_z;
			hz[nZAll - 1 - i] = h_z;
			zl += h_z;
			zr -= h_z;
		};

		for (int k = kMin; k < 0.5 * (kMax + kMin + 1); k++) {
			CoordinateZ[k] = zl;
			CoordinateZ[nZAll - 1 - k] = zr;
			hz[k] = h_z;
			hz[nZAll - 1 - k] = h_z;
			zl += 0.5 * h_z * (1.0 + gridInfo.qz);
			zr -= 0.5 * h_z * (1.0 + gridInfo.qz);
			h_z *= gridInfo.qz;
		};

		if (nDims < 2) hy[0] = 1.0;
		if (nDims < 3) hz[0] = 1.0;

		//Initialize gas model parameters and Riemann solver
		gamma = config.Gamma;
		viscosity = config.Viscosity;
		thermalConductivity = config.ThermalConductivity;
		if(config.IsViscousFlow == true) {
			isViscousFlow = true;
			isGradientRequired = true;
		} else {
			isViscousFlow = false;
			isGradientRequired = false;
		};
		if(config.IsExternalForceRequared == true) {
			isExternalForce = true;
		} else {
			isExternalForce = false;
		};
		if(config.IsUnifromAccelerationRequared == true) {
			isExternalAccelaration = true;
		} else {
			isExternalAccelaration = false;
		};
		ContinueComputation = config.ContinueComputation;

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
		InitBoundaryConditions(config);

		//External forces
		Sigma = config.Sigma;
		UniformAcceleration = config.UniformAcceleration;

		isSensorEnable = false;	//by default
		
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Initialize boundary conditions
	virtual void InitBoundaryConditions(KernelConfiguration& config) {
		if (!IsPeriodicX) {
			xLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yRightBC = std::unique_ptr<BoundaryConditions::BCGeneral>(new BoundaryConditions::BCGeneral());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};
		if ((!IsPeriodicZ) && (nDims > 2)) {
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
		
		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = hx[i] * hy[j] * hz[k];
					//Update cell values
					for(int nv = 0; nv < nVariables; nv++) {
						values[idx * nVariables + nv] += residual[idx * nVariables + nv] * dt / volume;
						//Compute total residual
						stepInfo.Residual[nv] += abs(residual[idx * nVariables + nv]);
					};
				};
			};
		};

		//Aggregate
		for(int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = pManager->Sum(stepInfo.Residual[nv]);
	};

	//Update Residual
	void UpdateResiduals() {
		stepInfo.Residual.resize(nVariables);
		for (int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = 0;
		
		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = hx[i - iMin + GetDummyCellLayerSize()] * hy[j - jMin + GetDummyCellLayerSize()] * hz[k - kMin + GetDummyCellLayerSize()];
					//Update cell values
					for(int nv = 0; nv < nVariables; nv++) {
						//Compute total residual
						stepInfo.Residual[nv] += abs(residual[idx * nVariables + nv]);
					};
				};
			};
		};

		//Aggregate
		for(int nv = 0; nv < nVariables; nv++) stepInfo.Residual[nv] = pManager->Sum(stepInfo.Residual[nv]);
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
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<rank<<
				"; rankL = "<<rankL<<
				"; rankR = "<<rankR<<
				std::endl<<std::flush;
		};

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
						double *U = getCellValues(iSend, j, k);
						for (int nv = 0; nv < nVariables; nv++) bufferToSend[idxBuffer * nVariables + nv] = U[nv];
					};
				};
			};

			if (DebugOutputEnabled) {
				std::cout<<"rank = "<<rank<<
				"; iSend = "<<iSend<<
				"; iRecv = "<<iRecv<<
				std::endl<<std::flush;
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
					double *U = getCellValues(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) U[nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			iSend = iMax - layer + 1; // layer index to send
			iRecv = iMin - layer; // layer index to recv

			if (DebugOutputEnabled) {
				std::cout<<"rank = "<<rank<<
				"; iSend = "<<iSend<<
				"; iRecv = "<<iRecv<<
				std::endl<<std::flush;
			};

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (j = jMin; j <= jMax; j++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
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
			for (j = jMin; j <= jMax; j++) {
				for (k = kMin; k <= kMax; k++) {
					int idxBuffer = (j-jMin) + (k-kMin)* nY; //Exclude x index
					double *U = getCellValues(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) U[nv] = bufferToRecv[idxBuffer * nVariables + nv];
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
		/*std::cout<<"rank = "<<rank<<
			"; rankL = "<<rankL<<
			"; rankR = "<<rankR<<
			std::endl<<std::flush;*/

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
		//rankL = pManager->GetRankByCartesianIndex(0, 0, -1);
		//rankR = pManager->GetRankByCartesianIndex(0, 0, 1);

		if (DebugOutputEnabled) {
			std::cout<<"Exchange Z-direction started"<<std::endl<<std::flush;
			std::cout<<"rank = "<<rank<<
				"; rankL = "<<rankL<<
				"; rankR = "<<rankR<<
				std::endl<<std::flush;
		};
		

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

			if (DebugOutputEnabled) {
				std::cout<<"rank = "<<rank<<
				"; kSend = "<<kSend<<
				"; kRecv = "<<kRecv<<
				std::endl<<std::flush;
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

			if (DebugOutputEnabled) {
				std::cout<<"rank = "<<rank<<
				"; kSend = "<<kSend<<
				"; kRecv = "<<kRecv<<
				std::endl<<std::flush;
			};

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

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<rank<<"Exchange values finished."<<			
			std::endl<<std::flush;
		};

	}; // function

	//Compute dummy cell values as result of boundary conditions and interprocessor exchange communication
	void ComputeDummyCellValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//Interprocessor exchange
		ExchangeValues();

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<rank<<"Dummy cell processing started."<<			
			std::endl<<std::flush;
		};

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

					if (DebugOutputEnabled) {
						std::cout<<"rank = "<<rank<<
						"; f.y = "<<faceCenter.y<<
						"; f.z = "<<faceCenter.z<<
						std::endl<<std::flush;
					};

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
		if (!IsPeriodicY) {
			for (i = iMin; i <= iMax; i++) {
				for (k = kMin; k <= kMax; k++) {				
					//Inner cell
					cellCenter.x = CoordinateX[i];
					cellCenter.z = CoordinateZ[k];

					//And face
					faceCenter.x = CoordinateX[i];
					faceCenter.z = CoordinateZ[k];

					for (int layer = 1; layer <= dummyCellLayersY; layer++) {
						if (pManager->rankCart[1] == 0) {
							//Left border
							j = jMin - layer; // layer index
							int jIn = jMin + layer - 1; // opposite index
							cellCenter.y = CoordinateY[jIn];
							faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;	//for uniform dummy cells		TO DO MODIFY
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

		if (nDims < 3) return;
		//Z direction		
		faceNormalL = Vector(0.0, 0.0, -1.0);
		faceNormalR = Vector(0.0, 0.0, 1.0);
		if (!IsPeriodicZ) {
			for (i = iMin; i <= iMax; i++) {
				for (j = jMin; j <= jMax; j++) {				
					//Inner cell
					cellCenter.x = CoordinateX[i];
					cellCenter.y = CoordinateY[j];

					//And face
					faceCenter.x = CoordinateX[i];
					faceCenter.y = CoordinateY[j];

					for (int layer = 1; layer <= dummyCellLayersY; layer++) {

						if (pManager->rankCart[2] == 0) {
							//Left border
							k = kMin - layer; // layer index
							int kIn = kMin + layer - 1; // opposite index
							cellCenter.z = CoordinateZ[kIn];
							faceCenter.z = (CoordinateZ[kIn] + CoordinateZ[k]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, j, kIn);
					
							//Apply left boundary conditions						
							std::vector<double> dValues = zLeftBC->getDummyValues(&values[idxIn * nVariables], faceNormalL, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};

						if (pManager->rankCart[2] == pManager->dimsCart[2]-1) {
							//Right border
							k = kMax + layer; // layer index
							int kIn = kMax - layer + 1; // opposite index
							cellCenter.z = CoordinateZ[kIn];
							faceCenter.z = (CoordinateZ[kIn] + CoordinateZ[k]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, j, kIn);
					
							//Apply right boundary conditions						
							std::vector<double> dValues = zRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
							for (int nv = 0; nv < nVariables; nv++) {
								values[idx * nVariables + nv] = dValues[nv];
							};
						};
					};
				};
			};
		}; //Z direction

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Dummy values Z-direction computed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//compute source terms
	virtual void ProcessExternalForces() {
		//right handside part - mass foces
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = hx[i] * hy[j] * hz[k];

					double *V = getCellValues(i, j, k);
					double u = V[1]/V[0];
					double v = V[2]/V[0];
					double w = V[3]/V[0];
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
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					int idx = getSerialIndexLocal(i, j, k);
					//Compute cell volume
					double volume = hx[i] * hy[j] * hz[k];
					double *V = getCellValues(i, j, k);
					double u = V[1]/V[0];
					double v = V[2]/V[0];
					double w = V[3]/V[0];
					Vector velocity = Vector(u, v, w); //Cell velocity

					//Mass forces represents body accelerations (per unit mass)
					double ro = V[0];
					residual[idx * nVariables + 1] += ro * volume * UniformAcceleration.x;				//rou
					residual[idx * nVariables + 2] += ro * volume * UniformAcceleration.y;				//rov
					residual[idx * nVariables + 3] += ro * volume * UniformAcceleration.z;				//row
					residual[idx * nVariables + 4] += ro * UniformAcceleration * velocity * volume ;	//roE
				};
			};
		};
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

	//Save solution TecPlot Format 
	void SaveSliceToTecplot(std::string fname, int I, int J, int K) {
		
		std::ofstream ofs;
		ofs<<std::scientific;
		int rank = pManager->getRank();

		//Open file
		if (pManager->IsMaster()) {
			//std::cout<<"File created"<<std::endl<<std::flush;
			ofs.open(fname, std::ios_base::out);
			ofs<<"VARIABLES = ";
			int dims = 0;
			if (I == -1) {
				ofs<<"\""<<"X"<<"\" ";
				dims++;		
			};
			if (J == -1) {
				ofs<<"\""<<"Y"<<"\" ";
				dims++;
			};
			if (K == -1) {
				ofs<<"\""<<"Z"<<"\" ";
				dims++;
			};
							
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			ofs<<"\""<<"v"<<"\" ";
			ofs<<"\""<<"w"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"e"<<"\" ";
			ofs<<std::endl;

			ofs << "ZONE T=\"1\"";	//zone name
			ofs<<std::endl;

			if (I == -1) ofs << "I=" << nX;
			else ofs << "I=" << 1;
			if (J == -1) ofs << "J=" << nY;
			else ofs << "J=" << 1;
			if (K == -1) ofs << "K=" << nZ;
			else ofs << "K=" << 1;
			ofs << "F=POINT\n";
			ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
		} else {
			//Wait for previous process to finish writing
			pManager->Wait(rank - 1);

			//Open file for modification
			//std::cout<<"File opened"<<std::endl<<std::flush;
			ofs.open(fname, std::ios_base::app);
		};

		//Solution output
		for (int i = iMin; i <= iMax; i++) {
			if ((I != -1) && (i != I)) continue;
			for (int j = jMin; j <= jMax; j++) {
				if ((J != -1) && (j != J)) continue;
				for (int k = kMin; k <= kMax; k++) {
					if ((K != -1) && (k != K)) continue;
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];					
					double* U = getCellValues(i,j,k);
					double ro = U[0];
					double u = U[1] / ro;
					double v = U[2] / ro;
					double w = U[3] / ro;
					double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
					double P = (gamma - 1.0) * ro * e;

					//Write to file
					if (I == -1) ofs<<x<<" ";
					if (J == -1) ofs<<y<<" ";
					if (K == -1) ofs<<z<<" ";					
					ofs<<ro<<" ";
					ofs<<u<<" ";
					ofs<<v<<" ";
					ofs<<w<<" ";
					ofs<<P<<" ";
					ofs<<e<<" ";
					ofs<<std::endl;
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

	void SaveSolution(std::string fname) {
		SaveSliceToTecplot(fname, -1, -1, 0);

		return;

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

			if (nDims > 1) ofs<<"J = "<<nY<<", ";
			ofs<<"I = "<<nX<<", ";
			
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
			MPI_Send(&values[0], values.size(), MPI_LONG_DOUBLE, masterRank, 0, pManager->_commCart);
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
							recvBuff = std::vector<double>(std::begin(values), std::end(values));								
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

	virtual void SaveSolutionSega(std::string fname) {
		//Tecplot version
		std::ofstream ofs(fname);

		//1D tecplot style
		if(nDims == 1)
		{
			//Header				
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"X"<<"\" ";
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			ofs<<"\""<<"v"<<"\" ";
			ofs<<"\""<<"w"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"e"<<"\" ";
			ofs<<std::endl;
			
			//Solution
			for (int i = iMin; i <= iMax; i++) {
				//Obtain cell data
				double x = CoordinateX[i];
				double* U = getCellValues(i, jMin, kMin);
				double ro = U[0];
				double u = U[1] / ro;
				double v = U[2] / ro;
				double w = U[3] / ro;
				double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
				double P = (gamma - 1.0) * ro * e;

				//Write to file
				ofs<<x<<" ";
				ofs<<ro<<" ";
				ofs<<u<<" ";
				ofs<<v<<" ";
				ofs<<w<<" ";
				ofs<<P<<" ";
				ofs<<e<<" ";
				ofs<<std::endl;
			};	//end cycle

			ofs.close();
			return;
		};	//end if

		//2D/3D tecplot style
		if(nDims > 1)
		{
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
		
			//Solution
			for (int k = kMin; k <= kMax; k++) {
				for (int j = jMin; j <= jMax; j++) {
					for (int i = iMin; i <= iMax; i++) {
						//Obtain cell data
						double x = CoordinateX[i];
						double y = CoordinateY[j];
						double z = CoordinateZ[k];
						double* U = getCellValues(i,j,k);
						double ro = U[0];
						double u = U[1] / ro;
						double v = U[2] / ro;
						double w = U[3] / ro;
						double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
						double P = (gamma - 1.0) * ro * e;

						//Write to file
						ofs<<x<<" ";
						ofs<<y<<" ";
						ofs<<z<<" ";
						ofs<<ro<<" ";
						ofs<<u<<" ";
						ofs<<v<<" ";
						ofs<<w<<" ";
						ofs<<P<<" ";
						ofs<<e<<" ";
						ofs<<std::endl;
					};
				};
			};	//end for ijk
		};	//end if

		ofs.close();
		return;
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

		//Start timers
		Timer workTimer;
		workTimer.Start();
		workTimer.Resume();

		Timer iterTimer;
		iterTimer.Start();

		//Calc loop		
		if (pManager->IsMaster()) {
			std::cout<<"Calculation started!"<<std::endl<<std::flush;
		};
		pManager->Barrier();

		//Processor info output
		rank = pManager->getRank();
		if (!pManager->IsFirstNode()) pManager->Wait(rank - 1);
		std::cout<<"Info for node rank : "<<rank<<""<<std::endl<<
			"rankX = "<<pManager->rankCart[0]<<
			", iMin = "<<iMin<<
			", iMax = "<<iMax
			<<std::endl;
			if (nDims > 1) {
				std::cout<<
				"rankY = "<<pManager->rankCart[1]<<
				", jMin = "<<jMin<<
				", jMax = "<<jMax
				<<std::endl;
			};
			if (nDims > 2) {
				std::cout<<
				"rankZ = "<<pManager->rankCart[2]<<
				", jMin = "<<jMin<<
				", jMax = "<<jMax
				<<std::endl;
			};
		if (!pManager->IsLastNode()) pManager->Signal(rank + 1);
		pManager->Barrier();

		for (int iteration = 0; iteration <= MaxIteration; iteration++) {
			//Exchange and boundary conditions
			ComputeDummyCellValues();

			//Calculate one time step
			IterationStep();

			//Output step information						
			if (pManager->IsMaster() && (ResidualOutputIterations != 0) && (stepInfo.Iteration % ResidualOutputIterations) == 0) {
				iterTimer.Stop();
				double iterationsPackTimeAverage = iterTimer.ElapsedTimeMilliseconds();
				iterationsPackTimeAverage /= ResidualOutputIterations;
				iterTimer.Start();
				std::cout<<"Iteration = "<<stepInfo.Iteration<<
					"; Total time = "<< stepInfo.Time <<
					"; Time step = " <<stepInfo.TimeStep <<
					"; RMSrou = "<<stepInfo.Residual[1]<<
					"; avgIterTime = "<<iterationsPackTimeAverage<<" ms"<<
					std::endl;
			};		

			//Sensors recording
			if(isSensorEnable == true)
			{
				for(auto &r : Sensors) {
					r->NewRecord();
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
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".dat";						
				SaveSolutionSega(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout<<"Solution has been written to file \""<<snapshotFileName.str()<<"\""<<std::endl;
				};
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (stepInfo.NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);								
				snapshotFileName<<"dataT"<<stepInfo.Time<<".dat";							
				SaveSolutionSega(snapshotFileName.str());

				if (pManager->IsMaster()) {
					std::cout<<"Solution has been written to file \""<<snapshotFileName.str()<<"\""<<std::endl;
				};

				//Adjust next snapshot time
				stepInfo.NextSnapshotTime += SaveSolutionSnapshotTime;
			};

			//Save convergence history		
			if (pManager->IsMaster() && (stepInfo.Iteration % ResidualOutputIterations)) {
				ofs<<stepInfo.Time<<" ";			
				ofs<<stepInfo.Residual[0]<<" ";							
				ofs<<stepInfo.Residual[1]<<" ";			
				if (nDims > 1) ofs<<stepInfo.Residual[2]<<" ";			
				if (nDims > 2) ofs<<stepInfo.Residual[3]<<" ";		
				if (nVariables > 4) ofs<<stepInfo.Residual[4]<<" ";				
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

			//Now master node outputs connectivity information
			if (pManager->IsMaster()) {
			};
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