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
#include "BoundaryConditions/BCGeneral.h"

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
	int nCellsLocal; //Total number of local cells including dummy
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

	//Boundary conditions
	std::unique_ptr<BCGeneral> xLeftBC;
	std::unique_ptr<BCGeneral> xRightBC;
	std::unique_ptr<BCGeneral> yLeftBC;
	std::unique_ptr<BCGeneral> yRightBC;
	std::unique_ptr<BCGeneral> zLeftBC;
	std::unique_ptr<BCGeneral> zRightBC;

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
		//std::cout<<rank<<" "<<iMin<<" "<<iMax<<"\n";
		//exit(0);
		pManager->Barrier();

		double Lx = config.LX;
		double Ly = config.LY;
		double Lz = config.LZ;		

		//Generate cells (uniform grid)
		nXAll = nX + 2 * dummyCellLayersX;
		nYAll = nY + 2 * dummyCellLayersY;
		nZAll = nZ + 2 * dummyCellLayersZ;
		nlocalXAll = nlocalX + 2 * dummyCellLayersX;
		nlocalYAll = nlocalY + 2 * dummyCellLayersY;
		nlocalZAll = nlocalZ + 2 * dummyCellLayersZ;
		nCellsLocal = nlocalXAll * nlocalYAll * nlocalZAll;
		CoordinateX.resize(nXAll);
		CoordinateY.resize(nYAll);
		CoordinateZ.resize(nZAll);
		double dx = Lx / nX;
		hx = dx;
		double xMin = 0;
		xMin = xMin - (dummyCellLayersX * dx) + 0.5 * dx;
		double xMax = Lx;
		xMax = xMax + (dummyCellLayersX * dx) - 0.5 * dx;
		for (int i = iMin - dummyCellLayersX; i < iMax + dummyCellLayersX; i++) {
			double x = xMin + (xMax - xMin) * 1.0 * i / (nlocalXAll - 1);
			CoordinateX[i] = x;
		};
		CoordinateX[iMax + dummyCellLayersX] = xMax;

		double dy = Ly / nY;
		hy = dy;
		double yMin = 0;
		yMin = yMin - (dummyCellLayersY * dy) + 0.5 * dy;
		double yMax = Ly;
		yMax = yMax + (dummyCellLayersY * dy) - 0.5 * dy;
		for (int j = jMin - dummyCellLayersY; j < jMax + dummyCellLayersY; j++) {
			double y = yMin + (yMax - yMin) * 1.0 * j / (nlocalYAll - 1);					
			CoordinateY[j] = y;
		};
		CoordinateY[jMax + dummyCellLayersY] = yMax;

		double dz = Lz / nZ;
		double zMin = 0;
		zMin = zMin - (dummyCellLayersZ * dz) + 0.5 * dz;
		double zMax = Lz;
		zMax = zMax + (dummyCellLayersZ * dz) - 0.5 * dz;
		for (int k = kMin - dummyCellLayersZ; k < kMax + dummyCellLayersZ; k++) {
			double z = zMin + (zMax - zMin) * 1.0 * k / (nlocalZAll - 1);
			CoordinateZ[k] = z;
		};
		CoordinateZ[kMax + dummyCellLayersZ] = zMax;

		//Cell linear sizes
		hx = hy = hz = 1.0;
		hx = dx;
		if (nDims > 1) hy = dy;
		if (nDims > 2) hz = dz;

		//Initialize gas model parameters and riemann solver
		gamma = config.gamma;
		nVariables = config.nVariables;	

		//Allocate data structures
		values.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);	
		residual.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);	

		//Initialize calculation parameters
		MaxTime = config.MaxTime;
		MaxIteration = config.MaxIteration;
		SaveSolutionSnapshotTime = config.SaveSolutionSnapshotTime;	
		SaveSolutionSnapshotIterations = config.SaveSolutionSnapshotIterations;
		stepInfo.Time = 0;
		stepInfo.Iteration = 0;
		stepInfo.NextSnapshotTime = stepInfo.Time;

		//Initialize boundary conditions
		if (!gridInfo.IsPeriodicX) {
			xLeftBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			xRightBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!gridInfo.IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			yRightBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};

		if ((!gridInfo.IsPeriodicZ) && (nDims > 2)) {
			zLeftBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			zRightBC = std::unique_ptr<BCGeneral>(new BCGeneral());
			zLeftBC->loadConfiguration(config.zLeftBoundary);
			zRightBC->loadConfiguration(config.zRightBoundary);
		};
	};

	//Update solution
	void UpdateSolution(double dt) {
		//Compute cell volume
		double volume = hx * hy * hz;
		stepInfo.Residual.resize(5, 0);

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
						int idxBuffer = (j + k * nY); //Exclude x index
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
					int idxBuffer = (j + k * nY); //Exclude x index
					int idxValues = getSerialIndexLocal(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			iSend = iMax - layer; // layer index to send
			iRecv = iMin - layer; // layer index to recv

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (j = jMin; j <= jMax; j++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (j + k * nY); //Exclude x index
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
					int idxBuffer = (j + k * nY); //Exclude x index
					int idxValues = getSerialIndexLocal(iRecv, j, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 

		//
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
						int idxBuffer = (i + k * nX); //Exclude x index
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
					int idxBuffer = (i + k * nX); //Exclude x index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

			// Plus direction exchange
			nSend = 0; //
			nRecv = 0; //
			jSend = jMax - layer; // layer index to send
			jRecv = jMin - layer; // layer index to recv

			// Prepare values to send
			if (rankR != -1) {
				nSend = layerSize * nVariables;
				for (i = iMin; i <= iMax; i++) {
					for (k = kMin; k <= kMax; k++) {
						int idxBuffer = (i + k * nX); //Exclude x index
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
					int idxBuffer = (i + k * nX); //Exclude x index
					int idxValues = getSerialIndexLocal(i, jRecv, k);
					for (int nv = 0; nv < nVariables; nv++) values[idxValues * nVariables + nv] = bufferToRecv[idxBuffer * nVariables + nv];
				};
			};

		}; // 


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
		for (j = jMin; j <= jMax; j++) {
			for (k = kMin; k <= kMax; k++) {				
				//Inner cell
				cellCenter.y = CoordinateY[j];
				cellCenter.z = CoordinateZ[k];

				//And face
				faceCenter.y = CoordinateY[j];
				faceCenter.z = CoordinateZ[k];

				for (int layer = 1; layer <= dummyCellLayersX; layer++) {
					if (!IsPeriodicX) {
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

						//Right border
						i = iMax + layer; // layer index
						iIn = iMax - layer + 1; // opposite index
						cellCenter.x = CoordinateX[iIn];
						faceCenter.x = (CoordinateX[iIn] + CoordinateX[i]) / 2.0;
						idx = getSerialIndexLocal(i, j, k);
						idxIn = getSerialIndexLocal(iIn, j, k);
					
						//Apply right boundary conditions						
						dValues = xRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
						for (int nv = 0; nv < nVariables; nv++) {
							values[idx * nVariables + nv] = dValues[nv];
						};
					};
				};
			};
		}; //X direction

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
				faceCenter.x = CoordinateY[i];
				faceCenter.z = CoordinateZ[k];

				for (int layer = 1; layer <= dummyCellLayersY; layer++) {
					if (!IsPeriodicY) {
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

						//Right border
						j = jMax + layer; // layer index
						jIn = jMax - layer + 1; // opposite index
						cellCenter.y = CoordinateY[jIn];
						faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;
						idx = getSerialIndexLocal(i, j, k);
						idxIn = getSerialIndexLocal(i, jIn, k);
					
						//Apply right boundary conditions						
						dValues = yRightBC->getDummyValues(&values[idxIn * nVariables], faceNormalR, faceCenter, cellCenter);
						for (int nv = 0; nv < nVariables; nv++) {
							values[idx * nVariables + nv] = dValues[nv];
						};
					};
				};
			};
		}; //Y direction

	};

	//Gather 
	//void GatherField(std::vector<double> result, std::vector<double> local, std::function<double(double *)> func) {
	//	if (pManager->IsMaster()) {
	//	for (int i = iMin; i <= iMax; i++) {
	//		for (int j = jMin; j <= jMax; j++) {
	//			for (int k = kMin; k <= kMax; k++) {
	//				//Obtain cell data
	//				int indexGlobal = getSerialIndexGlobal(i, j, k);
	//				int indexLocal = getSerialIndexLocal(i, j, k);
	//				double* U = getCellValues(i,j,k);
	//				double value = func(U);					
	//			};
	//		};
	//	};
	//};

	//Save solution
	void SaveSolution(std::string fname) {
		std::ofstream ofs;

		{
			//Header
		if (pManager->IsMaster()) {
			ofs.open(fname, std::ios_base::out);
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"X"<<"\" ";
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"e"<<"\" ";
			ofs<<std::endl;
		};

		//Send solution to master
		if (!pManager->IsMaster()) {			
			MPI_Send(&values.front(), values.size(), MPI_LONG_DOUBLE, 0, 0, pManager->_commCart);
		};

		//On master output each part
		if (pManager->IsMaster()) {
			std::vector<double> recvBuff;			
			for (int rI = 0; rI < pManager->dimsCart[0]; rI++) {
				for (int rJ = 0; rJ < pManager->dimsCart[1]; rJ++) {
					for (int rK = 0; rK < pManager->dimsCart[2]; rK++) {
						if ((rI == pManager->rankCart[0]) && (rJ == pManager->rankCart[1]) && (rK == pManager->rankCart[2])) {
							//Do not recieve
						} else {
							//Recieve
							recvBuff.resize(values.size());
							int rank = pManager->GetRankByCartesianIndex(rI, rJ, rK);
							MPI_Status status;
							MPI_Recv(&recvBuff.front(), values.size(), MPI_LONG_DOUBLE, rank, 0, pManager->_commCart, &status);
						};
						
						//Solution output
						for (int i = iMin; i <= iMax; i++) {
							for (int j = jMin; j <= jMax; j++) {
								for (int k = kMin; k <= kMax; k++) {
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
									ofs<<ro<<" ";
									ofs<<u<<" ";
									ofs<<P<<" ";
									ofs<<e<<" ";
									ofs<<std::endl;
								};
							};
						}; // Solution output
					}; // for rK
				}; // for rJ
			}; // for rI

		};// if	

		//Gather solution
		if (pManager->IsMaster()) ofs.close();

		//Sync
		pManager->Barrier();
		};

		//2D tecplot style
		if(nDims>1)
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
			};
		};
		std::cout << "Solution is written\n";	
	};

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = 1 - std::floor(std::log10(SaveSolutionSnapshotTime));
		};

		//Start timer
		Timer workTimer;
		workTimer.Start();

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
			if (pManager->IsMaster()) {
				std::cout<<"Iteration = "<<stepInfo.Iteration<<"; Total time = "<< stepInfo.Time << "; Time step = " <<stepInfo.TimeStep << "; RMSrou = "<<stepInfo.Residual[1]<<"\n";			
			};			

			//Solution snapshots
			//Every few iterations
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".cgns";						
				SaveSolution(snapshotFileName.str());
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (stepInfo.NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);								
				snapshotFileName<<"dataT"<<stepInfo.Time<<".cgns";							
				SaveSolution(snapshotFileName.str());

				//Adjust next snapshot time
				stepInfo.NextSnapshotTime += SaveSolutionSnapshotTime;
			};

			//Save history			

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
		if (pManager->IsMaster()) {
			std::cout<<"Calculation finished!\n";
		};

		//Stop timer
		workTimer.Stop();
		std::cout<<"Total work time = " << workTimer.ElapsedTimeMilliseconds() / 1e3 << " seconds.\n";			
	};

	//Finalize kernel
	void Finilaze() {
		//Finalize MPI
	};

};

#endif