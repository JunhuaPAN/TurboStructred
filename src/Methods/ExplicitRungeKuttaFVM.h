#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "Methods\Method.h"
#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "RiemannSolvers\RoeSolverPerfectGasEOS.h"
#include "kernel.h"

class Kernel;

//Base class for all solution methods that desribe iterations process in detail
class ExplicitRungeKuttaFVM : public Method {
	//Current riemann problem solver
	std::unique_ptr<RiemannSolver> _riemannSolver;
public:
	//Method parameters
	int RungeKuttaOrder;
	int nVariables;
	double CFL;

	//Parallel run information
	int rank; //Rank of process
	int rankCart[3]; //Position of process in cartesian mpi grid
	int dimsCart[3]; //Dimensions of cartesian mpi grid


	//Additional data
	std::vector<double> spectralRadius;		//array for storing spectral radiuses

	//Required parameters from kernel for faster work of inline functions
	int iMin; 
	int iMax;
	int jMin;
	int jMax;
	int kMin;
	int kMax;
	int nlocalXAll;
	int nlocalYAll;
	int nlocalZAll;
	int dummyCellLayersX;
	int dummyCellLayersY;
	int dummyCellLayersZ;

	//Grid information
	int nDims;
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
	double hx, hy, hz;	//cell sizes
	bool IsPeriodicX; //X periodicity
	bool IsPeriodicY; //Y periodicity
	bool IsPeriodicZ; //Z periodicity

	//Number of required dummy layers
	virtual int GetDummyCellLayerSize() override {
		return 1;
	};

	//Get serial index
	inline int getSerialIndexLocal(int i, int j, int k) {
		int sI =  (k - kMin + dummyCellLayersZ) * nlocalXAll * nlocalYAll + (j - jMin + dummyCellLayersY) * nlocalXAll + (i - iMin + dummyCellLayersX);
		return sI;
	};

	//Get cell values
	inline double* getCellValues(int i, int j, int k) {
		int sBegin = getSerialIndexLocal(i, j, k) * nVariables;		
		return &values[sBegin];
	};

	//Initizalization
	void Init(KernelConfiguration& kernelConfig) {
		MethodConfiguration config = kernelConfig.methodConfiguration;
		_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.gamma, 0.05, 0.0));
		CFL = config.CFL;
		nVariables = kernelConfig.nVariables;

		//Initialize MPI
		rank = 0;
		rankCart[0] = 0;
		rankCart[1] = 0;
		rankCart[2] = 0;
		dimsCart[0] = 1;
		dimsCart[1] = 1;
		dimsCart[2] = 1;

		//Init params for fast memory access
		nDims = kernelConfig.nDims;

		//Initialize local grid
		int dummyCellLayers = GetDummyCellLayerSize(); //number of dummy cell layers				
		nX = kernelConfig.nX;
		IsPeriodicX = kernelConfig.isPeriodicX;
		dummyCellLayersX = dummyCellLayers;
		nY = 1;
		IsPeriodicY = true;
		dummyCellLayersY = 0;
		nZ = 1;
		IsPeriodicZ = true;
		dummyCellLayersZ = 0;
		if (nDims > 1) {
			nY = kernelConfig.nY;
			IsPeriodicY = kernelConfig.isPeriodicY;
			dummyCellLayersY = dummyCellLayers;
		};
		if (nDims > 2) {
			nZ = kernelConfig.nZ;
			IsPeriodicZ = kernelConfig.isPeriodicZ;
			dummyCellLayersZ = dummyCellLayers;
		};
		
		nlocalX = nX / dimsCart[0];
		nlocalY = nY / dimsCart[1];
		nlocalZ = nZ / dimsCart[2];	
		iMin = rankCart[0] * nlocalX + dummyCellLayersX;
		iMax = (rankCart[0]+1) * nlocalX + dummyCellLayersX - 1;
		jMin = rankCart[1] * nlocalY + dummyCellLayersY;
		jMax = (rankCart[1]+1) * nlocalY + dummyCellLayersY - 1;
		kMin = rankCart[2] * nlocalZ + dummyCellLayersZ;
		kMax = (rankCart[2]+1) * nlocalZ + dummyCellLayersZ - 1;
		double Lx = kernelConfig.LX;
		double Ly = kernelConfig.LY;
		double Lz = kernelConfig.LZ;		

		//Generate cells (uniform grid)
		nlocalXAll = nlocalX + 2 * dummyCellLayersX;
		nlocalYAll = nlocalY + 2 * dummyCellLayersY;
		nlocalZAll = nlocalZ + 2 * dummyCellLayersZ;
		nCellsLocal = nlocalXAll * nlocalYAll * nlocalZAll;
		double dx = Lx / nX;
		double dy = Ly / nY;
		double dz = Lz / nZ;

		//Cell linear sizes
		hx = hy = hz = 1.0;
		hx = dx;
		if (nDims > 1) hy = dy;
		if (nDims > 2) hz = dz;

		//Allocate memory for values and residual
		spectralRadius.resize(nCellsLocal);
		values.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);	
		residual.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);	
	};

	//Compute residual
	void ComputeResidual(const std::vector<double> values, std::vector<double>& residual, std::vector<double>& spectralRadius) {
		//  init spectral radius storage for each cell
		for (double& sr : spectralRadius) sr = 0; //Nullify

		// Stencil size and values references
		int stencilSizeX = 2;
		std::vector<double*> U(stencilSizeX);
		
		// Fluxes temporary storage
		std::vector<double> fl(5,0); //left flux -1/2
		std::vector<double> fr(5,0); //right flux +1/2

		// I step
		Vector fn = Vector(1.0, 0.0, 0.0); // x direction
		for (int k = kMin; k <= kMax; k++) {
			for (int j = jMin; j <= jMax; j++) {
				//Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = hy;
				if (nDims == 3) fS = hy * hz;

				//Initial load of values we assume symmetric stencil
				for (int s = 0; s < stencilSizeX - 1; s++) {
					U[s] = getCellValues(iMin - stencilSizeX / 2 + s, j, k);
				};

				for (int i = iMin; i <= iMax + 1; i++)
				{
					//Update stencil values
					U[stencilSizeX - 1] = getCellValues(i + stencilSizeX / 2 - 1, j, k);

					//Compute flux
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
					fr = result.Fluxes;
	
					//Update residuals
					if(i > iMin)
					{
						//Fluxes difference equals residual 
						int idx = getSerialIndexLocal(i - 1, j, k);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] = - (fr[nv] - fl[nv]) * fS;

						//Add up spectral radius estimate
						spectralRadius[idx] += fS * result.MaxEigenvalue;
					};
          
					//Shift stencil
					fl = fr;
					for (int s = 1; s < stencilSizeX; s++) U[s - 1] = U[s];
				 };
			};
		};
		

		// J step
		if (nDims < 2) return;

		// K step
		if (nDims < 3) return;
	};

	//Compute global time step
	double ComputeTimeStep(std::vector<double>& spectralRadius) {
		//Compute cell volume
		double volume = hx * hy * hz;
		double dt = std::numeric_limits<double>::max();
		for (int cellIndex = 0; cellIndex< nCellsLocal; cellIndex++)
		{		
			double sR = spectralRadius[cellIndex];
			double localdt = CFL * volume / sR; //Blazek f. 6.20

			//Find minimum
			if (dt > localdt) dt = localdt;
		}

		dt = std::min(stepInfo->NextSnapshotTime - stepInfo->Time, dt);

		return dt;
	};
	
	//Update solution
	void UpdateSolution(double dt) {
		//Compute cell volume
		double volume = hx * hy * hz;
		stepInfo->Residual.resize(5, 0);

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
					stepInfo->Residual[0] += abs(residual[idx * nVariables]);			//ro
					stepInfo->Residual[1] += abs(residual[idx * nVariables + 1]);		//rou
					stepInfo->Residual[2] += abs(residual[idx * nVariables + 2]);		//rov
					stepInfo->Residual[3] += abs(residual[idx * nVariables + 3]);		//row
					stepInfo->Residual[4] += abs(residual[idx * nVariables + 4]);		//roE
				};
			};
		};
	};

	//Explicit time step
	virtual void IterationStep(StepInfo& _stepInfo) override {	
		stepInfo = &_stepInfo;

		//Compute residual
		ComputeResidual(values, residual, spectralRadius);

		//Determine timestep
		stepInfo->TimeStep = ComputeTimeStep(spectralRadius);

		//Compute residuals and update solution
		UpdateSolution(stepInfo->TimeStep);
		
		// Update time
		stepInfo->Time += stepInfo->TimeStep;
		stepInfo->Iteration++;
	};	

};

#endif