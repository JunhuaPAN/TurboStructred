#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "Methods\Method.h"
#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "RiemannSolvers\RoeSolverPerfectGasEOS.h"
#include "kernel.h"

//Base class for all solution methods that desribe iterations process in detail
class ExplicitRungeKuttaFVM : public Kernel {
	//Current riemann problem solver
	std::unique_ptr<RiemannSolver> _riemannSolver;
public:
	//Method parameters
	int RungeKuttaOrder;
	double CFL;

	//Additional data
	std::vector<double> spectralRadius;		//array for storing spectral radiuses

	//Number of required dummy layers
	virtual int GetDummyCellLayerSize() override {
		return 1;
	};

	//Initizalization
	void Init(KernelConfiguration& kernelConfig) {
		//Invoke base class init
		Kernel::Init(kernelConfig);

		//Method specific part
		MethodConfiguration config = kernelConfig.methodConfiguration;
		_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.gamma, 0.05, 0.0));
		CFL = config.CFL;
		RungeKuttaOrder = config.RungeKuttaOrder;
		
		//Allocate memory for values and residual
		spectralRadius.resize(nCellsLocal);		
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

		dt = std::min(stepInfo.NextSnapshotTime - stepInfo.Time, dt);

		return dt;
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

	//Explicit time step
	virtual void IterationStep() override {			
		//Compute residual
		ComputeResidual(values, residual, spectralRadius);

		//Determine timestep
		stepInfo.TimeStep = ComputeTimeStep(spectralRadius);

		//Compute residuals and update solution
		UpdateSolution(stepInfo.TimeStep);
		
		// Update time
		stepInfo.Time += stepInfo.TimeStep;
		stepInfo.Iteration++;
	};	

};

#endif