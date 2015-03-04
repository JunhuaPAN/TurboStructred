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
	//Inherit constructor
	//using Kernel::Kernel;
	ExplicitRungeKuttaFVM(int* argc, char **argv[]) : Kernel(argc, argv) {};

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
		for (double& r : residual) r = 0; //Nullify residual

		// array of pointers to cell values
		std::vector<double*> U(2);
		
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

				//Initial load of values
				U[0] = getCellValues(iMin - 1, j, k);

				for (int i = iMin; i <= iMax + 1; i++) {
					//Update stencil values
					U[1] = getCellValues(i, j, k);

					//Compute flux
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
					fr = result.Fluxes;
	
					//Update residuals
					if(i > iMin)
					{
						//Fluxes difference equals residual
						int idx = getSerialIndexLocal(i - 1, j, k);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

						//Add up spectral radius estimate
						spectralRadius[idx] += fS * result.MaxEigenvalue;
					};
          
					//Shift stencil
					fl = fr;
					U[0] = U[1];
				 };
			};
		};
		
		// J step
		if (nDims > 1)
		{
			fn = Vector(0.0, 1.0, 0.0); // y direction
			for (int k = kMin; k <= kMax; k++) {
				for (int i = iMin; i <= iMax; i++) {
					//Compute face square
					double fS = 0;
					if (nDims == 2) fS = hx;
					if (nDims == 3) fS = hx * hz;

					//Initial load of values
					U[0] = getCellValues(i, jMin - 1, k);

					for (int j = jMin; j <= jMax + 1; j++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;
	
						//Update residuals
						if(j > jMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j - 1, k);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
						};
          
						//Shift stencil
						fl = fr;
						U[0] = U[1];
					};
				};
			};
		};
		
		// K step
		if (nDims > 2)
		{
			fn = Vector(0.0, 0.0, 1.0); // z direction
			for (int i = iMin; i <= iMax; i++) {
				for (int j = jMin; j <= jMax; j++) {
					//Face square
					double fS = hx * hy;

					//Initial load of values
					U[0] = getCellValues(i, j, kMin - 1);

					for (int k = kMin; k <= kMax + 1; k++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;
	
						//Update residuals
						if(k > kMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j, k - 1);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
						};
          
						//Shift stencil
						fl = fr;
						U[0] = U[1];
					};
				};
			};
		};


		//Right hand side
		double mu = viscosity;
		double sigma = 0.0001;
		double volume = hx * hy * hz;		
		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);

					//Cell values
					double *V = getCellValues(i, j, k);
					double *VxL = getCellValues(i - 1, j, k);
					double *VxR = getCellValues(i + 1, j, k);
					double *VyL = getCellValues(i, j - 1, k);
					double *VyR = getCellValues(i, j + 1, k);
					//double *vzL = getCellValues(i, j, k - 1);
					//double *vzR = getCellValues(i, j, k + 1);

					
					double u = V[1]/V[0];
					double uxL = VxL[1] / VxL[0];
					double uxR = VxR[1] / VxR[0];
					double dudx = (uxR - uxL) / (2 * hx);
					double d2udx2 = (uxR - 2*u + uxL) / (hx * hx);

					double v = V[2]/V[0];
					double vyL = VyL[1] / VyL[0];
					double vyR = VyR[1] / VyR[0];
					double dvdy = (vyR - vyL) / (2 * hy);
					double d2vdy2 = (vyR - 2*v + vyL) / (hy * hy);

					Vector velocity = Vector(u, v, V[3]/V[0]); //Cell velocity

					//Friction forces

					//Mass forces
					double mForce = 0;

					//Potential forces
					Vector pForce = Vector(sigma, 0, 0);

					//Compute total residual
					residual[idx * nVariables];			//ro
					residual[idx * nVariables + 1] += hx * pForce.x;	//rou
					residual[idx * nVariables + 2] += hy * pForce.y;	//rov
					residual[idx * nVariables + 3] += hz * pForce.z;		//row
					residual[idx * nVariables + 4] += pForce * velocity;	//roE
				};
			};
		};
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