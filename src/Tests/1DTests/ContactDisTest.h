#ifndef TurboStructured_Tests_1DTests_ContactDisTest
#define TurboStructured_Tests_1DTests_ContactDisTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h"
#include "Methods\ExplicitRungeKuttaFVM.h"
#include "TestsUtility.h"

// single contact discontinuity test
// roL = 1.0 || roR = 0.125
// uL = 1.0  || uR = 1.0
// pL = 1.0  || pR = 1.0

namespace ContactDisTest
{
	// IC structure
	struct ShockTubeParameters {
		double roL;
		double roR;
		double P;
		double uL;
		double uR;
		double x0;
		double ref_velocity;
	};

	// Default settings and default state
	KernelConfiguration DefaultSettings() {
		KernelConfiguration conf;

		// Values by default
		conf.nDims = 1;
		conf.nX = 200;
		conf.LX = 1.0;
		conf.isPeriodicX = false;
		conf.isUniformAlongX = true;
		conf.qx = 1.00;
		conf.Gamma = TestsUtility::gamma1 + 1;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;

		conf.DummyLayerSize = 1;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;

		conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
		conf.xLeftBoundary.Gamma = 1.4;
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = 1.4;

		conf.MaxTime = 0.2;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 0;
		conf.SaveSolutionSnapshotIterations = 0;
		conf.ResidualOutputIterations = 100;

		return conf;
	};

	ShockTubeParameters DefaultState() {
		ShockTubeParameters params;
		params.roL = 1.0;
		params.roR = 0.125;
		params.P = 1.0;
		params.uL = 1.0;
		params.uR = 1.0;
		params.x0 = 0.2;		// initial discontinuity position
		params.ref_velocity = 0;

		return params;
	};

	std::valarray<double> ComputeExactSolution(ShockTubeParameters& pars, Grid& g, double time) {
		double u = 1.0 - pars.ref_velocity;
		double x_disc = pars.x0 + u * time;		// uL = uR = 1
		double p = pars.P;
		std::valarray<double> res(g.nlocalX * TestsUtility::nVar);
		for (int i = 0; i < g.nlocalX; i++) {
			double rho = pars.roR;
			if (g.CoordinateX[g.iMin + i] <= x_disc) rho = pars.roL;
			double e = p / (TestsUtility::gamma1 * rho);
			res[i * TestsUtility::nVar] = rho;
			res[i * TestsUtility::nVar + 1] = rho*u;
			res[i * TestsUtility::nVar + 2] = 0;
			res[i * TestsUtility::nVar + 3] = 0;
			res[i * TestsUtility::nVar + 4] = rho * e + 0.5 * rho * u * u;
		};
		return res;
	};

	// Collision of two media
	std::vector<double> RunSingleExperiment(int argc, char *argv[], int Nx, double MaxTime, Reconstruction RecType, RPSolver rpSolver) {
		// Result
		std::vector<double> errors;

		// Use default settings
		KernelConfiguration conf = DefaultSettings();
		ShockTubeParameters params = DefaultState();
		conf.nX = Nx;
		conf.MaxTime = MaxTime;
		conf.methodConfiguration.ReconstructionType = RecType;
		conf.methodConfiguration.RiemannProblemSolver = rpSolver;

		// Solution file
		std::ostringstream fname;

		// Init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
			fname << "PWConstant";
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
			fname << "ENO2";
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::WENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
			fname << "WENO2";
		};
		kernel->Init(conf);

		// IC
		auto initD = [&conf, &params](Vector r) {
			double p = params.P;
			double u = 0;
			double ro = 0;
			if (r.x < params.x0) {
				u = params.uL - params.ref_velocity;
				ro = params.roL;
			}
			else {
				u = params.uR - params.ref_velocity;
				ro = params.roR;
			};

			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = 0;
			res[3] = 0;
			res[4] = p / TestsUtility::gamma1 + 0.5 * ro * u * u;
			return res;
		};
		kernel->SetInitialConditions(initD);

		// Run computation
		kernel->Run();

		// Compute exact solution
		TestsUtility::exact_solution = ComputeExactSolution(params, kernel->g, conf.MaxTime);

		// Compute accuracy
		std::valarray<double> inner_values(&(kernel->values[TestsUtility::nVar * kernel->g.dummyCellLayersX]), kernel->g.nlocalX * kernel->nVariables);
		std::vector<double> L2 = TestsUtility::ComputeL2Error(inner_values, kernel->g, *(kernel->pManager));
		std::vector<double> L1 = TestsUtility::ComputeL1Error(inner_values, kernel->g, *(kernel->pManager));

		// Write in errors
		errors = L2;
		for (int i = 0; i < TestsUtility::nVar; i++) errors.push_back(L1[i]);

		// Show the errors
		if (kernel->pManager->IsMaster()) {
			std::cout << "Contact discontinuity test was finished. ";
			std::cout << std::endl;
			std::cout << "ReconstructionType: " << fname.str();
			std::cout << ", Nx = " << Nx << std::endl;
			std::cout << "L2_rho = " << L2[0];
			std::cout << ", L1_rho = " << L1[0];
			std::cout << std::endl;
		};
		
		// Save both solution in TecPlot
		fname << ", Nx=";
		fname << Nx;
		fname << ", t=";
		fname << MaxTime;
		fname << ".dat";
		std::string filename = fname.str();
		kernel->SaveSolution(filename);
		TestsUtility::SaveExactSolution(filename, kernel->g, *(kernel->pManager));

		// Save the errors
		std::string fname2("Test-ContactDiscontinuity.dat");
		TestsUtility::WriteErrors(fname2, errors, kernel->g, *(kernel->pManager));

		// Finalize kernel
		kernel->Finalize();

		return errors;
	};

	void RunExperiment(int argc, char *argv[]) {
		int Nx = 400;
		double MaxTime = 0.2;

		// Reconstruction type
		//Reconstruction RecType{ Reconstruction::PiecewiseConstant };
		//Reconstruction RecType{ Reconstruction::ENO2PointsStencil };
		Reconstruction RecType{ Reconstruction::WENO2PointsStencil };

		// RP solver
		RPSolver rSolver{ RPSolver::GodunovSolver };

		// collect all errors in file
		std::vector<double> err;
		err = RunSingleExperiment(argc, argv, Nx, MaxTime, RecType, rSolver);

		return;
	};

}	//end of namespace area


#endif
