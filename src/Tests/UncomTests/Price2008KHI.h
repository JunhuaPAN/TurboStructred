#ifndef TurboStructured_Tests_UncomTests_Price2008KHI
#define TurboStructured_Tests_UncomTests_Price2008KHI


#include <iostream>
#include <vector>
#include "math.h"
#include "kernel.h";
#include "Methods/ExplicitRungeKuttaFVM.h"
#include "utility/Vector.h"

// 2D Kelvin-Helmholtz instability from [Daniel J. Price - 2008]
//
// http://arxiv.org/pdf/0709.2772v3.pdf

namespace Price2008KHI {
	double PI = 3.14159265359;
	
	// struct for main parameters of the test (default settings)
	struct Parameters {
		double gamma{ 5.0 / 3.0 };	// specific heat ratio
		double Lx{ 1.0 };			// domain size
		double Ly{ 0.5 };			// domain size
		double p{ 2.5 };			// initial pressure
		double rho_top{ 2.0 };		// initial densities
		double rho_bot{ 1.0 };		// time of computation
		double Vsh{ 0.5 };			// shear velocity of up and down part
		double lambda{ 1.0 / 6.0 };	// perturbation wave length
		double ampl{ 0.025 };		// perturbation amplitude
	} par;

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		
		// Create config file
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 160;
		conf.nY = 160;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = true;
		conf.isPeriodicY = false;
		conf.Gamma = par.gamma;
		conf.MolarMass = 1.0;

		// BC
		conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
		conf.yLeftBoundary.Gamma = 1.4;
		conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
		conf.yRightBoundary.Gamma = 1.4;

		// Method configuration
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.4;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = 20.0;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 0.1;
		conf.ResidualOutputIterations = 50;

		// Init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		kernel->Init(conf);


		NumericQuadrature Integ(8, 2);
		auto Init = [&conf](Vector r) {
			// initial perturbation
			double v = par.ampl * sin(2 * PI * r.x / par.lambda);
			double w = 0;
			double u = par.Vsh;
			double rho = par.rho_top;
			if (r.y < 0.5 * conf.LY) {
				u = -par.Vsh;
				rho = par.rho_bot;
			};

			double rhoe = par.p / (conf.Gamma - 1);
			std::vector<double> res(5);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = rhoe + 0.5 * rho * (u * u + v * v + w * w);
			return res;
		};
		kernel->SetInitialConditions(Init, Integ);

		//save solution
		kernel->SaveSolution("init.dat");

		//run computation
		kernel->Run();

		//finalize kernel
		kernel->Finalize();
	};

	// Run Computation Experiment
	void RunExperiment(int argc, char *argv[]) {
		// here you can change any parameters

		// Run experiments
		RunSingleExperiment(argc, argv);

		//end of experiments
		std::cout << "Lid-driven cavity test is completed";
		std::cout << std::endl;
		return;
	};
};



#endif