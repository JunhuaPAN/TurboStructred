#ifndef TurboStructured_Tests_2DTests_ExactEulerSolution
#define TurboStructured_Tests_2DTests_ExactEulerSolution

#include "kernel.h";

// Test with smooth solution of 2D Euler system of equations
//
// from RICHARD LISKA and BURTON WENDROFF 2002 TO DO ACCOMPLISH
//
// http://math.unm.edu/~bbw/sisc.pdf

namespace ExactEulerSolution {
	double PI = 3.14159265359;

	// struct for main parameters of the test
	struct Parameters {
		double gamma;		// specific heat ratio
		double amplitude;	// amplitude of density perturbations: a*sin (pi * ((x + y) - t * (u + v)) )
		double Lx;			// domain size
		double Ly;			// domain size
		double u;			// velocity x-component
		double v;			// velocity y-component
		double p;			// pressure is constant
		double rho;			// constant part of density
		double comp_time;	// time of computation
	} par;

	// Default parameters
	void DefaultSettings() {
		par.gamma = 1.4;
		par.amplitude = 0.2;
		par.Lx = 1;
		par.Ly = 1;
		par.u = 1;
		par.v = -0.5;
		par.p = 1.0;
		par.rho = 1.0;
	};

	// compute total density at the point r 		
	inline double ComputeDensity(Vector r, double t) {
		return par.rho + par.amplitude * sin(PI * (r.x + r.y - t * (par.u + par.v) ) );
	};

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 100;
		conf.nY = 100;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = true;
		conf.isPeriodicY = true;
		conf.Gamma = par.gamma;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.methodConfiguration.CFL = 0.4;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.DummyLayerSize = 1;

		conf.MaxTime = par.comp_time;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionTime = 0.1 * par.comp_time;
		conf.ResidualOutputIterations = 20;

		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::WENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2CharactVars) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2CharactVars>(&argc, &argv));
		};
		kernel->Init(conf);

		auto initD = [&conf](Vector r) {
			double rho = ComputeDensity(r, 0);
			double ro;
			double p;
			double u;

			// compute ro_e and write conservative variables
			double roe = p / (par.gamma - 1.0);
			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = 0.0;
			res[3] = 0.0;
			res[4] = roe + 0.5 * ro * u * u;		//total energy equals internal one because a motion is absent

			return res;
		};
		kernel->SetInitialConditions(initD);

		//save solution
		kernel->SaveSolution("init.dat");

		// Set sensors if needed
		auto GetInEnergy = [](std::valarray<double> vals) {
			double roe = vals[4] - 0.5 * (vals[1] * vals[1] + vals[2] * vals[2]) / vals[0];
			return roe / vals[0];
		};

		kernel->isSensorEnable = true;
		kernel->SaveSensorRecordIterations = 1;

		//run computation
		kernel->Run();

		//finalize kernel
		kernel->Finalize();
	};

	// Run Computation Experiment
	void RunExperiment(int argc, char *argv[]) {
		// Fill parameters structure (SI system)
		DefaultSettings();

		// Run experiments
		RunSingleExperiment(argc, argv);

		//end of experiments
		std::cout << "Smooth Solution Test is completed";
		std::cout << std::endl;
		return;
	};
};



#endif