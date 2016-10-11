#ifndef TurboStructured_Tests_2DTests_ExactEulerSolutionTest
#define TurboStructured_Tests_2DTests_ExactEulerSolutionTest

#include "kernel.h";

// Test with smooth solution of 2D Euler system of equations
//
// from RICHARD LISKA and BURTON WENDROFF 2002
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
		par.comp_time = 2.0;
	};

	// compute total density at the point r 		
	inline double ComputeDensity(Vector r, double t) {
		return par.rho + par.amplitude * sin(2.0 * PI * (r.x + r.y - t * (par.u + par.v) ) );
	};

	// apply compression to the center with coeff q
	void ApplyCompression(KernelConfiguration& conf, double q) {
		BlockNode nleft, ncenter;

		// X direction first
		nleft.N_cells = conf.nX / 2;
		nleft.q_com = q;
		ncenter.pos = 0.5 * conf.LX;
		ncenter.N_cells = conf.nX - nleft.N_cells;
		ncenter.q_com = 1.0 / nleft.q_com;
		conf.CompressionX[0] = nleft;
		conf.CompressionX.push_back(ncenter);

		// Y direction now
		nleft.N_cells = conf.nY / 2;
		nleft.q_com = q;
		ncenter.pos = 0.5 * conf.LY;
		ncenter.N_cells = conf.nY - nleft.N_cells;
		ncenter.q_com = 1.0 / nleft.q_com;
		conf.CompressionY[0] = nleft;
		conf.CompressionY.push_back(ncenter);

		return;
	}

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 40;
		conf.nY = 40;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = true;
		conf.isPeriodicY = true;
		conf.Gamma = par.gamma;

		// apply grid compression to the center
		ApplyCompression(conf, 0.95);

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2PointsStencil;
		conf.methodConfiguration.CFL = 0.4;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.DummyLayerSize = 1;

		conf.MaxTime = par.comp_time;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionTime = 0.2 * par.comp_time;
		conf.SaveSliceTime = par.comp_time;
		conf.ResidualOutputIters = 20;

		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2PointsStencil) {
				kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<Linear2PointsStencil>(&argc, &argv));
		};
		kernel->Init(conf);

		NumericQuadrature In(5, 2);
		auto initD = [&conf](Vector r) {
			double ro = ComputeDensity(r, 0);
			double p = par.p;
			double u = par.u;
			double v = par.v;

			// compute ro_e and write conservative variables
			double roe = p / (par.gamma - 1.0);
			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = ro * v;
			res[3] = 0.0;
			res[4] = roe + 0.5 * ro * (u * u + v * v);		//total energy equals internal one because a motion is absent

			return res;
		};
		kernel->SetInitialConditions(initD, In);

		// Add slice
		kernel->slices.push_back(Slice(-1, conf.nY / 2, 0));

		//save solution
		kernel->SaveSolutionToTecplot("init.dat");
		kernel->SaveSliceToTecplot("slice_init.dat", kernel->slices[0]);

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
		std::cout << "Smooth Solution Test is complited";
		std::cout << std::endl;
		return;
	};
};



#endif