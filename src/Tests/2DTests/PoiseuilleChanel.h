/*!
\file
\brief 2D Test: Poiseuille flow beatween two walls 

Laminar flow for viscous gas or fluid in boundary area around flat plane.
Analitical solution for uncompressible flow was obtained by Blasius.
For validation perpose we consider subsonic flow with main parameters are:
M = 0.2, P = 10^5 ...
*/

#ifndef TurboStructured_Tests_2DTests_PoiseuilleChanelTest
#define TurboStructured_Tests_2DTests_PoiseuilleChanelTest

#include "kernel.h";


namespace PoiseuilleChanel {
	
	// struct for main parameters of the test
	struct Parameters {
		double gamma;
		double Lx;			// domain size
		double Ly;			// domain size
		double ro;			// init density
		double Pin;			// inlet pressure
		double sigma;		// (Pin - Pout) / Lx 
		double viscosity;	// dynamic viscosity
	} par;

	// Default parameters
	void DefaultSettings() {
		par.gamma = 1.4;
		par.Lx = 1.0;
		par.Ly = 0.5;
		par.Pin = 101579;
		par.sigma = 2.0;
		par.ro = 1.0;
		par.viscosity = 0.01;
	};

	// compute velocity at the center of the channel
	double ComputeCenterVelocity() {
		return 0.125 * par.sigma * par.Ly * par.Ly / par.viscosity;
	};

	// compute Re number
	double ComputeRe() {
		auto U = ComputeCenterVelocity();
		return par.ro * U * par.Ly / par.viscosity;
	};

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		
		// Init config structure
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 4;
		conf.nY = 80;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicY = false;
		conf.Gamma = par.gamma;
		conf.IsViscousFlow = true;		
		conf.Viscosity = par.viscosity;
		conf.Sigma = Vector(par.sigma, 0, 0);
		conf.IsExternalForceRequared = true;

		// Discribe grid compression
		BlockNode y_bot, y_cen;

		// Y direction first
		y_bot.N_cells = 0.5 * conf.nY;
		y_bot.q_com = 1.05;
		conf.CompressionY[0] = y_bot;

		y_cen.pos = 0.5 * par.Ly;
		y_cen.q_com = 1.0 / y_bot.q_com;
		y_cen.N_cells = conf.nY - y_bot.N_cells;
		conf.CompressionY.push_back(y_cen);

		// Discribe boundary conditions
		BoundaryConditionConfiguration Wall(BoundaryConditionType::Wall);
		conf.MyConditions[1] = Wall;

		// Describe conditions on the domain sides
		conf.yLeftBoundary.SetMarker(1);
		conf.yRightBoundary.SetMarker(1);
		
		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
		conf.DummyLayerSize = 1;

		// compute velocity at center of the channel for laminar regime
		auto Uc = ComputeCenterVelocity();
		auto Re = ComputeRe();

		// Computational settings
		conf.MaxTime = 10 * par.Lx / Uc;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.1;
		conf.SaveSliceTime = 0.05;
		conf.SaveBinarySolIterations = 10000;
		conf.ResidualOutputIterations = 500;
		
		// init kernel
		auto kernel = CreateKernel(conf, argc, argv);
		kernel->Init(conf);

		// init distributions
		NumericQuadrature Integ(5, 2);
		auto InitTrivial = [](Vector r) {
			// Compute primitive values
			double rho = par.ro;
			double u = 0;
			double v = 0;
			double w = 0;
			double roe = par.Pin / (par.gamma - 1.0);

			// Compute local initial values
			std::vector<double> res(5, 0);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = roe + 0.5 * rho * (u * u + v * v + w * w);
			return res;
		};
		auto InitExact = [Uc](Vector r) {
			auto y = 1.0 - abs( 2 * r.y / par.Ly - 1.0);
			auto vel = Uc * y * (2.0 - y);

			// Compute primitive values
			double rho = par.ro;
			double u = vel;
			double v = 0;
			double w = 0;
			double roe = par.Pin / (par.gamma - 1.0);


			// Compute local initial values
			std::vector<double> res(5, 0);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = roe + 0.5 * rho * (u * u + v * v + w * w);
			return res;
		};
		auto InitHigher = [Uc](Vector r) {
			auto y = 1.0 - abs(2 * r.y / par.Ly - 1.0);
			auto vel = 1.2 * Uc * y * (2.0 - y);

			// Compute primitive values
			double rho = par.ro;
			double u = vel;
			double v = 0;
			double w = 0;
			double roe = par.Pin / (par.gamma - 1.0);


			// Compute local initial values
			std::vector<double> res(5, 0);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = roe + 0.5 * rho * (u * u + v * v + w * w);
			return res;
		};
		auto InitLower = [Uc](Vector r) {
			auto y = 1.0 - abs(2 * r.y / par.Ly - 1.0);
			auto vel = 0.9 * Uc * y * (2.0 - y);

			// Compute primitive values
			double rho = par.ro;
			double u = vel;
			double v = 0;
			double w = 0;
			double roe = par.Pin / (par.gamma - 1.0);


			// Compute local initial values
			std::vector<double> res(5, 0);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = roe + 0.5 * rho * (u * u + v * v + w * w);
			return res;
		};
		kernel->SetInitialConditions(InitLower, Integ);

		// Create slices
		kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));
		kernel->SaveSliceToTecplot("slice_init.dat", kernel->slices[0]);

		//save init solution and run the test
		kernel->SaveSolutionToTecplot("init.dat");

		// Run test
		if (kernel->pManager->IsMaster()) std::cout << "Poiseuille test runs." << std::endl <<
			"U_center: " << Uc << ", Re: " << Re << std::endl;
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
		std::cout << "Poiseuille test is completed";
		std::cout << std::endl;
		return;
	};
};

#endif