/*!
\file
\brief 2D Test: laminar boundary layer on the flat plate

Laminar flow for viscous gas or fluid in boundary area around flat plane.
Analitical solution for uncompressible flow was obtained by Blasius.
For validation perpose we consider subsonic flow with main parameters are:
M = 0.2, P = 10^5 ...
*/

#ifndef TurboStructured_Tests_2DTests_BlasiusFlowTest
#define TurboStructured_Tests_2DTests_BlasiusFlowTest

#include "kernel.h";


namespace BlasiusFlowTest {
	
	// struct for main parameters of the test
	struct Parameters {
		double gamma;		// specific heat ratio
		double Lx;			// domain size
		double Ly;			// domain size
		double Xplate;		// position of the plane left end
		double ro;			// init density
		double M;			// inlet Mach number
		double Pin;			// inlet pressure
		double Pout;		// outlet pressure
		double viscosity;	// dynamic viscosity
	} par;

	// Default parameters
	void DefaultSettings() {
		par.gamma = 1.4;
		par.Lx = 1.2;
		par.Ly = 0.5;
		par.Xplate = 0.2;
		par.M = 0.2;
		par.Pin = 1.0e5;
		par.Pout = par.Pin;
		par.ro = 1.2;
		par.viscosity = 1.0e-3;
	};

	// compute viscosity value
	inline double ComputeInletVelocity() {
		double sound_speed = sqrt(par.gamma * par.Pin / par.ro);
		return par.M * sound_speed;
	};

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		
		// Init config structure
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 40;
		conf.nY = 40;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.Gamma = 1.4;
		conf.IsViscousFlow = true;
		conf.Viscosity = par.viscosity;

		// BC
		conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
		conf.xLeftBoundary.Gamma = 1.4;
		conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
		conf.xRightBoundary.Gamma = 1.4;
		conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
		conf.yLeftBoundary.Gamma = 1.4;
		conf.yRightBoundary.BCType = BoundaryConditionType::MovingWall;
		conf.yRightBoundary.Gamma = 1.4;
		conf.yRightBoundary.Velocity = Vector(0, 0, 0);

		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = 10.0;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.1;
		conf.SaveSliceTime = 1.0;
		conf.ResidualOutputIterations = 20;
		
		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		kernel->Init(conf);

		// Compute driven velocity
		double Udr = ComputeInletVelocity();

		// init distributions
		NumericQuadrature Integ(3, 2);
		auto InitDriven = [Udr, &conf](Vector r) {
			// Density
			double rho = par.ro;

			// Velocity
			double u = Udr;
			double v = 0;
			double w = 0;

			// Energy
			double roe = par.Pin / (par.gamma - 1.0);

			// Compute local initial values
			std::vector<double> res(5, 0);
			res[0] = rho;
			res[1] = rho * u;
			res[2] = rho * v;
			res[3] = rho * w;
			res[4] = roe + rho * (u * u + v * v + w * w);
			return res;
		};
		kernel->SetInitialConditions(InitDriven, Integ);

		// Create slices
		kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));
		//kernel->SaveSliceToTecplot("ySlice_init.dat", kernel->slices[0]);

		//save init solution and run the test
		kernel->SaveSolution("init.dat");

		// Run test
		if (kernel->pManager->IsMaster()) std::cout << "Flat plate test runs." << std::endl <<
			"Inlet velocity: " << Udr <<  std::endl;
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
		std::cout << "Flat plate test is completed";
		std::cout << std::endl;
		return;
	};
};



#endif