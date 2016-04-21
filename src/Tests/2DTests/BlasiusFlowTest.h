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
		par.Lx = 1.25;
		par.Ly = 0.5;
		par.Xplate = 0.25;
		par.M = 0.2;
		par.Pin = 1.0e5;
		par.Pout = par.Pin;
		par.ro = 1.2;
		par.viscosity = 1.0e-4;
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
		conf.Gamma = par.gamma;
		conf.IsViscousFlow = true;
		conf.Viscosity = par.viscosity;

		// Boundary conditions Subsonic inlet
		conf.xLeftBoundary.BCType = BoundaryConditionType::SubsonicInlet;
		conf.xLeftBoundary.Gamma = par.gamma;
		conf.xLeftBoundary.Density = par.ro;
		conf.xLeftBoundary.Pstatic = par.Pin;
		conf.xLeftBoundary.Vdirection = Vector(1, 0, 0);

		// Supersonic outlet
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = par.gamma;

		// No-slip condition on the bottom
		conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
		conf.yLeftBoundary.Gamma = par.gamma;

		// Symmetry or somthing else on the top border
		conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.yRightBoundary.Gamma = par.gamma;

		// Symmetry in front of the plate
		// Compute driven velocity
		double Udr = ComputeInletVelocity();
		conf.yLeftSpecialBoundary.BCType = BoundaryConditionType::MovingWall;
		conf.yLeftSpecialBoundary.Gamma = par.gamma;
		conf.yLeftSpecialBoundary.Velocity = Vector(Udr, 0, 0);

		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = 10.0;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.1;
		conf.SaveSliceTime = 0.1;
		conf.ResidualOutputIterations = 100;
		
		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		kernel->Init(conf);

		// init distributions
		NumericQuadrature Integ(3, 2);
		auto InitDriven = [Udr](Vector r) {
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
			res[4] = roe + 0.5 * rho * (u * u + v * v + w * w);
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