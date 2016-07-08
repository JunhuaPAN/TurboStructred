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
		par.Pin = 101579;
		par.Pout = par.Pin;
		par.ro = par.Pin * par.gamma / (1006.43 * 300.214 * (par.gamma - 1));
		par.viscosity = 1.7894e-03;
		par.M = 10.0 / sqrt(par.gamma * par.Pin / par.ro);		// Udriven = 10
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
		conf.nX = 60;
		conf.nY = 30;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.Gamma = par.gamma;
		conf.IsViscousFlow = true;
		conf.Viscosity = par.viscosity;

		// Compute driven velocity
		double Udr = ComputeInletVelocity();

		// Discribe Boundary conditions
		// Subsonic inlet
		BoundaryConditionConfiguration Inlet(BoundaryConditionType::SubsonicInlet);
		Inlet.Density = par.ro;
		Inlet.Pstatic = par.Pin;
		Inlet.Velocity = Vector(Udr, 0, 0);
		// Supersonic outlet
		BoundaryConditionConfiguration Outlet(BoundaryConditionType::Natural);
		// Plate condition
		BoundaryConditionConfiguration Wall(BoundaryConditionType::Wall);
		// Symmetry sides
		BoundaryConditionConfiguration Symmetry(BoundaryConditionType::Symmetry);

		// Create mapping Marker -> BCondition
		conf.MyConditions[1] = Inlet;
		conf.MyConditions[2] = Outlet;
		conf.MyConditions[3] = Wall;
		conf.MyConditions[4] = Symmetry;

		// Describe conditions on the domain sides
		conf.xLeftBoundary.SetMarker(1);
		conf.xRightBoundary.SetMarker(2);
		conf.yRightBoundary.SetMarker(2);
		// complex case
		conf.yLeftBoundary.isComplex = true;
		auto getMarkerBottom = [](Vector r) {
			if (r.x < par.Xplate) return 4;
			return 3;
		};
		conf.yLeftBoundary.getMarker = getMarkerBottom;
		
		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = 10 * par.Lx / Udr;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.02;
		conf.SaveSliceTime = par.Lx / Udr;
		conf.ResidualOutputIterations = 100;
		
		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2psLim) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limBarsJespersen> >(&argc, &argv));
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
		auto InitLinearBL = [Udr](Vector r) {
			// Compute laminar layer height
			double h = 4.91 * sqrt( par.viscosity * (r.x - par.Xplate) / (Udr * par.ro) );

			// Density
			double rho = par.ro;

			// Velocity
			double u = Udr;
			if (r.y < h) u *= r.y / h;
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
		kernel->slices.push_back(Slice((int)(0.9 * conf.nX), -1, 0));
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