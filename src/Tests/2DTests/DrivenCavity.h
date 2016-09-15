/*!
\file
\brief 2D test: driven cavity problem

Standart experiment with viscous vortex in the cavity see Ghia et al (1982) to compare result
http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem
*/

#ifndef TurboStructured_Tests_2DTests_DrivenCavityTest
#define TurboStructured_Tests_2DTests_DrivenCavityTest

#include "kernel.h"


namespace DrivenCavity {
	
	// struct for main parameters of the test
	struct Parameters {
		double gamma;		// specific heat ratio
		double Lx;			// domain size
		double Ly;			// domain size
		double U_dr;		// driven velocity
		double p;			// initial pressure
		double rho;			// initial density
		double comp_time;	// time of computation
		double Re;			// Reinolds number
	} par;

	// Default parameters
	void DefaultSettings() {
		par.gamma = 1.4;
		par.Lx = 1.0;
		par.Ly = 1.0;
		par.U_dr = 1.0;
		par.p = 1.0e5;
		par.rho = 1.0;
		par.comp_time = 10.0 * par.U_dr;
		par.Re = 400;
	};

	// compute viscosity value
	inline double ComputeViscosity() {
		return par.U_dr * par.rho * par.Lx / par.Re;
	};

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int argc, char *argv[]) {
		
		// Init config structure
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 80;
		conf.nY = 80;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.Gamma = par.gamma;
		conf.IsViscousFlow = true;
		conf.Viscosity = ComputeViscosity();

		// Discribe grid compression
		double qx = 1.08;
		double qy = 1.05;

		BlockNode xl, xc, yl, yc;
		// X
		xl.N_cells = conf.nX / 2;
		xl.q_com = qx;
		xc.N_cells = conf.nX - xl.N_cells;
		xc.q_com = 1.0 / qx;
		xc.pos = 0.5 * conf.LX;
		// Y
		yl.N_cells = conf.nY / 2;
		yl.q_com = qy;
		yc.N_cells = conf.nY - yl.N_cells;
		yc.q_com = 1.0 / qy;
		yc.pos = 0.5 * conf.LY;
		// write in config
		conf.CompressionX[0] = xl;
		conf.CompressionX.push_back(xc);
		conf.CompressionY[0] = yl;
		conf.CompressionY.push_back(yc);

		// BC
		BoundaryConditionConfiguration TopPlate(BoundaryConditionType::MovingWall);
		TopPlate.Velocity = Vector(par.U_dr, 0, 0);
		conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Wall);
		conf.MyConditions[2] = TopPlate;
		conf.xLeftBoundary.SetMarker(1);
		conf.xRightBoundary.SetMarker(1);
		conf.yLeftBoundary.SetMarker(1);
		conf.yRightBoundary.SetMarker(2);

		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = par.comp_time;
		conf.MaxIteration = 30000000;
		conf.SaveSolutionTime = 0.1;
		conf.SaveSliceTime = 0.1;
		conf.SaveSolutionIterations = 100000;
		conf.SaveBinarySolIterations = 100000;
		conf.ResidualOutputIterations = 1000;
		
		// init kernel
		auto kernel = CreateKernel(conf, argc, argv);
		kernel->Init(conf);

		// init distributions
		NumericQuadrature Int(1, 2);
		auto init = [&](Vector r) {
			double roe = par.p / (par.gamma - 1.0);
			std::vector<double> res(5, 0);
			res[0] = par.rho;
			res[4] = roe;
			return res;
		};
		kernel->SetInitialConditions(init, Int);

		// Create slices
		kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));
		//kernel->SaveSliceToTecplot("ySlice_init.dat", kernel->slices[0]);

		//save init solution and run the test
		kernel->LoadSolution("sol_to_load.sol");
		kernel->SaveSolutionToTecplot("init.dat");

		// Run test
		if (kernel->pManager->IsMaster()) std::cout << "Lid-driven cavity test runs." << std::endl <<
			"Dynamic Viscosity: " << conf.Viscosity << "; Re = " << par.Re << std::endl;
		kernel->Run();

		// Save solution
		kernel->SaveSolution("final_sol.sol");

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
		std::cout << "Lid-driven cavity test is completed";
		std::cout << std::endl;
		return;
	};
};



#endif