#ifndef TurboStructured_Tests_2DTests_DrivenCavityTest
#define TurboStructured_Tests_2DTests_DrivenCavityTest


#include <iostream>
#include <vector>
#include "math.h"
#include "kernel.h";
#include "Methods/ExplicitRungeKuttaFVM.h"
#include "utility/Vector.h"

// Standart experiment with viscous vortex in the cavity 
//
// see Ghia et al (1982) to compare result
//
// http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem

// TO DO FIX -- tecplot save solution function wasn't works on 6 or 8 comput. nodes
namespace DrivenCavityTest {
	
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
		par.U_dr = 1.5;
		par.p = 1.0e5;
		par.rho = 1.0;
		par.comp_time = 1.0;
		par.Re = 100;
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
		conf.nX = 40;
		conf.nY = 40;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.Gamma = 1.4;
		conf.IsViscousFlow = true;
		conf.Viscosity = ComputeViscosity();
		//conf.Viscosity = 0;

		// BC
		conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
		conf.xLeftBoundary.Gamma = 1.4;
		conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
		conf.xRightBoundary.Gamma = 1.4;
		conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
		conf.yLeftBoundary.Gamma = 1.4;
		conf.yRightBoundary.BCType = BoundaryConditionType::MovingWall;
		conf.yRightBoundary.Gamma = 1.4;
		conf.yRightBoundary.Velocity = Vector(par.U_dr, 0, 0);

		// Method settings
		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.CFL = 0.5;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.DummyLayerSize = 1;

		// Computational settings
		conf.MaxTime = 100.0;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.01;
		conf.SaveSolutionIterations = 0;
		conf.SaveSliceTime = 0.01;
		conf.ResidualOutputIterations = 100;
		
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

		// init distributions
		auto init = [&](Vector r) {
			double roe = par.p / (par.gamma - 1.0);
			std::vector<double> res(5, 0);
			res[0] = par.rho;
			res[4] = roe;
			return res;
		};
		kernel->SetInitialConditions(init);

		// Create slices
		kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));
		//kernel->SaveSliceToTecplot("ySlice_init.dat", kernel->slices[0]);

		//save init solution and run the test
		kernel->SaveSolution("init.dat");

		// Run test
		if (kernel->pManager->IsMaster()) std::cout << "Lid-driven cavity test runs." << std::endl <<
			"Dynamic Viscosity: " << conf.Viscosity << "; Re = " << par.Re << std::endl;
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
		std::cout << "Lid-driven cavity test is completed";
		std::cout << std::endl;
		return;
	};
};



#endif