#ifndef TurboStructured_Tests_Tests
#define TurboStructured_Tests_Tests

#include <iostream>
#include <vector>
#include <random>


#include "kernel/kernel.h"
#include "Methods/ExplicitRungeKuttaFVM.h"
#include "Tests/1DTests/testlist.h"
#include "Tests/2DTests/testlist.h"
#include "Tests/UncomTests/testlist.h"
#include "RiemannSolvers/RiemannSolversList.h"

//// Tests 1D
void RunSODTestRoe1D(int argc, char *argv[]) {
	
	// Make configuration file for SOD test
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.LX = 1.0;
	conf.nX = 400;
	conf.isPeriodicX = false;
	conf.DummyLayerSize = 1;

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.xLeftBoundary.SetMarker(1);
	conf.xRightBoundary.SetMarker(1);

	// Model settings
	conf.Gamma = 1.4;
	conf.IsViscousFlow = false;

	// Method settings
	conf.methodConfiguration.CFL = 0.3;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.0;
	conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;

	// Task settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.25;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 100;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2psLim) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limVenkatar> >(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial conditions
	struct ShockTubeParameters {
		double roL{ 1.0 };
		double PL{ 1.0 };
		double uL{ 0 };
		double roR{ 0.125 };
		double PR{ 0.1 };
		double uR{ 0 };
		double x0{ 0.5 };
	} params;
	auto initD = [&params, &conf](Vector r) {
		double ro, u, p;
		if (r.x < params.x0) {
			ro = params.roL;
			u = params.uL;
			p = params.PL;
		}
		else {
			ro = params.roR;
			u = params.uR;
			p = params.PR;
		};

		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = p / (conf.Gamma - 1.0) + 0.5 * ro * u * u;
		return res;
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

void RunContactDisconTest1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 400;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.Gamma = 1.4;

	// Method settings
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.45;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.xLeftBoundary.SetMarker(1);
	conf.xRightBoundary.SetMarker(1);

	// Computation parameters
	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 50;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2psLim) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limBarsJespersen> >(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limVenkatar> >(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial conditions
	struct ShockTubeParameters {
		double roL{ 1.0 };
		double PL{ 1.0 };
		double uL{ 1.0 };
		double roR{ 0.125 };
		double PR{ 1.0 };
		double uR{ 1.0 };
		double x0{ 0.2 };
	} params;

	auto initD = [&params, &conf](Vector r) {
		double ro, u, p;
		if (r.x < params.x0) {
			ro = params.roL;
			u = params.uL;
			p = params.PL;
		}
		else {
			ro = params.roR;
			u = params.uR;
			p = params.PR;
		};

		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = p / (conf.Gamma - 1.0) + 0.5 * ro * u * u;
		return res;
	};
	kernel->SetInitialConditions(initD);

	// Save solution
	kernel->SaveSolution("init.dat");

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};


/*


// Y direction test
void RunSODYTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 1;
	conf.nY = 400;
	conf.LX = 1.0;
	conf.LY = 2.0;
	conf.isPeriodicX = true;
	conf.isUniformAlongX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = true;

	// BC
	conf.yLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.yRightBoundary.Gamma = 1.4;

	// Method parameters
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::GodunovSolver;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.Gamma = 1.4;

	// Task and output settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	//conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.1;
	conf.ResidualOutputIterations = 10;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	struct ShockTubeParameters {
		double gammaL = 1.4;
		double gammaR = 1.4;
		double roL = 1.0;
		double PL = 1.0;
		Vector uL = { 0, 0, 0 };
		double roR = 0.125;
		double PR = 0.1;
		Vector uR = { 0, 0, 0 };
	} params;

	double y0 = 0.5 * conf.LY;
	auto init = [&params, &conf, y0](Vector r) {
		std::vector<double> res(5);
		double rho, p;
		double gamma = params.gammaL;
		Vector V;
		if (r.y < y0) {
			rho = params.roL;
			p = params.PL;
			V = params.uL;
		}
		else {
			rho = params.roR;
			p = params.PR;
			V = params.uR;
		};

		res[0] = rho;
		res[1] = rho * V.x;
		res[2] = rho * V.y;
		res[3] = rho * V.z;
		res[4] = p / (gamma - 1) + 0.5 * rho * (V.x * V.x + V.y * V.y + V.z * V.z);

		return res;
	};

	// IC
	kernel->SetInitialConditions(init);

	// Init the slice and save initial solution
	//kernel->SaveSolution("init.dat");
	kernel->slices.push_back(Slice(1, -1, 0));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

void RunSODInverseYTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 1;
	conf.nY = 400;
	conf.LX = 1.0;
	conf.LY = 2.0;
	conf.isPeriodicX = true;
	conf.isUniformAlongX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = true;

	// BC
	conf.yLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.yRightBoundary.Gamma = 1.4;

	// Method parameters
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::GodunovSolver;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.Gamma = 1.4;

	// Task and output settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	//conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.1;
	conf.ResidualOutputIterations = 10;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	struct ShockTubeParameters {
		double gammaL = 1.4;
		double gammaR = 1.4;
		double roL = 1.0;
		double PL = 1.0;
		Vector uL = { 0, 0, 0 };
		double roR = 0.125;
		double PR = 0.1;
		Vector uR = { 0, 0, 0 };
	} params;

	double y0 = 0.5 * conf.LY;
	auto init = [&params, &conf, y0](Vector r) {
		std::vector<double> res(5);
		double rho, p;
		double gamma = params.gammaL;
		Vector V;
		if (r.y < y0) {
			rho = params.roL;
			p = params.PL;
			V = params.uL;
		}
		else {
			rho = params.roR;
			p = params.PR;
			V = params.uR;
		};

		res[0] = rho;
		res[1] = rho * V.x;
		res[2] = rho * V.y;
		res[3] = rho * V.z;
		res[4] = p / (gamma - 1) + 0.5 * rho * (V.x * V.x + V.y * V.y + V.z * V.z);

		return res;
	};

	// IC
	kernel->SetInitialConditions(init);

	// Init the slice and save initial solution
	//kernel->SaveSolution("init.dat");
	kernel->slices.push_back(Slice(1, -1, 0));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};


//// Multidimension cases
void Run2DComparisonTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 2;
	conf.nY = 2;
	conf.LX = 1.0;
	conf.LY = 2.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	conf.Gamma = 1.4;

	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 10;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSolutionIterations = 1;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	auto initD = [&conf](Vector r) {
		double ro;
		double p;
		double v;

		// domain divided into four areas
		if (r.x < 0.5 * conf.LX) {
			if (r.y < 0.5 * conf.LY) {
				ro = 1.0;
				p = 1.0;
				v = 0.5;
			}
			else {
				ro = 2.0;
				p = 1.0;
				v = 1.0;
			};
		}
		else if (r.y < 0.5 * conf.LY) {
			ro = 1.0;
			p = 2.0;
			v = 0.5;
		}
		else {
			ro = 2.0;
			p = 2.0;
			v = 1.0;
		};

		// compute ro_e and write conservative variables
		double roe = p / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = 0.0;
		res[2] = ro * v;
		res[3] = 0.0;
		res[4] = roe + 0.5 * ro * v * v;		//total energy equals internal one because a motion is absent

		return res;
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();

};

//test for comparison of fluxes in different codes
void RunFluxesTest2D(int argc, char *argv[]) {
	double viscosity = 2.0;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 400;
	conf.nY = 4;
	conf.LX = 2.0;
	conf.LY = 1.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0;
	conf.SaveSolutionIterations = 1;
	conf.ResidualOutputIterations = 1;

	conf.Viscosity = viscosity;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [](Vector r) {
		double u = r.x + r.y;
		double v = r.x - r.y;
		double w = 0;
		double ro = 1.5*(1.0 + r.y*0.01);
		double roe = ro * 100;
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro*u;
		res[2] = ro*v;
		res[3] = ro*w;
		res[4] = roe + 0.5*ro*(u*u + v*v);
		return res;
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// 2D Shock Tube test http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
void RunSODTestRoe2D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 50;
	conf.nY = 50;
	conf.LX = 1.0;
	conf.LY = 1.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::SymmetryX;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = { 0, 0, 0 };
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = { 0, 0, 0 };
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.4, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// 2D Noh Problem [Boscheri - 2015]
void RunNohProblem2D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 100;
	conf.nY = 100;
	conf.LX = 2.0;
	conf.LY = 2.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::SymmetryX;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 0.6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 20;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// Test parameters
	double ro1 = 1.0;
	double ro2 = 1.0;
	double ro3 = 0.125;
	double p1 = 1;
	double p2 = 0.1;
	double p3 = 0.1;

	auto initD = [ro1, ro2, ro3, p1, p2, p3, &conf](Vector r) {
		double ro;
		double p;

		// domain divided into three areas
		if (r.x < conf.LX / 7.0) {
			ro = ro1;
			p = p1;
		}
		else if (r.y < 0.5 * conf.LY) {
			ro = ro2;
			p = p2;
		}
		else {
			ro = ro3;
			p = p3;
		};

		// compute ro_e and write conservative variables
		double roe = p / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = 0.0;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] = roe;		//total energy equals internal one because a motion is absent

		return res;
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// 2D Mono-Material Triple Point Problem from [Boscheri - 2015]
void RunTriplePointRoe2D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 200;
	conf.nY = 150;
	conf.LX = 7.0;
	conf.LY = 3.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 5.5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.5;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// Test parameters
	double ro1 = 1.0;
	double ro2 = 1.0;
	double ro3 = 0.125;
	double p1 = 1;
	double p2 = 0.1;
	double p3 = 0.1;

	auto initD = [ro1, ro2, ro3, p1, p2, p3, &conf](Vector r) {
		double ro;
		double p;

		// domain divided into three areas
		if (r.x < conf.LX / 7.0) {
			ro = ro1;
			p = p1;
		}
		else if (r.y < 0.5 * conf.LY) {
			ro = ro2;
			p = p2;
		}
		else {
			ro = ro3;
			p = p3;
		};

		// compute ro_e and write conservative variables
		double roe = p / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = 0.0;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] = roe;		//total energy equals internal one because a motion is absent

		return res;
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// 2D Kelvin-Helmholtz instability
void RunKHI2D(int argc, char *argv[]) {

	// Test parameters
	double ro_top{ 2.0 };
	double ro_bot{ 1.0 };
	double Vel{ 0.5 };
	double P{ 2.5 };
	double A{ 0.025 };				// amplitude of initial perturbation
	double lambda{ 1.0 / 6.0 };		// init. pert. wave length

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 120;
	conf.nY = 120;
	conf.LX = 1.0;
	conf.LY = 0.5;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.Gamma = 1.4;
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
	conf.SaveSolutionTime = 0.1;
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
	auto Init = [ro_top, ro_bot, Vel, P, A, lambda, &conf](Vector r) {
		// perturbation terms
		double v = A * sin(2 * PI * r.x / lambda) ;
		double w = 0;
		double u = Vel;
		double ro = ro_top;
		if (r.y < 0.5 * conf.LY) {
			u = -Vel;
			ro = ro_bot;
		};

		double roe = P / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = ro * v;
		res[3] = ro * w;
		res[4] = roe + 0.5 * ro * (u * u + v * v + w * w);
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

// 2D tests with smooth solution (Test 4.1 from Liska, Wendroff)
void RunExactEulerTest2D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 40;
	conf.nY = 40;
	conf.LX = 2.0;
	conf.LY = 2.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = true;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	conf.Gamma = 1.4;
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	//conf.methodConfiguration.RiemannProblemSolver = RPSolver::GodunovSolver;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 4.0;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 1.0;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 50;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// Test parameters
	NumericQuadrature IntVolume(5, 2);

	auto initD = [&IntVolume, &conf](Vector r) {

		// initialize density
		double ro{ 1.0 };
		ro += 0.2 * sin(PI * (r.x + r.y));
		//if((r.x * r.x + r.y * r.y) > 1) ro += 0.2 * sin(PI * (r.x + r.y)) };

		// velocity
		double u{ 1.0 };
		double v{ -0.5 };

		// pressure
		double p{ 1.0 };

		// compute ro_e and write conservative variables
		double roe = p / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = ro * v;
		res[3] = 0.0;
		res[4] = roe + 0.5 * ro * (u * u + v * v);

		return res;
	};
	kernel->SetInitialConditions(initD, IntVolume);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// Pouseuille's flow
void RunPoiseuille2D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double molarmass = 2.9e-2;
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	// Fill configuration structure
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 140;
	conf.nY = 120;
	conf.LX = 0.6;
	conf.LY = 0.1;
	conf.isPeriodicY = false;

	// Gas model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;
	conf.Viscosity = viscosity;
	conf.Sigma = Vector(sigma, 0, 0);
	conf.MolarMass = molarmass;

	// Boundary conditions
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 10.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.1;
	conf.ResidualOutputIterations = 50;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Init Conditions
	NumericQuadrature Integ(4, 2);
	auto ExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.5 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto NotExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.55 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(NotExactSol, Integ);
	kernel->SaveSolution("init.dat");

	// Velocity at center of the channel
	if (kernel->pManager->IsMaster()) std::cout << "U_max in laminar flow: " << 0.125 * conf.Sigma.x * conf.LY * conf.LY / conf.Viscosity << std::endl;

	// Create slices
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));
	kernel->SaveSliceToTecplot("slice_init.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

void RunPoiseuille3D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double molarmass = 2.9e-2;
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	//Test parameters
	//ro_init = 1.0;
	//Pave = 20.0;
	//sigma = 1.0;
	//viscosity = 0.25;

	// Fill configuration structure
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 40;
	conf.nY = 80;
	conf.nZ = 20;
	conf.LX = 0.6;
	conf.LY = 0.1;
	conf.LZ = 0.3;
	conf.isPeriodicY = false;

	// Gas model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;
	conf.Viscosity = viscosity;
	conf.Sigma = Vector(sigma, 0, 0);
	conf.MolarMass = molarmass;

	// Boundary conditions
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	//conf.SaveSolutionTime = 0.001;
	conf.SaveSliceTime = 0.001;
	conf.ResidualOutputIterations = 50;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Init Conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.5 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto NotExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.55 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(NotExactSol, Integ);
	//kernel->SaveSolution("init.dat");

	// Velocity at center of the channel
	if (kernel->pManager->IsMaster()) std::cout << "U_max in laminar flow: " << 0.125 * conf.Sigma.x * conf.LY * conf.LY / conf.Viscosity << std::endl;

	// Create slices
	//kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1 , (int)(0.25 * conf.nZ)));
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, (int)(0.5 * conf.nZ)));
	//kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, (int)(0.75 * conf.nZ)));
	//kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 1));
	kernel->SaveSliceToTecplot("slice_init.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// Shear flow
void RunShearFlow2D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double ro_init = 1.225;		// normal density
	double p_init = 1.0e5;		// normal pressure
	double u_top = 5.0;			// top plane velocity

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 40;
	conf.nY = 20;
	conf.LX = 1.0;
	conf.LY = 0.1;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = true;
	conf.qy = 1.0;

	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;

	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::MovingWall;
	conf.yRightBoundary.Gamma = 1.4;
	conf.yRightBoundary.Velocity = Vector(u_top, 0, 0);

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 20.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.05;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 1000;
	conf.DebugOutputEnabled = false;

	conf.Viscosity = viscosity;

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

	NumericQuadrature Integ(8, 2);
	auto ExactSol = [ro_init, p_init, u_top, &conf](Vector r) {
		// exact solution
		double v = 0;
		double w = 0;
		double u = r.y * u_top / conf.LY;

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto PeakWiseVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double u = r.y * u_top / conf.LY;
		if (r.y > 0.5 * conf.LY) {
			u = u_top - u;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = 0.0;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v);
		return res;
	};
	auto ContinInitVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double u = 0;
		if (r.y > 0.5 * conf.LY) {
			u = (r.y / conf.LY - 0.5) * 2 * u_top;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = 0.0;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v);
		return res;
	};
	kernel->SetInitialConditions(ContinInitVelocity, Integ);

	//save solution
	kernel->SaveSolution("init.dat");

	//Set sensors
	//auto GetXVel = [](std::valarray<double> vals) {
	//	return vals[1] / vals[0];
	//};
	//kernel->isSensorEnable = true;
	//std::unique_ptr<CellSensor> sen1 = std::make_unique<CellSensor>("Ysensor(0.5, 0.75).dat", GetXVel, kernel->nVariables);
	//sen1->SetSensor((int)(conf.nX * 0.5), (int)(conf.nY * 0.75), 0);

	//std::unique_ptr<CellSensor> sen2 = std::make_unique<CellSensor>("Ysensor(0.5, 0.5).dat", GetXVel, kernel->nVariables);
	//sen2->SetSensor((int)(conf.nX * 0.5), (int)(conf.nY * 0.5), 0);

	//std::unique_ptr<CellSensor> sen3 = std::make_unique<CellSensor>("Ysensor(0.5, 0.25).dat", GetXVel, kernel->nVariables);
	//sen3->SetSensor((int)(conf.nX * 0.5), (int)(conf.nY * 0.25), 0);
	//kernel->Sensors.push_back(std::move(sen1));
	//kernel->Sensors.push_back(std::move(sen2));
	//kernel->Sensors.push_back(std::move(sen3));
	//for(auto& r : kernel->Sensors) r->Process(kernel->values);		//initial recording

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

void RunShearFlow3D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double ro_init = 1.225;		// normal density
	double p_init = 1.0e5;		// normal pressure
	double u_top = 5.0;			// top plane velocity

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 40;
	conf.nY = 40;
	conf.nZ = 20;
	conf.LX = 1.0;
	conf.LY = 0.1;
	conf.LZ = 0.2;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = true;
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;

	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::MovingWall;
	conf.yRightBoundary.Gamma = 1.4;
	conf.yRightBoundary.Velocity = Vector(u_top, 0, 0);

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 20.0;
	conf.MaxIteration = 10000000;
	//conf.SaveSolutionTime = 0.05;
	//conf.SaveSolutionIterations = 0;
	conf.SaveSliceTime = 0.05;
	conf.ResidualOutputIterations = 100;
	conf.DebugOutputEnabled = false;

	conf.Viscosity = viscosity;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, p_init, u_top, &conf](Vector r) {
		// exact solution
		double v = 0;
		double w = 0;
		double u = r.y * u_top / conf.LY;

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto PeakWiseVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double w = 0;
		double u = r.y * u_top / conf.LY;
		if (r.y > 0.5 * conf.LY) {
			u = u_top - u;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto ContinInitVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double w = 0;
		double u = 0;
		if (r.y > 0.5 * conf.LY) {
			u = (r.y / conf.LY - 0.5) * 2 * u_top;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(ExactSol, Integ);

	// Save solution
	kernel->SaveSolution("init.dat");

	// Crete slices
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, (int)(0.5 * conf.nZ)));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// Shear flow (Z bounded for code testing)
void RunShearFlow3DZ(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double ro_init = 1.225;		// normal density
	double p_init = 1.0e5;		// normal pressure
	double u_top = 5.0;			// top plane velocity

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 20;
	conf.nY = 20;
	conf.nZ = 40;
	conf.LX = 1.0;
	conf.LY = 0.2;
	conf.LZ = 0.1;
	conf.isPeriodicZ = false;
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;
	conf.Viscosity = viscosity;

	conf.zLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.zLeftBoundary.Gamma = 1.4;
	conf.zRightBoundary.BCType = BoundaryConditionType::MovingWall;
	conf.zRightBoundary.Gamma = 1.4;
	conf.zRightBoundary.Velocity = Vector(u_top, 0, 0);

	// Method configuration
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
	conf.DummyLayerSize = 1;

	// Output info
	conf.MaxTime = 2.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.0005;
	conf.SaveSliceTime = 0.0005;
	conf.ResidualOutputIterations = 10;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, p_init, u_top, &conf](Vector r) {
		// exact solution
		double v = 0;
		double w = 0;
		double u = r.z * u_top / conf.LZ;

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto PeakWiseVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double w = 0;
		double u = r.z * 2 * u_top / conf.LZ;
		if (r.z > 0.5 * conf.LZ) {
			u = u_top;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto ContinInitVelocity = [ro_init, p_init, u_top, &conf](Vector r) {

		double v = 0;
		double u = 0;
		double w = 0;
		if (r.z > 0.5 * conf.LZ) {
			u = (r.z / conf.LZ - 0.5) * 2 * u_top;
		};

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(ExactSol, Integ);

	// Save solution
	kernel->SaveSolution("init.dat");

	// Create slices
	kernel->slices.push_back(Slice((int)(0.4 * conf.nX), (int)(0.4 * conf.nY), -1));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// 2D classic Releigh-Teylor test
void RunRTI2D(int argc, char *argv[]) {

	// Test parameters
	double ro_top = 2.0;
	double ro_bot = 1.0;
	double P_bot = 100.0;		//	bottom pressure
	double dist_amp = 0.05;
	double g = 10;				//	external uniform acceleration field

								// Domain and grid parameters
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 50;
	conf.nY = 200;
	conf.LX = 0.5;
	conf.LY = 2.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
	conf.isUniformAlongY = true;

	// Model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = false;

	// BC
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	// Method Settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

	// Computation Settings
	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSolutionIterations = 0;
	conf.ResidualOutputIterations = 50;

	// Source terms
	conf.IsUnifromAccelerationRequared = true;
	conf.UniformAcceleration = { 0, -g, 0 };

	// Kernel initializing
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// IC
	auto initD = [&conf, ro_top, ro_bot, P_bot, dist_amp, g](Vector r) {
		double ro, p;

		// define top/bottom part we consider
		double h_border = 0.5 * conf.LY;
		h_border += dist_amp * cos(2 * PI * r.x / conf.LX);

		if (r.y > h_border) {
			ro = ro_top;
			p = P_bot - ro_bot * g * h_border - ro_top * g * (r.y - h_border);
		}
		else {
			ro = ro_bot;
			p = P_bot - ro_bot * g * r.y;
		};

		double roe = p / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = 0.0;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] = roe;
		return res;
	};
	kernel->SetInitialConditions(initD);

	// Save initial solution
	kernel->SaveSolution("init.dat");

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};


////// Not widespread tests now /////////


// for Troshkin mixing simulation
void RunTurbulentMixing(int argc, char *argv[]) {
	double ro_init = 1.225;				// normal density
	double p_init = 1.0e5;				// normal pressure
	double u_top = 5.0;					// top plane velocity
	double shear_width_ratio = 0.2;		// the width of shear layer with respect to full height of domain

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 40;
	conf.nY = 40;
	conf.nZ = 40;
	conf.LX = 1.0;
	conf.LY = 0.2;
	conf.LZ = 1.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = true;
	conf.Gamma = 1.4;
	conf.MolarMass = 2.9e-2;
	
	// BC
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 20.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.001;
	//conf.SaveSolutionIterations = 0;
	//conf.SaveSliceTime = 0.05;
	conf.ResidualOutputIterations = 10;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial conditions
	double sdv = 0.1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::normal_distribution<double> normal_dist(0.0, sdv);  // N(mean, stddeviation)
	auto Init = [ro_init, p_init, u_top, shear_width_ratio, &normal_dist, &mt, &conf](Vector r) {
		double h = r.y / conf.LY;

		// Velosities
		double v = normal_dist(mt);
		double w = normal_dist(mt);
		double u = 2.0 * (h - 0.5) * u_top / shear_width_ratio;
		if (h > 0.5 * (1.0 + shear_width_ratio)) u = u_top;
		if (h < 0.5 * (1.0 - shear_width_ratio)) u = -u_top;

		double roe = p_init / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(Init);

	// Save solution
	kernel->SaveSolution("init.dat");

	// Crete slices
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, (int)(0.5 * conf.nZ)));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};
void RunKonuhovMixing(int argc, char *argv[]) {

	// Test parameters
	double ro_top{ 2.0 };
	double ro_bot{ 1.0 };
	double Vel{ 0.4 };
	double P{ 2.5 };
	double A{ 0.025 };				// amplitude of initial perturbation
	double lambda{ 1.0 / 6.0 };		// init. pert. wave length

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 120;		//120
	conf.nY = 120;		//120
	conf.nZ = 100;		//100
	conf.LX = 1.0;
	conf.LY = 0.5;
	conf.LZ = 1.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = true;
	conf.Gamma = 1.4;
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
	conf.MaxTime = 10.0;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0;
	conf.ResidualOutputIterations = 20;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);
	
	NumericQuadrature Integ(2, 3);
	auto Init = [ro_top, ro_bot, Vel, P, A, lambda, &conf](Vector r) {
		// perturbation terms
		double v = A * sin(2 * PI * r.x / lambda);
		v *= 1.0 + 0.2 * sin(2 * PI * r.z / (3 * lambda));
		double w = 0;
		double u = Vel;
		double ro = ro_top;
		if (r.y < 0.5 * conf.LY) {
			u = -Vel;
			ro = ro_bot;
		};

		double roe = P / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = ro * v;
		res[3] = ro * w;
		res[4] = roe + 0.5 * ro * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(Init, Integ);

	//save solution
	kernel->SaveSolution("init.dat");

	//Set sensor at center
	auto GetTemp = [&conf](std::valarray<double> vals) {
		auto e = vals[4] - 0.5 * (vals[1] * vals[1] + vals[2] * vals[2] + vals[3] * vals[3]) / vals[0];
		e /= vals[0];
		return conf.MolarMass * (conf.Gamma - 1) * e / (conf.UniversalGasConstant);
	};
	kernel->isSensorEnable = true;
	kernel->SaveSensorRecordIterations = 10;
	std::unique_ptr<CellSensor> sen = std::make_unique<CellSensor>("temp_centr.dat", *kernel->pManager, kernel->grid, GetTemp);
	sen->SetSensor((int)(conf.nX * 0.5), (int)(conf.nY * 0.5), (int)(conf.nZ * 0.5), 5);
	kernel->Sensors.push_back(std::move(sen));
	for (auto& r : kernel->Sensors) r->Process(kernel->values);		//initial recording

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

//// Viscous part tests	////

// Z bounded channel
void RunPoiseuille3DZ(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	//Test parameters
	ro_init = 1.0;
	Pave = 20.0;
	sigma = 1.0;
	viscosity = 0.25;

	// Fill configuration structure
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 2;
	conf.nY = 30;
	conf.nZ = 40;
	conf.LX = 0.3;
	conf.LY = 0.6;
	conf.LZ = 0.1;
	conf.isPeriodicZ = false;

	// Gas model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;
	conf.Viscosity = viscosity;
	conf.Sigma = Vector(0, sigma, 0);		// Flow in Y direction

	// Boundary conditions
	conf.zLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.zLeftBoundary.Gamma = 1.4;
	conf.zRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.zRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::ENO2PointsStencil;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.005;
	//conf.SaveSliceTime = 0.005;
	conf.ResidualOutputIterations = 50;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Init Conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.0;
		double v = 0.5 * conf.Sigma.y * r.z * (conf.LZ - r.z) / conf.Viscosity;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto NotExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.0;
		double v = 0.55 * conf.Sigma.y * r.z * (conf.LZ - r.z) / conf.Viscosity;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(NotExactSol, Integ);
	kernel->SaveSolution("init.dat");

	// Print velocity at center of the channel
	if (kernel->pManager->IsMaster()) std::cout << "U_max in laminar flow: " << 0.125 * conf.Sigma.z * conf.LX * conf.LX / conf.Viscosity << std::endl;

	// Create slices
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), (int)(0.5 * conf.nY), -1));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// X bounded channel
void RunPoiseuille3DX(int argc, char *argv[]) {
	// Air
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	// Test parameters
	ro_init = 1.0;
	Pave = 20.0;
	sigma = 1.0;
	viscosity = 0.25;

	// Fill configuration structure
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 40;
	conf.nY = 2;
	conf.nZ = 30;
	conf.LX = 0.1;
	conf.LY = 0.3;
	conf.LZ = 0.6;
	conf.isPeriodicX = false;

	// Gas model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;
	conf.Viscosity = viscosity;
	conf.Sigma = Vector(0, 0, sigma);

	// Boundary conditions
	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	//conf.SaveSolutionIterations = 1;
	//conf.SaveSolutionTime = 0.005;
	conf.SaveSliceTime = 0.005;
	conf.ResidualOutputIterations = 10;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Init Conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.0;
		double v = 0.0;
		double w = 0.5 * conf.Sigma.z * r.x * (conf.LX - r.x) / conf.Viscosity;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto NotExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.0;
		double v = 0.0;
		double w = 0.6 * conf.Sigma.z * r.x * (conf.LX - r.x) / conf.Viscosity;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(NotExactSol, Integ);
	kernel->SaveSolution("init.dat");

	if (kernel->pManager->IsMaster()) std::cout << "U_max in laminar flow: " << 0.125 * conf.Sigma.z * conf.LX * conf.LX / conf.Viscosity << std::endl;

	// Create slices
	kernel->slices.push_back(Slice(-1, (int)(0.5 * conf.nY), (int)(0.5 * conf.nZ)));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// Y bounded channel (standart channel test)
void RunPoiseuille3DY(int argc, char *argv[]) {
	// Air
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	// Test parameters
	ro_init = 1.0;
	Pave = 20.0;
	sigma = 1.0;
	viscosity = 0.25;

	// Fill configuration structure
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 30;
	conf.nY = 40;
	conf.nZ = 2;
	conf.LX = 0.6;
	conf.LY = 0.1;
	conf.LZ = 0.3;
	conf.isPeriodicY = false;

	// Gas model
	conf.Gamma = 1.4;
	conf.IsViscousFlow = false;		// to do change
	conf.Viscosity = viscosity;
	conf.Sigma = Vector(sigma, 0, 0);

	// Boundary conditions
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	//conf.SaveSolutionIterations = 1;
	//conf.SaveSolutionTime = 0.001;
	conf.SaveSliceTime = 0.001;
	conf.ResidualOutputIterations = 20;
	conf.DebugOutputEnabled = false;

	// Init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	kernel->Init(conf);

	// Init Conditions
	NumericQuadrature Integ(8, 3);
	auto ExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.5 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	auto NotExactSol = [ro_init, Pave, &conf](Vector r) {
		double u = 0.6 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;
		double w = 0.0;

		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = ro_init * w;
		res[4] = roe + 0.5 * ro_init * (u * u + v * v + w * w);
		return res;
	};
	kernel->SetInitialConditions(NotExactSol, Integ);
	kernel->SaveSolution("init.dat");

	if (kernel->pManager->IsMaster()) std::cout << "U_max in laminar flow: " << 0.125 * conf.Sigma.x * conf.LY * conf.LY / conf.Viscosity << std::endl;

	// Create slices
	kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, (int)(0.5 * conf.nZ)));
	kernel->SaveSliceToTecplot("test_slice.dat", kernel->slices[0]);

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};
*/

#endif