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
	conf.nX = 100;
	conf.isPeriodicX = false;
	conf.DummyLayerSize = 1;

	// Describe grid compression here
	BlockNode nleft, ncenter;
	nleft.N_cells = conf.nX / 2;
	nleft.q_com = 1.0 / 1.05;
	ncenter.pos = 0.5 * conf.LX;
	ncenter.N_cells = conf.nX - nleft.N_cells;
	ncenter.q_com = 1.0 / nleft.q_com;
	conf.CompressionX[0] = nleft;
	conf.CompressionX.push_back(ncenter);

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.xLeftBoundary.SetMarker(1);
	conf.xRightBoundary.SetMarker(1);

	// Model settings
	conf.Gamma = 1.4;
	conf.IsViscousFlow = false;

	// Method settings
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;

	// Task settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionTime = 0.25;
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 100;

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

	//save init solution and run the test
	kernel->SaveSolutionToTecplot("init.dat");

	//save solution
	kernel->SaveSolution("sol.sol");
	kernel->LoadSolution("sol.sol");
	kernel->SaveSolutionToTecplot("init2.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 50;

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
	kernel->SaveSolutionToTecplot("init.dat");

	// Run computation
	kernel->Run();

	// Finalize kernel
	kernel->Finalize();
};

// YZ direction test
void RunSODXTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 100;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.DummyLayerSize = 1;

	// Describe grid compression here
	BlockNode nleft, ncenter;
	nleft.N_cells = conf.nX / 2;
	nleft.q_com = 1.0 / 1.05;
	ncenter.pos = 0.5 * conf.LX;
	ncenter.N_cells = conf.nX - nleft.N_cells;
	ncenter.q_com = 1.0 / nleft.q_com;
	conf.CompressionX[0] = nleft;
	conf.CompressionX.push_back(ncenter);

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.xLeftBoundary.SetMarker(1);
	conf.xRightBoundary.SetMarker(1);

	// Method parameters
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.3;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	// Model
	conf.Gamma = 1.4;

	// Task and output settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	//conf.SaveSolutionTime = 0.05;
	conf.SaveSliceTime = 0.25;
	conf.ResidualOutputIters = 100;

	// Init kernel
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

	// initial conditions
	struct ShockTubeParameters {
		double roL = 1.0;
		double PL = 1.0;
		Vector uL = { 0, 0, 0 };
		double roR = 0.125;
		double PR = 0.1;
		Vector uR = { 0, 0, 0 };
	} params;

	double x0 = 0.5 * conf.LX;
	auto init = [&params, &conf, x0](Vector r) {
		std::vector<double> res(5);
		double rho, p;
		double gamma = conf.Gamma;
		Vector V;
		if (r.x < x0) {
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
	kernel->SetInitialConditions(init);
	
	// Init the slice and save initial solution
	kernel->slices.push_back(Slice(-1, 0, 0));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};
void RunSODYTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 1;
	conf.nY = 100;
	conf.LX = 1.0;
	conf.LY = 1.0;
	conf.isPeriodicY = false;

	// Describe grid compression here
	BlockNode nleft, ncenter;
	nleft.N_cells = conf.nY / 2;
	nleft.q_com = 1.0 / 1.05;
	ncenter.pos = 0.5 * conf.LY;
	ncenter.N_cells = conf.nY - nleft.N_cells;
	ncenter.q_com = 1.0 / nleft.q_com;
	conf.CompressionY[0] = nleft;
	conf.CompressionY.push_back(ncenter);

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.yLeftBoundary.SetMarker(1);
	conf.yRightBoundary.SetMarker(1);

	// Method parameters
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.3;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	// Model
	conf.Gamma = 1.4;

	// Task and output settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	//conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.25;
	conf.ResidualOutputIters = 100;

	// Init kernel
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

	// initial conditions
	struct ShockTubeParameters {
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
		double gamma = conf.Gamma;
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
	kernel->SetInitialConditions(init);


	// Init the slice and save initial solution
	kernel->slices.push_back(Slice(1, -1, 0));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};
void RunSODZTest(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 1;
	conf.nY = 1;
	conf.nZ = 100;
	conf.LX = 1.0;
	conf.LY = 1.0;
	conf.LZ = 1.0;
	conf.isPeriodicZ = false;

	// Describe grid compression here
	BlockNode nleft, ncenter;
	nleft.N_cells = conf.nZ / 2;
	nleft.q_com = 1.0 / 1.05;
	ncenter.pos = 0.5 * conf.LZ;
	ncenter.N_cells = conf.nZ - nleft.N_cells;
	ncenter.q_com = 1.0 / nleft.q_com;
	conf.CompressionZ[0] = nleft;
	conf.CompressionZ.push_back(ncenter);

	// BC
	conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);
	conf.zLeftBoundary.SetMarker(1);
	conf.zRightBoundary.SetMarker(1);

	// Method parameters
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.3;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	// Model
	conf.Gamma = 1.4;

	// Task and output settings
	conf.MaxTime = 0.25;
	conf.MaxIteration = 1000000;
	//conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.25;
	conf.ResidualOutputIters = 100;

	// Init kernel
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

	// initial conditions
	struct ShockTubeParameters {
		double roL = 1.0;
		double PL = 1.0;
		Vector uL = { 0, 0, 0 };
		double roR = 0.125;
		double PR = 0.1;
		Vector uR = { 0, 0, 0 };
	} params;

	double z0 = 0.5 * conf.LZ;
	auto init = [&params, &conf, z0](Vector r) {
		std::vector<double> res(5);
		double rho, p;
		double gamma = conf.Gamma;
		Vector V;
		if (r.z < z0) {
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
	kernel->SetInitialConditions(init);


	// Init the slice and save initial solution
	kernel->slices.push_back(Slice(1, 1, -1));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finalize();
};

// Poiseuille flow (1D)
void RunPoiseuille1D(int argc, char *argv[]) {

	// Define test parameters and set all fields
	struct parameters {
		double gamma;
		double Lx;			// domain size
		double Ly;			// domain size
		double ro;			// init density
		double Pin;			// inlet pressure
		double sigma;		// (Pin - Pout) / Lx 
		double viscosity;	// dynamic viscosity
	} par;
	par.gamma = 1.4;
	par.Lx = 0.5;
	//par.Ly = 0.1;
	par.Pin = 101579;
	par.sigma = 2.0;
	par.ro = 1.0;
	par.viscosity = 0.01;

	// Init config structure
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 400;
	conf.LX = par.Lx;
	//conf.LY = par.Ly;
	conf.isPeriodicX= false;
	conf.Gamma = par.gamma;
	conf.IsViscousFlow = true;
	conf.Viscosity = par.viscosity;
	conf.Sigma = Vector(0, par.sigma, 0);
	conf.IsExternalForceRequared = true;

	// Discribe grid compression
	BlockNode x_bot, x_cen;
	x_bot.N_cells = 0.5 * conf.nX;
	x_bot.q_com = 1.01;
	conf.CompressionX[0] = x_bot;

	x_cen.pos = 0.5 * par.Lx;
	x_cen.q_com = 1.0 / x_bot.q_com;
	x_cen.N_cells = conf.nX - x_bot.N_cells;
	conf.CompressionX.push_back(x_cen);

	// Discribe boundary conditions
	BoundaryConditionConfiguration Wall(BoundaryConditionType::Wall);
	conf.MyConditions[1] = Wall;

	// Describe conditions on the domain sides
	conf.xLeftBoundary.SetMarker(1);
	conf.xRightBoundary.SetMarker(1);

	// Method settings
	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.45;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
	conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
	conf.DummyLayerSize = 1;

	// compute velocity at center of the channel for laminar regime
	auto Uc = 0.125 * par.sigma * par.Lx * par.Lx / par.viscosity;
	auto Re = par.ro * Uc * par.Lx / par.viscosity;

	// Computational settings
	conf.MaxTime = 10;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionTime = 0.1;
	conf.SaveSliceTime = 0.1;
	conf.SaveBinarySolIters = 100000;
	conf.ResidualOutputIters = 10000;

	// init kernel
	auto kernel = CreateKernel(conf, argc, argv);
	kernel->Init(conf);

	// init distributions
	NumericQuadrature Integ(5, 2);
	auto InitTrivial = [&par](Vector r) {
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
	auto InitExact = [&par, Uc](Vector r) {
		auto x = 1.0 - abs(2 * r.x / par.Lx - 1.0);
		auto vel = Uc * x * (2.0 - x);

		// Compute primitive values
		double rho = par.ro;
		double u = 0;
		double v = vel;
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
	auto InitHigher = [&par, Uc](Vector r) {
		auto x = 1.0 - abs(2 * r.x / par.Lx - 1.0);
		auto vel = 1.1 * Uc * x * (2.0 - x);

		// Compute primitive values
		double rho = par.ro;
		double u = 0;
		double v = vel;
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
	auto InitLower = [&par, Uc](Vector r) {
		auto x = 1.0 - abs(2 * r.x / par.Lx - 1.0);
		auto vel = 0.9 * Uc * x * (2.0 - x);

		// Compute primitive values
		double rho = par.ro;
		double u = 0;
		double v = vel;
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
	kernel->SetInitialConditions(InitExact, Integ);

	// Create slices
	//kernel->slices.push_back(-1, 0, 0));

	//save init solution and run the test
	kernel->SaveSolution("solution.sol");
	//kernel->LoadSolution("sol_to_load.sol");
	kernel->SaveSolutionToTecplot("init.dat");

	// Run test
	if (kernel->pManager->IsMaster()) std::cout << "Poiseuille test runs." << std::endl <<
		"U_center: " << Uc << ", Re: " << Re << std::endl;
	kernel->Run();


	//end of experiments
	std::cout << "Poiseuille test is completed";
	std::cout << std::endl;
	kernel->Finalize();
};

//// Test 2D

// Temporal test for code2code comparison on coarse grid
namespace BlasiusFlowTestDebug {

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
		auto Udr = 10.0;

		par.gamma = 1.4;
		par.Lx = 0.4;
		par.Ly = 0.5;
		par.Xplate = 0.2;
		par.Pin = 101579;
		par.Pout = par.Pin;
		par.ro = par.Pin * par.gamma / (1006.43 * 300.214 * (par.gamma - 1));
		par.viscosity = 1.7894e-03;
		par.M = Udr / sqrt(par.gamma * par.Pin / par.ro);		// Udriven = 10
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
		conf.nX = 4;
		conf.nY = 2;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.Gamma = par.gamma;
		conf.IsViscousFlow = false;		// TO DO change
		//conf.Viscosity = par.viscosity;

		// Describe grid compression
		BlockNode beforePlate, startPlate, bottomNode;

		// X direction first
		beforePlate.N_cells = 2;
		beforePlate.q_com = 0.5;
		conf.CompressionX[0] = beforePlate;

		startPlate.pos = par.Xplate;
		startPlate.q_com = 1.5;
		startPlate.N_cells = conf.nX - beforePlate.N_cells;
		conf.CompressionX.push_back(startPlate);

		// Y
		bottomNode.q_com = 2.0;
		bottomNode.N_cells = conf.nY;
		conf.CompressionY[0] = bottomNode;

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
		conf.MaxIteration = 25;
		conf.SaveSolutionTime = 0.1;
		conf.SaveSolutionIters = 1;
		conf.ResidualOutputIters = 1;

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
		auto Init = [Udr](Vector r) {
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
		kernel->SetInitialConditions(Init, Integ);

		//save init solution and run the test
		kernel->SaveSolutionToTecplot("init.dat");

		// Run test
		if (kernel->pManager->IsMaster()) std::cout << "Flat plate test runs." << std::endl <<
			"Inlet velocity: " << Udr << std::endl;
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

// Temporal PF
namespace PoiseuilleFlowTemp {

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
		par.Lx = 10.0;
		par.Ly = 5.0;
		par.Pin = 101579;
		par.sigma = 20.0;
		par.ro = 1.0;
		par.viscosity = 1.0;
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
		conf.nX = 80;
		conf.nY = 400;
		//conf.nY = 200;
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
		y_bot.q_com = 1.01;
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
		conf.MaxTime = 2.0;
		conf.MaxIteration = 10000000;
		conf.SaveSolutionTime = 0.01;
		conf.SaveSliceTime = 0.01;
		conf.SaveBinarySolIters = 100000;
		conf.ResidualOutputIters = 100;
		conf.ChangeHistoryFileIters = 1000000;
		conf.SaveHistoryFileIters = 50;
		//conf.ContinueComputation = true;

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
			auto y = 1.0 - abs(2 * r.y / par.Ly - 1.0);
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
			auto vel = 1.1 * Uc * y * (2.0 - y);

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
		if (kernel->ContinueComputation == true) kernel->LoadSolution("sol_to_load.sol");
		else kernel->SetInitialConditions(InitExact, Integ);

		// Create slices
		kernel->slices.push_back(Slice((int)(0.5 * conf.nX), -1, 0));

		//save init solution and run the test
		kernel->SaveSolution("solution.sol");
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


/*
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
	conf.ResidualOutputIters = 10;

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
	//kernel->SaveSolutionToTecplot("init.dat");
	kernel->slices.push_back(Slice(1, -1, 0));
	kernel->SaveSliceToTecplot("init_slice.dat", kernel->slices[0]);

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 10;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 20;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 10;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.ResidualOutputIters = 50;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 50;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.ResidualOutputIters = 50;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.ResidualOutputIters = 50;
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
	//kernel->SaveSolutionToTecplot("init.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 1000;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	//conf.SaveSolutionIters = 0;
	conf.SaveSliceTime = 0.05;
	conf.ResidualOutputIters = 100;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.ResidualOutputIters = 10;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.SaveSolutionIters = 0;
	conf.ResidualOutputIters = 50;

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
	kernel->SaveSolutionToTecplot("init.dat");

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
	//conf.SaveSolutionIters = 0;
	//conf.SaveSliceTime = 0.05;
	conf.ResidualOutputIters = 10;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	conf.ResidualOutputIters = 20;

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
	kernel->SaveSolutionToTecplot("init.dat");

	//Set sensor at center
	auto GetTemp = [&conf](std::valarray<double> vals) {
		auto e = vals[4] - 0.5 * (vals[1] * vals[1] + vals[2] * vals[2] + vals[3] * vals[3]) / vals[0];
		e /= vals[0];
		return conf.MolarMass * (conf.Gamma - 1) * e / (conf.UniversalGasConstant);
	};
	kernel->isSensorEnable = true;
	kernel->SaveSensorRecordIters = 10;
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
	conf.ResidualOutputIters = 50;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	//conf.SaveSolutionIters = 1;
	//conf.SaveSolutionTime = 0.005;
	conf.SaveSliceTime = 0.005;
	conf.ResidualOutputIters = 10;
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
	kernel->SaveSolutionToTecplot("init.dat");

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
	//conf.SaveSolutionIters = 1;
	//conf.SaveSolutionTime = 0.001;
	conf.SaveSliceTime = 0.001;
	conf.ResidualOutputIters = 20;
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
	kernel->SaveSolutionToTecplot("init.dat");

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