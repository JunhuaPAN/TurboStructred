#ifndef TurboStructured_Tests_Tests
#define TurboStructured_Tests_Tests

#include <iostream>
#include <vector>
#include <random>
#include "kernel\kernel.h";
#include "Methods\ExplicitRungeKuttaFVM.h"
#include "Methods\HybridFVM.h"
#include "Methods\GeneralEosMethods\HybridGeneralEOSOnePhase.h"

#define PI 3.14159265359

//collect here all tests now

//some usefull structures
struct ShockTubeParameters {
	double gammaL;
	double roL;
	double PL; 
	double uL; 
	double roR;
	double PR;
	double uR;
	double gammaR;
};

std::vector<double> SODinitialDistribution(Vector r, double xI, ShockTubeParameters params) {
	std::vector<double> U(5);	
	double ro;
	double u;
	double gamma;
	double P;

	//Left
	if (r.x <= xI) {
		ro = params.roL;
		u = params.uL;
		P = params.PL;
		gamma = params.gammaL;
	} else {
		ro = params.roR;
		u = params.uR;
		P = params.PR;
		gamma = params.gammaR;
	};

	U[0] = ro;
	U[1] = u * ro;
	U[2] = 0;
	U[3] = 0;
	U[4] = P / (gamma - 1.0) + ro * u * u / 2.0;

	return U;
};

std::vector<double> SODinitialDistributionY(Vector r, double yI, ShockTubeParameters params) {
	std::vector<double> U(5);	
	double ro;
	double u;
	double gamma;
	double P;

	//Left
	if (r.y <= yI) {
		ro = params.roL;
		u = params.uL;
		P = params.PL;
		gamma = params.gammaL;
	} else {
		ro = params.roR;
		u = params.uR;
		P = params.PR;
		gamma = params.gammaR;
	};

	U[0] = ro;
	U[1] = 0;
	U[2] = u * ro;
	U[3] = 0;
	U[4] = P / (gamma - 1.0) + ro * u * u / 2.0;

	return U;
};

//Tests
void RunSODTestRoe1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 10;
	conf.nY = 10;
	conf.LX = 1.0;
	conf.LY = 1.0;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunSODTestHybrid1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 1000;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.03;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 100;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	//params.uL = 0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};
void RunSODTestHybrid1DY(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 1;
	conf.nY = 200;
	conf.LX = 0.5;
	conf.LY = 1.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::SymmetryX;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::SymmetryX;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 100;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	//params.uL = 0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistributionY, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};
void RunSODTestHybrid1DX(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 200;
	conf.nY = 1;
	conf.LX = 1.0;
	conf.LY = 0.5;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 100;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	//params.uL = 0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

//SOD Test for general EOS
void RunSODTestHybrid1DGeneral(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 200;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	conf.Gamma = 1.4;
	IdealGasEOS eos(conf.Gamma);
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = conf.Gamma;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = conf.Gamma;

	conf.SolutionMethod = KernelConfiguration::Method::HybridGeneralEOSOnePhase;
	conf.methodConfiguration.CFL = 0.1;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.UseExactPressureDerivative = false;
	conf.methodConfiguration.eos = &eos;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridGeneralEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridGeneralEOSOnePhase(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	//params.uL = 0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunSODTestHybrid2D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 500;
	conf.nY = 10;
	conf.LX = 1.0;
	conf.LY = 2.0;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.05;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	//params.uL = 0.75;
	params.uL = 0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunSODTestRoe2DX(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 100;
	conf.nY = 100;
	conf.LX = 1.0;
	conf.LY = 1.0;	
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

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
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.0;
	params.roR = 0.125;
	params.PR = 0.1;
	params.uR = 0.0;
	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [](Vector r) { 
		std::vector<double> res(5);
		res[0] = 1.0;
		res[1] = 0.0;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] = 1.0 / (0.4);
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

//test for comparison of fluxes in different codes
void RunFluxesTest2D(int argc, char *argv[]) {
	double viscosity = 2.0;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 4;
	conf.nY = 4;
	conf.LX = 2.0;
	conf.LY = 1.0;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

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
	conf.SaveSolutionSnapshotTime = 0;
	conf.SaveSolutionSnapshotIterations = 1;
	conf.ResidualOutputIterations = 1;

	conf.Viscosity = viscosity;
	conf.Sigma = 0;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [](Vector r) {
		double u = r.x + r.y;
		double v = r.x - r.y;
		double w = 0;
		double ro = 1.5*(1.0 + r.y*0.01);
		double roe = ro*100;
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
	kernel->Finilaze();
};

void RunPoiseuille2DFVM(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// -dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure		

	//Test
	ro_init = 1.0;
	Pave = 20.0;
	sigma = 1.0;
	viscosity = 0.000025;	

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 100;
	conf.nY = 100;
	conf.LX = 0.2;
	conf.LY = 0.1;	
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;
	conf.IsViscousFlow = true;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	//Example of moving wall BC
	//conf.yLeftBoundary.BCType = BoundaryConditionType::MovingWall;
	//conf.yLeftBoundary.Gamma = 1.4;
	//conf.yLeftBoundary.Velocity = Vector(5, 0, 0);

	double uShear = std::sqrt(sigma * conf.LY);
	double Re = uShear * conf.LY * ro_init / viscosity;
	std::cout<<"Reynolds Number = "<<Re<<std::endl;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 0.01;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1;
	conf.SaveSolutionSnapshotIterations = 100;
	conf.ResidualOutputIterations = 10;

	conf.Viscosity = viscosity;
	conf.Sigma = sigma;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);
	
	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	double sdv = 0.001;
	std::random_device rd;
    std::mt19937 mt(rd());
	std::normal_distribution<double> normal_dist(0.0, sdv);  // N(mean, stddeviation)
	
	auto initD = [ro_init, Pave, &conf, &normal_dist, &mt](Vector r) {
		double u = 0.5*conf.Sigma*r.y*(conf.LY - r.y)/conf.Viscosity;
		double v = 0.0 + normal_dist(mt);
		double w = 0.0;
		//u = 0.01;
		double roe = Pave/(conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = res[0] * u;
		res[2] = res[0] * v;
		res[3] = res[0] * w;
		res[4] = roe + 0.5 * res[0] * u * u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunPoiseuille3D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 10;
	conf.nY = 40;
	conf.nZ = 6;
	conf.LX = 0.2;
	conf.LY = 0.1;
	conf.LZ = 0.1;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = true;

	conf.Gamma = 1.4;
	conf.nVariables = 5;

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

	conf.MaxTime = 0.5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.01;
	conf.SaveSolutionSnapshotIterations = 1000;
	conf.ResidualOutputIterations = 1;

	conf.Viscosity = viscosity;
	conf.Sigma = sigma;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [ro_init, Pave, &conf](Vector r) {
		double u = 0.5*conf.Sigma*r.y*(conf.LY - r.y)/conf.Viscosity;
		u = 0;
		double roe = Pave/(conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = res[0]*u;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] = roe + 0.5*res[0]*u*u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunShearFlow2D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.0;		// dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure
	double Utop = 1.0;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 10;
	conf.nY = 40;
	conf.LX = 0.2;
	conf.LY = 0.1;	
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;
	conf.IsViscousFlow = true;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::MovingWall;
	conf.yRightBoundary.Gamma = 1.4;
	conf.yRightBoundary.Velocity = Vector(Utop, 0, 0);

	//Example of moving wall BC
	//conf.yLeftBoundary.BCType = BoundaryConditionType::MovingWall;
	//conf.yLeftBoundary.Gamma = 1.4;
	//conf.yLeftBoundary.Velocity = Vector(5, 0, 0);

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.01;
	conf.SaveSolutionSnapshotIterations = 1000;
	conf.ResidualOutputIterations = 50;

	conf.Viscosity = viscosity;
	conf.Sigma = sigma;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [ro_init, Pave, Utop, &conf](Vector r) {
		double u = r.y * Utop / conf.LY;
		//u = 0;
		double roe = Pave/(conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init*u;
		res[2] = 0.0;
		res[3] = 0.0;
		res[4] =  roe + 0.5*ro_init*u*u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

//Demchenko test

//2D & 3D collision. First medium (steel) collides with second medium (leed)
void RunDemchenkoTest2D(int argc, char *argv[]) {

	//Collision parameters
	double uLeft = 250;		//Collision speed of lest material
	double uRight = -250;	// and right one
	double roLeft = 7900;	// SI for Steel
	double roRight = 11340; // SI	for Pb
	double pressure = 4e8;	//initial common pressure	

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 120;
	conf.nY = 60;
	conf.LX = 0.006;
	conf.LY = 0.002;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;
	conf.IsViscousFlow = false;

	//Symmetry for top and bottom faces and Natural for left and right ones
	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 6.0e-6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);
	
	auto initD = [roLeft, roRight, uLeft, uRight, pressure, &conf](Vector r) {
		double x0 = 0.5 * conf.LX;		//position of distorbancesfree surface
		double bound_pos = 0;			//position of free surface
		double yCentr = 0.5 * conf.LY;	//a half of a tube height
		if (r.y < yCentr) {
			bound_pos = 1 - cos((r.y - yCentr) * PI / yCentr);
			bound_pos *= 0.5 * yCentr;
		};
		bound_pos += x0;

		std::vector<double> res(5);
		double ro = 0;
		double u = 0;
		if(r.x <= bound_pos) {
			ro = roLeft;
			u = uLeft;
		} else {
			ro = roRight;
			u = uRight;
		};
		double roe = pressure / (conf.Gamma - 1.0);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = roe + 0.5 * ro * u * u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunDemchenkoTest2DGeneralEOS(int argc, char *argv[]) {

	//Collision parameters
	double uLeft = 250;		//Collision speed of lest material
	double uRight = -250;	// and right one
	double roLeft = 7900;	// SI for Steel
	double roRight = 11340; // SI	for Pb
	double pressure = 4e8;	//initial common pressure	

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 200;
	conf.nY = 100;
	conf.LX = 0.006;
	conf.LY = 0.002;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	IdealGasEOS eos(conf.Gamma);
	conf.nVariables = 5;
	conf.IsViscousFlow = false;

	//Symmetry for top and bottom faces and Natural for left and right ones
	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridGeneralEOSOnePhase;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;
	conf.methodConfiguration.UseExactPressureDerivative = true;

	conf.MaxTime = 6.0e-6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridGeneralEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridGeneralEOSOnePhase(&argc, &argv));
	};
	kernel->Init(conf);
	
	auto initD = [&eos, roLeft, roRight, uLeft, uRight, pressure, &conf](Vector r) {
		double x0 = 0.5 * conf.LX;		//position of distorbancesfree surface
		double bound_pos = 0;			//position of free surface
		double yCentr = 0.5 * conf.LY;	//a half of a tube height
		if (r.y < yCentr) {
			bound_pos = 1 - cos((r.y - yCentr) * PI / yCentr);
			bound_pos *= 0.5 * yCentr;
		};
		bound_pos += x0;

		std::vector<double> res(5);
		double ro = 0;
		double u = 0;
		if(r.x <= bound_pos) {
			ro = roLeft;
			u = uLeft;
		} else {
			ro = roRight;
			u = uRight;
		};
		double roe = ro * eos.GetInternalEnergy(ro, pressure);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = roe + 0.5 * ro * u * u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunDemchenkoTest3D(int argc, char *argv[]) {

	//Collision parameters
	double uLeft = 250;		//Collision speed of lest material
	double uRight = -250;	// and right one
	double roLeft = 7900;	// SI for Steel
	double roRight = 11340; // SI	for Pb
	double pressure = 4e8;	//initial common pressure	

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 40;
	conf.nY = 80;
	conf.nZ = 40;
	conf.LX = 0.006;
	conf.LY = 0.002;
	conf.LZ = 0.002;
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = false;

	conf.Gamma = 1.4;
	conf.nVariables = 5;
	conf.IsViscousFlow = false;

	//Symmetry for top and bottom faces and Natural for left and right ones
	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;
	conf.zLeftBoundary.BCType = BoundaryConditionType::SymmetryZ;
	conf.zLeftBoundary.Gamma = 1.4;
	conf.zRightBoundary.BCType = BoundaryConditionType::SymmetryZ;
	conf.zRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 1.0e-5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 10000;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);
	
	auto initD = [roLeft, roRight, uLeft, uRight, pressure, &conf](Vector r) {
		double x0 = 0.5 * conf.LX;		//position of distorbancesfree surface
		double bound_pos = 0;			//position of free surface
		double radius = 0.5 * conf.LY;	//a half of a tube height
		double DistanceFromX = sqrt( r.y * r.y + r.z * r.z );
		if ( DistanceFromX < radius) {
			bound_pos = 1 - cos((DistanceFromX - radius) * PI / radius);
			bound_pos *= 0.5 * radius;
		};
		bound_pos += x0;

		std::vector<double> res(5);
		double ro = 0;
		double u = 0;
		if(r.x <= bound_pos) {
			ro = roLeft;
			u = uLeft;
		} else {
			ro = roRight;
			u = uRight;
		};
		double roe = pressure / (conf.Gamma - 1.0);
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = roe + 0.5 * ro * u * u;
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};



#endif