#ifndef TurboStructured_Tests_Tests
#define TurboStructured_Tests_Tests

#include <iostream>
#include <vector>
#include <random>
#include "kernel\kernel.h";
#include "Methods\ExplicitRungeKuttaFVM.h"
#include "Methods\HybridFVM.h"
#include "Methods\GeneralEosMethods\HybridGeneralEOSOnePhase.h"
#include "Methods\GeneralEosMethods\HybridBarotropicEOSOnePhase.h"
#include "Methods\GeneralEosMethods\HybridBarotropicEOSTwoPhase.h"


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
	conf.nX = 2000;
	//conf.nY = 10;
	conf.LX = 1.0;
	//conf.LY = 1.0;
	conf.isPeriodicX = false;
	//conf.isPeriodicY = false;
	conf.isUniformAlongX = false;
	conf.qx = 1.00;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
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

void RunSODTestReconstruction(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 4;
	//conf.nY = 10;
	conf.LX = 8.0;
	//conf.LY = 1.0;
	conf.isPeriodicX = false;
	//conf.isPeriodicY = false;
	conf.isUniformAlongX = false;
	conf.qx = 1.00;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
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
	auto initD = [&conf](Vector r) {
		double ro, u, p;
		if (r.x < 0.5 * conf.LX) {
			ro = r.x;
			u = 1.0;
		}
		else {
			ro = 2 * r.x;
			u = 2.0;
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
	kernel->Finilaze();
};

void RunSODTestHybrid1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 1000;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	conf.Gamma = 1.4;

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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
void RunCollisionHybrid1DBaratropic(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 500;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	//task parameters
	double E = 2.0e13;
	double p0 = 2.0e14;
	double ro0 = 11400.0;
	double uLeft = 250.0;
	double uRight = -250.0;
	BarotropicEOS eos(E, p0, ro0);

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = conf.Gamma;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = conf.Gamma;

	conf.SolutionMethod = KernelConfiguration::Method::HybridBarotropicEOSOnePhase;
	conf.methodConfiguration.CFL = 0.2;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;

	conf.MaxTime = 1.0e-5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridGeneralEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridGeneralEOSOnePhase(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridBarotropicEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridBarotropicEOSOnePhase(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	auto initD = [&kernel, &conf, &eos, uLeft, uRight](Vector r) {
		std::vector<double> res(kernel->nVariables);
		double u = 0;
		if(r.x < 0.5 * conf.LX) u = uLeft; 
		else u = uRight;
		res[0] = eos.Ro0;
		res[1] = eos.Ro0 * u;
		res[2] = 0;
		res[3] = 0;

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

void RunCollisionHybrid1DBaratropicTwoPhase(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 400;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	//task parameters
	double x0 = 0.5 * conf.LX;
	double E1 = 2.0e13;
	double E2 = E1;
	double p01 = 2.0e14;
	double p02 = p01;
	double ro01 = 7900.0;
	double ro02 = 11400;
	double uLeft = 250.0;
	double uRight = -250.0;
	BarotropicTwoPhaseEOS eos(E1, E2, p01, p02, ro01, ro02);

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = conf.Gamma;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = conf.Gamma;

	conf.SolutionMethod = KernelConfiguration::Method::HybridBarotropicEOSTwoPhase;
	conf.methodConfiguration.CFL = 0.01;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;
	conf.methodConfiguration.UseExactPressureDerivative = false;

	conf.MaxTime = 1.0e-5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel = std::unique_ptr<Kernel>(new HybridBarotropicEOSTwoPhase(&argc, &argv));
	kernel->Init(conf);

	// initial conditions
	auto initD = [&kernel, &conf, &eos, x0, uLeft, uRight, ro01, ro02](Vector r) {
		std::vector<double> res(kernel->nVariables);
		double u, ro, P, alpha;
		if(r.x < x0) {
			u = uLeft;
			ro = ro01;
			alpha = 1.0;
		} else {
			u = uRight;
			ro = ro02;
			alpha = 0;
		};
		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = ro*alpha;

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

//1D collision of identical materials by Baratropic EOS
void RunSODTestHybrid1DGeneral(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 200;
	conf.LX = 1.0;
	conf.isPeriodicX = false;

	conf.Gamma = 1.4;
	IdealGasEOS eos(conf.Gamma);

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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
	conf.SaveSolutionSnapshotTime = 0;
	conf.SaveSolutionSnapshotIterations = 1;
	conf.ResidualOutputIterations = 1;

	conf.Viscosity = viscosity;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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

	//viscosity = 2.0e-5;
	//sigma = 0.16;

	//Test parameters
	ro_init = 1.0;
	Pave = 20.0;
	sigma = 1.0;
	viscosity = 0.25;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 10;
	conf.nY = 16;
	conf.LX = 0.2;
	conf.LY = 0.1;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = false;
	conf.qy = 1.4;

	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;

	double uShear = std::sqrt(sigma * conf.LY);
	double Re = uShear * conf.LY * ro_init / viscosity;
	std::cout << "Reynolds Number = " << Re << std::endl;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.1;

	conf.MaxTime = 10.1;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 100;

	conf.Viscosity = viscosity;
	conf.IsExternalForceRequared = true;
	conf.Sigma = Vector(sigma, 0, 0);

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		double u = 0.45 * conf.Sigma.x * r.y * (conf.LY - r.y) / conf.Viscosity;
		double v = 0.0;// + normal_dist(mt);
		double w = 0.0;

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
	kernel->SaveSolutionSega("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunPoiseuille3D(int argc, char *argv[]) {
	double viscosity = 1.73e-5;	//Air
	double sigma = 0.14;		// absolute value of dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure

	KernelConfiguration conf;
	conf.nDims = 3;
	conf.nX = 10;
	conf.nY = 10;
	conf.nZ = 4;
	conf.LX = 0.2;
	conf.LY = 0.1;
	conf.LZ = 0.1;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isPeriodicZ = true;

	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.xRightBoundary.Gamma = 1.4;
	conf.yLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.yRightBoundary.Gamma = 1.4;
	conf.zLeftBoundary.BCType = BoundaryConditionType::Wall;
	conf.zLeftBoundary.Gamma = 1.4;
	conf.zRightBoundary.BCType = BoundaryConditionType::Wall;
	conf.zRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.IsExternalForceRequared = true;

	conf.MaxTime = 0.5;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.01;
	conf.SaveSolutionSnapshotIterations = 10;
	conf.ResidualOutputIterations = 1;
	conf.DebugOutputEnabled = false;

	conf.Viscosity = viscosity;
	conf.Sigma = Vector(-sigma, 0, 0);

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [ro_init, Pave, &conf](Vector r) {
		double u = -0.6*conf.Sigma.x*r.y*(conf.LY - r.y)/conf.Viscosity;
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
	kernel->SaveSolutionSega("init.dat");
	
	//run computation
	kernel->Run();		

	//finalize kernel
	kernel->Finilaze();
};

void RunShearFlow2D(int argc, char *argv[]) {
	//double viscosity = 1.73e-5;	//Air
	double viscosity = 1.0e-5;	
	double sigma = 0.0;			// dPdx
	double ro_init = 1.225;		//Air
	double Pave = 1.0e5;		//average pressure
	double Utop = 5.0;	
	double dispertion = 0.01;	//standart deviation

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 20;
	conf.nY = 100;
	conf.LX = 0.5;
	conf.LY = 0.1;	
	//conf.LY = 6.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = false;
	conf.qy = 1.01;

	conf.Gamma = 1.4;
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

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;

	conf.MaxTime = 1.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionSnapshotTime = 0.01;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	conf.Viscosity = viscosity;
	conf.Sigma = Vector(sigma, 0, 0);

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};	
	kernel->Init(conf);

	double sdv = dispertion;
	std::random_device rd;
    std::mt19937 mt(rd());
	std::normal_distribution<double> normal_dist(0.0, sdv);  // N(mean, stddeviation)
	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [ro_init, Pave, Utop, &conf, &normal_dist, &mt](Vector r) {

		double rnd = normal_dist(mt);
		double u = r.y * Utop * (1.0 + rnd) / conf.LY;
		rnd = normal_dist(mt);
		double v = r.y * Utop * rnd / conf.LY;
		//u = 0;
		double roe = Pave / (conf.Gamma - 1);
		std::vector<double> res(5);
		res[0] = ro_init;
		res[1] = ro_init * u;
		res[2] = ro_init * v;
		res[3] = 0.0;
		res[4] =  roe + 0.5 * ro_init * (u * u + v * v);
		return res; 
	};
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");

	//Set sensors
	auto GetXVel = [](std::valarray<double> vals) {
		return vals[1] / vals[0];
	};
	kernel->isSensorEnable = true;
	std::shared_ptr<CellSensor> sen1 = std::make_shared<CellSensor>("Ysensor(0.5, 0.75).dat", GetXVel, kernel->nVariables);
	sen1->SetIndex(kernel->getSerialIndexGlobal((int)(conf.nX * 0.5), (int)(conf.nY * 0.75), 0));

	std::shared_ptr<CellSensor> sen2 = std::make_shared<CellSensor>("Ysensor(0.5, 0.5).dat", GetXVel, kernel->nVariables);
	sen2->SetIndex(kernel->getSerialIndexGlobal((int)(conf.nX * 0.5), (int)(conf.nY * 0.5), 0));

	std::shared_ptr<CellSensor> sen3 = std::make_shared<CellSensor>("Ysensor(0.5, 0.25).dat", GetXVel, kernel->nVariables);
	sen3->SetIndex(kernel->getSerialIndexGlobal((int)(conf.nX * 0.5), (int)(conf.nY * 0.25), 0));
	kernel->Sensors.push_back(sen1);
	kernel->Sensors.push_back(sen2);
	kernel->Sensors.push_back(sen3);
	for(auto& r : kernel->Sensors) r->Process(kernel->values);		//initial recording

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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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

void RunDemchenkoTestBaratropicEOS(int argc, char *argv[]) {

	//collision parameters
	double E = 2.0e13;
	double p0 = 2.0e14;
	double ro0 = 11400.0;
	double uLeft = 250.0;
	double uRight = -250.0;
	BarotropicEOS eos(E, p0, ro0);

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 120;
	conf.nY = 80;
	conf.LX = 0.006;
	conf.LY = 0.002;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

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

	conf.SolutionMethod = KernelConfiguration::Method::HybridBarotropicEOSOnePhase;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;

	conf.MaxTime = 6.0e-6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridGeneralEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridGeneralEOSOnePhase(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridBarotropicEOSOnePhase) {
		kernel = std::unique_ptr<Kernel>(new HybridBarotropicEOSOnePhase(&argc, &argv));
	};
	kernel->Init(conf);
	
	auto initD = [&eos, ro0, uLeft, uRight, &conf](Vector r) {
		double roLeft = ro0;
		double roRight = ro0;
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

		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
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

void RunDemchenkoTestBaratropicEOSTwoPhase(int argc, char *argv[]) {

	//collision parameters
	double E = 2.0e13;
	double p0 = 2.0e14;
	double ro01 = 7900.0;
	double ro02 = 11400.0;
	double uLeft = 250.0;
	double uRight = -250.0;
	BarotropicTwoPhaseEOS eos(E, E, p0, p0, ro01, ro02);

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 120;
	conf.nY = 60;
	conf.LX = 0.006;
	conf.LY = 0.002;	
	conf.isPeriodicX = false;
	conf.isPeriodicY = false;

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

	conf.SolutionMethod = KernelConfiguration::Method::HybridBarotropicEOSTwoPhase;
	conf.methodConfiguration.CFL = 0.4;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;
	conf.methodConfiguration.UseExactPressureDerivative = true;

	conf.MaxTime = 6.0e-6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 1.0e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel>	kernel = std::unique_ptr<Kernel>(new HybridBarotropicEOSTwoPhase(&argc, &argv));
	kernel->Init(conf);
	
	auto initD = [&eos, ro01, ro02, uLeft, uRight, &conf](Vector r) {
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
		double alpha = 0;
		if(r.x <= bound_pos) {
			ro = ro01;
			u = uLeft;
			alpha = 1.0;
		} else {
			ro = ro02;
			u = uRight;
		};

		res[0] = ro;
		res[1] = ro * u;
		res[2] = 0;
		res[3] = 0;
		res[4] = ro * alpha;
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
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

//DoludenkoTest
void RunDoludenko2D(int argc, char *argv[]) {

	//collision parameters
	double E = 2.0e13;
	double p0 = 2.0e14;
	double ro01 = 4500.0;
	double ro02 = 11400.0;
	double pert_amplitude = 1.0e-4;
	double g = -1.0e6;
	BarotropicTwoPhaseEOS eos(E, E, p0, p0, ro01, ro02);

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 160;
	conf.nY = 80;
	conf.LX = 8.0e-3;
	conf.LY = 4.0e-3;	
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.IsUnifromAccelerationRequared = true;
	conf.UniformAcceleration = Vector(0, g, 0);

	conf.IsViscousFlow = false;

	//Symmetry for top and bottom faces and Natural for left and right ones
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridBarotropicEOSTwoPhase;
	conf.methodConfiguration.CFL = 0.2;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.eos = &eos;
	conf.methodConfiguration.UseExactPressureDerivative = true;

	conf.MaxTime = 1.0e-6;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1e-6;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel>	kernel = std::unique_ptr<Kernel>(new HybridBarotropicEOSTwoPhase(&argc, &argv));
	kernel->Init(conf);
	
	auto initD = [&eos, ro01, ro02, pert_amplitude, &conf](Vector r) {
		double yCentr = 0.5 * conf.LY;	//a half of a tube height
		double ro = 0;
		double alpha = 0;
		if (r.y < yCentr) {
			ro = ro01;
			alpha = 1;
		} else ro = ro02;

		double v = pert_amplitude * cos(1200 * (r.x - 0.5 * conf.LX)) * exp(-abs(r.y - yCentr));
		
		std::vector<double> res(5);
		res[0] = ro;
		res[1] = 0;
		res[2] = ro * v;
		res[3] = 0;
		res[4] = ro * alpha;
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

//2D and 3D releigh_taylor classic test
void RunReleighTaylor2D(int argc, char *argv[]) {

	//Collision parameters
	double ro_top = 2.0;		//Collision speed of lest material
	double ro_bottom = 1.0;	// and right one
	double perturbation_amplitude = 0.01;
	double pressure0 = 2.5;	//initial common pressure
	double g = -0.25;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 60;
	conf.nY = 180;
	conf.LX = 0.5;
	conf.LY = 1.5;	
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;

	conf.Gamma = 1.4;
	conf.IsViscousFlow = false;

	//Symmetry for top and bottom faces and Natural for left and right ones
	conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yLeftBoundary.Gamma = 1.4;
	conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
	conf.yRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::HybridFVM;
	conf.methodConfiguration.CFL = 0.2;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.UniformAcceleration = Vector(0, g, 0);

	conf.MaxTime = 10.0;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.5;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;
	conf.IsUnifromAccelerationRequared = true;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.SolutionMethod == KernelConfiguration::Method::HybridFVM) {
		kernel = std::unique_ptr<Kernel>(new HybridFVM(&argc, &argv));
	};
	kernel->Init(conf);
	
	auto initD = [ro_top, ro_bottom, pressure0, perturbation_amplitude, &conf](Vector r) {
		std::vector<double> res(5);

		double ro = 0;
		if(r.y <= 0.5 * conf.LY) ro = ro_bottom;
		else ro = ro_top;

		double v = (1.0 - cos(2.0 * PI * r.x / conf.LX)) * (1.0 - cos(2.0 * PI * r.y / conf.LY));
		v *= 0.25 * perturbation_amplitude;

		double P = pressure0 - ro_bottom * conf.UniformAcceleration.y * conf.LY + ro * conf.UniformAcceleration.y * r.y; 
		double roe = P / (conf.Gamma - 1.0);
		res[0] = ro;
		res[1] = 0;
		res[2] = ro * v;
		res[3] = 0;
		res[4] = roe + 0.5 * ro * v * v;
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