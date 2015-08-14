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
	conf.nX = 500;
	//conf.nY = 10;
	conf.LX = 1.0;
	//conf.LY = 1.0;
	conf.isPeriodicX = false;
	//conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
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
	conf.DummyLayerSize = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
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
	conf.nX = 2;
	//conf.nY = 10;
	conf.LX = 1.0;
	//conf.LY = 1.0;
	conf.isPeriodicX = false;
	//conf.isPeriodicY = false;
	conf.isUniformAlongX = true;
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
	conf.DummyLayerSize = 1;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// Initial Conditions
	auto initD = [&conf](Vector r) {
		double ro, u, p;
		if (r.x < 0.5 * conf.LX) {
			ro = 1.0;
			u = 1.0;
		}
		else {
			ro = 2 * r.x;
			u = 2.0;
		};
		p = 1.0;

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

void RunContactDisconTest1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 200;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.isUniformAlongX = true;
	conf.qx = 1.00;

	conf.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 2.0;
	params.PL = 1.0;
	params.uL = 1.0;
	params.roR = 1.0;
	params.PR = params.PL;
	params.uR = params.uL;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finilaze();
};

void RunShockWaves1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 100;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.isUniformAlongX = true;
	conf.qx = 1.00;

	conf.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 5;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 1.4;
	params.roL = 1.0;
	params.PL = 1.0;
	params.uL = 0.75;
	params.roR = params.roL;
	params.PR = params.PL;
	params.uR = -params.uL;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finilaze();
};

// Just one shock wave TO DO compute correct IC
void RunShockWave1D(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 100;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.isUniformAlongX = true;
	conf.qx = 1.00;

	conf.Gamma = 2.0;		// Specific Gamma

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.DummyLayerSize = 1;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.MaxTime = 0.2;
	conf.MaxIteration = 1000000;
	conf.SaveSolutionSnapshotTime = 0.1;
	conf.SaveSolutionSnapshotIterations = 5;
	conf.ResidualOutputIterations = 10;

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	// initial conditions
	double shock_vel = 0.25 * (sqrt(33.0) + 1.0);
	ShockTubeParameters params;
	params.gammaL = params.gammaR = 2.0;
	params.roL = 1.0;
	params.PL = 1.0 + shock_vel;
	params.uL = 0;
	params.roR = shock_vel / (1.0 + shock_vel);
	params.PR = 1.0;
	params.uR = -1.0;
	auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.3, params);
	kernel->SetInitialConditions(initD);

	//save solution
	kernel->SaveSolutionSega("init.dat");

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
	//ro_init = 1.0;
	//Pave = 20.0;
	//sigma = 1.0;
	//viscosity = 0.25;

	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 10;
	conf.nY = 100;
	conf.LX = 0.2;
	conf.LY = 0.1;
	conf.isPeriodicX = true;
	conf.isPeriodicY = false;
	conf.isUniformAlongY = true;
	conf.qy = 1.0;

	conf.Gamma = 1.4;
	conf.IsViscousFlow = true;

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
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

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
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
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

	//Set sensor at center
	auto GetXVel = [](std::valarray<double> vals) {
		return vals[1] / vals[0];
	};
	kernel->isSensorEnable = true;
	kernel->SaveSensorRecordIterations = 100;
	std::shared_ptr<CellSensor> sen = std::make_shared<CellSensor>("Ysensor(0.5, 0.5).dat", GetXVel, kernel->nVariables);
	sen->SetIndex(kernel->getSerialIndexGlobal((int)(conf.nX * 0.5), (int)(conf.nY * 0.5), 0));
	kernel->Sensors.push_back(sen);
	for (auto& r : kernel->Sensors) r->Process(kernel->values);		//initial recording
	
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
	conf.nY = 20;
	conf.LX = 0.5;
	conf.LY = 0.1;	
	//conf.LY = 6.0;
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
	conf.yRightBoundary.Velocity = Vector(Utop, 0, 0);

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.5;
	conf.methodConfiguration.RungeKuttaOrder = 1;
	conf.methodConfiguration.Eps = 0.05;
	conf.DummyLayerSize = 1;

	conf.MaxTime = 5.0;
	conf.MaxIteration = 10000000;
	conf.SaveSolutionSnapshotTime = 0.5;
	conf.SaveSolutionSnapshotIterations = 0;
	conf.ResidualOutputIterations = 100;

	conf.Viscosity = viscosity;
	conf.Sigma = Vector(sigma, 0, 0);

	//init kernel
	std::unique_ptr<Kernel> kernel;
	if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	kernel->Init(conf);

	double sdv = dispertion;
	std::random_device rd;
    std::mt19937 mt(rd());
	std::normal_distribution<double> normal_dist(0.0, sdv);  // N(mean, stddeviation)
	//auto initD = std::bind(SODinitialDistribution, std::placeholders::_1, 0.5, params);
	auto initD = [ro_init, Pave, Utop, &conf, &normal_dist, &mt](Vector r) {

		// solution with initial distortion
		double rnd = normal_dist(mt);
		double u = r.y * Utop * (1.0 + rnd) / conf.LY;
		rnd = normal_dist(mt);
		double v = r.y * Utop * rnd / conf.LY;
		
		// exact solution
		v = 0;
		u = r.y * Utop / conf.LY;

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



#endif