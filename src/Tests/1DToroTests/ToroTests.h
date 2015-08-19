#ifndef TurboStructured_Tests_1DToroTests_IToroTest
#define TurboStructured_Tests_1DToroTestsIToroTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h";
#include "Methods\ExplicitRungeKuttaFVM.h"

// some usefull structures
struct ShockTubeParameters {
	double gammaL;
	double roL;
	double PL;
	Vector uL;
	double roR;
	double PR;
	Vector uR;
	double gammaR;
};


// Error norm functions
double ComputeL2Error(std::valarray<double>& comp_val, std::valarray<double>& exac_val, int nV) {
	return 0;
}

double ComputeLinfError(std::valarray<double>& comp_val, std::valarray<double>& exac_val, int nV) {
	return 0;
}

std::pair<double, double> RunToroTest1(int argc, char *argv[]) {
	KernelConfiguration conf;
	conf.nDims = 1;
	conf.nX = 200;
	conf.LX = 1.0;
	conf.isPeriodicX = false;
	conf.isUniformAlongX = true;
	conf.qx = 1.00;
	conf.Gamma = 1.4;

	conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
	conf.xLeftBoundary.Gamma = 1.4;
	conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
	conf.xRightBoundary.Gamma = 1.4;

	conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
	conf.methodConfiguration.CFL = 0.4;
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
		kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
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

	//save solution
	kernel->SaveSolution("init.dat");

	//run computation
	kernel->Run();

	//finalize kernel
	kernel->Finilaze();

	return{};
};


#endif
