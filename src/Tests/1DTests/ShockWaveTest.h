#ifndef TurboStructured_Tests_1DToroTests_IToroTest
#define TurboStructured_Tests_1DToroTestsIToroTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h";
#include "Methods\ExplicitRungeKuttaFVM.h"

namespace ShockWaveTest
{
	enum ShockDirection { left, right };

	//collect here all tests now
	struct ShockTubeParameters {
		double gamma;
		double roL;
		double roR;
		double P;
		double uL;
		double uR;
		double x0;
		double ref_velocity;
	};

	KernelConfiguration DefaultSettings() {
		KernelConfiguration conf;

		// Values by default
		conf.nDims = 1;
		conf.nX = 200;
		conf.LX = 1.0;
		conf.isPeriodicX = false;
		conf.isUniformAlongX = true;
		conf.qx = 1.00;
		conf.Gamma = 1.4;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.ReconstructionType = KernelConfiguration::Reconstruction::PiecewiseConstant;
		//conf.ReconstructionType = KernelConfiguration::Reconstruction::ENO2PointsStencil;
		conf.DummyLayerSize = 1;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;

		conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
		conf.xLeftBoundary.Gamma = 1.4;
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = 1.4;

		conf.MaxTime = 2.0;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 0.1;
		conf.SaveSolutionSnapshotIterations = 0;
		conf.ResidualOutputIterations = 10;
		
		return conf;
	};

	ShockTubeParameters DefaultStates() {
		ShockTubeParameters params;
		params.roL = 1.0;
		params.roR = 2.0;
		params.P = 1.0;
		params.uL = 1.75;
		params.uR = -1.75;
		params.x0 = 0.2;		// initial discontinuity position
		params.gamma = 1.4;
		params.ref_velocity = 0;

		return params;
	};

	// return velocity of the Jump
	double ApproximateSolution(ShockDirection Dir, double &ro_st, double &u_st, double &p_st) {
		ShockTubeParameters par = DefaultStates();
		double ro0, u0, p0;
		if (Dir == ShockDirection::right) {
			ro_st = 6.62427997;
			u_st = -0.300254;
			p_st = 7.0214958;

			ro0 = par.roR;
			u0 = par.uR;
			p0 = par.P;
		};

		// Apply RH condition to find shock velocity
		double g1 = par.gamma - 1.0;

		double dU = ro0 - ro_st;
		double dF = ro0 * u0 - ro_st * u_st;
		double D1 = dF / dU;

		dU = ro0 * u0 - ro_st * u_st;
		dF = p0 - p_st + ro0 * u0 * u0 - ro_st * u_st * u_st;
		double D2 = dF / dU;

		double roe0 = p0 / g1;
		double roE0 = roe0 + 0.5 * ro0 * u0 * u0;
		double FE0 = u0 * (roE0 + p0);
		double roe_st = p_st / g1;
		double roE_st = roe_st + 0.5 * ro_st * u_st * u_st;
		double FE_st = u_st * (roE_st + p_st);
		dU = roE0 - roE_st;
		dF = FE0 - FE_st;
		double D3 = dF / dU;
		return (D1+ D2 + D3) / 3.0;
	};
	
	// collision of two media
	void RunExperiment(int argc, char *argv[]) {
		auto conf = DefaultSettings();
		auto params = DefaultStates();
		double rho_st, u_st, p_st;
		params.ref_velocity = ApproximateSolution(ShockDirection::right, rho_st, u_st, p_st);

		//init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
		};
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::WENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
		};
		kernel->Init(conf);

		// IC
		auto initD = [&conf, &params](Vector r) {
			params.x0 = 0.5;
			double p = params.P;
			double u = 0;
			double ro = 0;
			if (r.x < params.x0) {
				u = params.uL;
				ro = params.roL;
			} else {
				u = params.uR;
				ro = params.roR;
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
		kernel->SaveSolutionSega("init.dat");

		// Run computation
		kernel->Run();

		// compute exact solution

		// Finalize kernel
		kernel->Finalize();
	};
	
}	//end of namespace area


#endif
