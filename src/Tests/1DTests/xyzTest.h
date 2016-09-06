#ifndef TurboStructured_Tests_1DTests_xyzTest
#define TurboStructured_Tests_1DTests_xyzTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h"
#include "Methods\ExplicitRungeKuttaFVM.h"

// 1D SOD test in X Y and Z directions (3D representation)

namespace SODxyzTest
{
	// All what we need in these tests
	struct ShockTubeParameters {
		double roL{ 1.0 };
		double roR{ 0.125 };
		double PL{ 1.0 };
		double PR{ 0.1 };
		double uL{ 0.0 };
		double uR{ 0.0 };
		double r0{ 0.5 };
	} params;	  // IC structure (velocity is zero)
	KernelConfiguration conf;		// experiment configurator
	std::ostringstream fname;		// Solution file
	Slice my_slice;					// one slice for solution
	Direction my_dir;				// choosen direction
	
	// Default settings
	void DefaultSettings() {
		// Values by default
		conf.nDims = 3;
		conf.nX = 1;
		conf.nY = 1;
		conf.nZ = 1;
		conf.LX = 1.0;
		conf.LY = 1.0;
		conf.LZ = 1.0;
		conf.Gamma = 1.4;

		// BC
		conf.MyConditions[1] = BoundaryConditionConfiguration(BoundaryConditionType::Natural);

		// Method configuration
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::Linear2psLim;
		conf.DummyLayerSize = 1;
		conf.methodConfiguration.CFL = 0.3;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;

		// Computation parameters
		conf.MaxTime = 0.25;
		conf.MaxIteration = 10000;
		conf.SaveSliceTime = 0.25;
		conf.ResidualOutputIterations = 100;
		fname << "SODtest";
	};

	// Functions which depend on the direction of test
	void ApplyDirection(int Nc) {
		if (my_dir == Direction::XDirection) {
			conf.nX = Nc;
			conf.isPeriodicX = false;

			// Grid compression
			BlockNode nleft, ncenter;
			nleft.N_cells = conf.nX / 2;
			nleft.q_com = 1.0 / 1.05;
			ncenter.pos = 0.5 * conf.LX;
			ncenter.N_cells = conf.nX - nleft.N_cells;
			ncenter.q_com = 1.0 / nleft.q_com;
			conf.CompressionX[0] = nleft;
			conf.CompressionX.push_back(ncenter);

			// BC
			conf.xLeftBoundary.SetMarker(1);
			conf.xRightBoundary.SetMarker(1);

			// filename and slice
			fname << "Xdir";
			my_slice = Slice(-1, 1, 1);
		};
		if (my_dir == Direction::YDirection) {
			conf.nY = Nc;
			conf.isPeriodicY = false;

			// Grid compression
			BlockNode nleft, ncenter;
			nleft.N_cells = conf.nY / 2;
			nleft.q_com = 1.0 / 1.05;
			ncenter.pos = 0.5 * conf.LY;
			ncenter.N_cells = conf.nY - nleft.N_cells;
			ncenter.q_com = 1.0 / nleft.q_com;
			conf.CompressionY[0] = nleft;
			conf.CompressionY.push_back(ncenter);

			// BC
			conf.yLeftBoundary.SetMarker(1);
			conf.yRightBoundary.SetMarker(1);

			// filename and slice
			fname << "Ydir";
			my_slice = Slice(1, -1, 1);
		};
		if (my_dir == Direction::ZDirection) {
			conf.nZ = Nc;
			conf.isPeriodicZ = false;

			// Grid compression
			BlockNode nleft, ncenter;
			nleft.N_cells = conf.nZ / 2;
			nleft.q_com = 1.0 / 1.05;
			ncenter.pos = 0.5 * conf.LZ;
			ncenter.N_cells = conf.nZ - nleft.N_cells;
			ncenter.q_com = 1.0 / nleft.q_com;
			conf.CompressionZ[0] = nleft;
			conf.CompressionZ.push_back(ncenter);

			// BC
			conf.zLeftBoundary.SetMarker(1);
			conf.zRightBoundary.SetMarker(1);

			// filename and slice
			fname << "Zdir";
			my_slice = Slice(1, 1, -1);
		};
	};
	bool isLeftState(Vector r) {
		if (my_dir == Direction::XDirection) return (r.x < params.r0);
		if (my_dir == Direction::YDirection) return (r.y < params.r0);
		if (my_dir == Direction::ZDirection) return (r.z < params.r0);
	};

	// Run Experiment
	void RunSingleExperiment(int argc, char *argv[]) {
		// Use default settings
		DefaultSettings();
		
		// Grid size
		auto Nc = 100;
		ApplyDirection(Nc);

		// Init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
			fname << "PWConstant";
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
			fname << "ENO2";
		};
		if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2psLim) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<Linear2psLim<limBarsJespersen> >(&argc, &argv));
			fname << "LinearRecLim";
		};
		kernel->Init(conf);

		// IC
		auto initD = [](Vector r) {
			double p, ro;
			if (isLeftState(r)) {
				ro = params.roL;
				p = params.PL;
			} else {
				ro = params.roR;
				p = params.PR;
			};

			std::vector<double> res(5, 0);
			res[0] = ro;
			res[4] = p / (conf.Gamma - 1.0);
			return res;
		};
		kernel->SetInitialConditions(initD);

		// push slice
		kernel->slices.push_back(my_slice);

		// Run computation
		kernel->Run();

		// Finalize kernel
		kernel->Finalize();
	};

	void RunExperiment(int argc, char *argv[], Direction dir) {
		
		// choose the direction and run the test
		my_dir = dir;
		RunSingleExperiment(argc, argv);
		return;
	};

}	//end of namespace area


#endif
