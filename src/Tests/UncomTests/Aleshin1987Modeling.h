#ifndef TurboStructured_Tests_UncomTests_AleshinExp
#define TurboStructured_Tests_UncomTests_AleshinExp


#include <iostream>
#include <vector>
#include "math.h"
#include "kernel.h";
#include "Methods/ExplicitRungeKuttaFVM.h"
#include "Sensors/Sensors.h"

namespace AleshinExp {
	double PI = 3.14159265359;

	// struct for main parameters of that test
	struct Parameters {
		double gamma;		// specific heat ratio
		double amplitude;	// amplitude of perturbations a*sin (kx)
		double Ly;			// width of tube
		double Lx;			// length of the domain
		double xShock;		// initial shock position
		double xInterface;	// initial position of the contact discontinuity
		double roL;			// argon density
		double roR;			// xenon density
		double p0;			// initial pressure around interface
		double D;			// shock wave speed
		double uReff;		// velocity of frame of refference
	};

	// Default parameters
	void DefaultSettings(Parameters &par) {
		par.amplitude = 0.01;
		par.D = -600;			//shock moves in left direction
		par.gamma = 1.67;
		par.Lx = 0.08;
		par.xInterface = 0.015;
		par.xShock = 0.029;
		par.uReff = -400;		// frame of refference moves in left direction
		par.Ly = 0.5 * 7.2e-2;
		par.p0 = 5.0e4;			// a half of one barr
		par.roR = 1.784;		// argon
		par.roL = 5.894;		// xenon
	};

	// compute state after shock		
	void ComputeStateAfterShock(Parameters &par, double &p, double &ro, double &u ) {

		// compute velocity after the shock
		double D = par.D;
		double g = par.gamma;

		// compute velocity pressure and denscity after the shock
		u = -2.0 * (g * par.p0 - D * D * par.roR) / (D * par.roR * (1.0 + g));
		ro = par.roR * D / (D - u);
		p = par.p0 + ro * u * (D - u);

		// TO DO DELETE
		// Test for permutations of states 
		
		// Rankine - Hugoniot conditions cheking
		/*double ro0 = par.roR;
		double p0 = par.p0;
		double u0 = 0;
		double g1 = g - 1.0;
		
		double dU1 = ro0 - ro;
		double dF1 = ro0 * u0 - ro * u;
		dU1 *= D;
		
		double dU2 = ro0 * u0 - ro * u;
		double dF2 = p0 - p + ro0 * u0 * u0 - p - ro * u * u;
		dF2 /= D;

		double roe0 = p0 / g1;
		double roE0 = roe0 + 0.5 * ro0 * u0 * u0;
		double FE0 = u0 * (roE0 + p0);
		double roe = p / g1;
		double roE = roe + 0.5 * ro * u * u;
		double FE = u * (roE + p);
		double dU3 = roE0 - roE;
		double dF3 = FE0 - FE;
		dF3 /= D;*/

		return;
	};

	// Run one experiment ( parameters is as input data )
	void RunSingleExperiment(int modeNumber, double TotalTime, Parameters& par, int argc, char *argv[]) {
		KernelConfiguration conf;
		conf.nDims = 2;
		conf.nX = 400;
		conf.nY = 200;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.isPeriodicX = false;
		conf.isPeriodicY = false;
		conf.isUniformAlongX = true;
		conf.isUniformAlongY = true;

		conf.Gamma = par.gamma;

		conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
		conf.xLeftBoundary.Gamma = conf.Gamma;
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = conf.Gamma;
		conf.yLeftBoundary.BCType = BoundaryConditionType::SymmetryY;
		conf.yLeftBoundary.Gamma = conf.Gamma;
		conf.yRightBoundary.BCType = BoundaryConditionType::SymmetryY;
		conf.yRightBoundary.Gamma = conf.Gamma;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.methodConfiguration.CFL = 0.4;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.DummyLayerSize = 1;

		conf.MaxTime = TotalTime;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 5.0e-6;
		conf.SaveSolutionSnapshotIterations = 0;
		conf.ResidualOutputIterations = 10;

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

		double u3, ro3, p3;		// state after the shock
		ComputeStateAfterShock(par, p3, ro3, u3);

		auto initD = [&par, &conf, modeNumber, p3, ro3, u3](Vector r) {
			double ro;
			double p;
			double u;

			// domain is divided into three areas

		    // on the right of the shock
			if (r.x > par.xShock) {
				ro = ro3;
				p = p3;
				u = u3;
			}	// on the left of the shock   
			else {
				double xborder = par.xInterface;
				xborder -= par.amplitude * cos(2 * PI * r.y * modeNumber / par.Ly);

				// left part (xenon)
				if (r.x < xborder) {
					ro = par.roL;
					p = par.p0;
					u = 0;
				}
				// right part (argon)
				else {
					ro = par.roR;
					p = par.p0;
					u = 0;
				};
			};

			// compute ro_e and write conservative variables
			double roe = p / (par.gamma - 1.0);
			u -= par.uReff;
			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = 0.0;
			res[3] = 0.0;
			res[4] = roe + 0.5 * ro * u * u;		//total energy equals internal one because a motion is absent

			return res;
		};
		kernel->SetInitialConditions(initD);

		//save solution
		kernel->SaveSolution("init.dat");

		// Set sensors if needed
		auto GetInEnergy = [](std::valarray<double> vals) {
			double roe =  vals[4] - 0.5 * (vals[1] * vals[1] + vals[2] * vals[2]) / vals[0];
			return roe / vals[0];
		};

		kernel->isSensorEnable = true;
		kernel->SaveSensorRecordIterations = 1;

		// create a sensor
		std::unique_ptr<MValuePosXSensor2> sen1 = std::make_unique<MValuePosXSensor2>("border_pos.dat", *kernel->pManager, kernel->g, GetInEnergy);
		sen1->SetSensor((int)(0.5 * conf.nY / modeNumber + 1), 0, kernel->nVariables);
		kernel->Sensors.push_back(std::move(sen1));

		//run computation
		kernel->Run();

		//finalize kernel
		kernel->Finalize();
	};

	// Run Computation Experiment
	void RunExperiment(int argc, char *argv[]) {
		int modeNumber = 1;		// initial perturbation modes number
		double TotalTime = 40.0e-5;

		// Fill parameters structure (SI system)
		Parameters par;
		DefaultSettings(par);

		// Run experiments
		RunSingleExperiment(modeNumber, TotalTime, par, argc, argv);

		//end of experiments
		std::cout << "Aleshin's experiment simulation is completed";
		std::cout << std::endl;
		return;
	};

	// Additional part - 3D case
	// Run one experiment ( parameters is as input data)
	void Run3DExperiment(int argc, char *argv[]) {
		// experiment parameters
		int modeNumber = 2;			// initial perturbation modes number
		double TotalTime = 2.0e-4;

		// Fill parameters structure (SI system)
		Parameters par;
		DefaultSettings(par);

		KernelConfiguration conf;
		conf.nDims = 3;
		conf.nX = 160;
		conf.nY = 80;
		conf.nZ = 40;
		conf.LX = par.Lx;
		conf.LY = par.Ly;
		conf.LZ = conf.LY / modeNumber;
		conf.isPeriodicX = false;
		conf.isPeriodicY = true;
		conf.isPeriodicZ = true;
		conf.isUniformAlongX = true;
		conf.isUniformAlongY = true;
		conf.isUniformAlongZ = true;

		conf.Gamma = par.gamma;

		conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
		conf.xLeftBoundary.Gamma = conf.Gamma;
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = conf.Gamma;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;
		conf.methodConfiguration.CFL = 0.4;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.DummyLayerSize = 1;

		conf.MaxTime = TotalTime;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 1.0e-5;
		conf.SaveSolutionSnapshotIterations = 0;
		conf.ResidualOutputIterations = 10;

		// init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.SolutionMethod == KernelConfiguration::Method::ExplicitRungeKuttaFVM) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
			//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		kernel->Init(conf);

		double u3, ro3, p3;		// state after the shock
		ComputeStateAfterShock(par, p3, ro3, u3);

		auto initD = [&par, &conf, modeNumber, p3, ro3, u3](Vector r) {
			double ro;
			double p;
			double u;

			// domain is divided into three areas

			// on the right of the shock
			if (r.x > par.xShock) {
				ro = ro3;
				p = p3;
				u = u3;
			}	// on the left of the shock   
			else {
				double dr_z = abs(r.z - 0.5 * conf.LZ);		//dr - is a Vector (in YZ plane) from nearest center of perturbation to projection of r in YZ plane
				double lyambda = conf.LY / modeNumber;
				int n_mode = (int)(r.y / lyambda);
				double dr_y = abs(r.y - (n_mode + 0.5) * lyambda);
				double dr = sqrt(dr_y * dr_y + dr_z * dr_z);
				dr /= 0.5 * lyambda;
				double xborder = par.xInterface;
				if (dr >= 1.0) {
					xborder -= par.amplitude;
				}
				else {
					xborder -= par.amplitude * cos(PI * (1.0 - dr));
				};

				// left part (xenon)
				if (r.x < xborder) {
					ro = par.roL;
					p = par.p0;
					u = 0;
				}
				// right part (argon)
				else {
					ro = par.roR;
					p = par.p0;
					u = 0;
				};
			};

			// compute ro_e and write conservative variables
			double roe = p / (par.gamma - 1.0);
			u -= par.uReff;
			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = 0.0;
			res[3] = 0.0;
			res[4] = roe + 0.5 * ro * u * u;		//total energy equals internal one because a motion is absent

			return res;
		};
		kernel->SetInitialConditions(initD);

		//save solution
		kernel->SaveSolution("init.dat");

		// Set sensors if needed
		kernel->isSensorEnable = true;
		kernel->SaveSensorRecordIterations = 1;

		// Init target function
		auto GetInEnergy = [](std::valarray<double> vals) {
			double roe = vals[4] - 0.5 * (vals[1] * vals[1] + vals[2] * vals[2]) / vals[0];
			return roe / vals[0];
		};
		
		// Create a sensor
		std::unique_ptr<MValuePosXSensor2> sen1 = std::make_unique<MValuePosXSensor2>("border_pos.dat", *kernel->pManager, kernel->g, GetInEnergy);
		sen1->SetSensor((int)(0.5 * conf.nY / modeNumber + 1), (int)(0.5 * conf.nZ + 1), kernel->nVariables);
		kernel->Sensors.push_back(std::move(sen1));

		//run computation
		kernel->Run();

		//finalize kernel
		kernel->Finalize();

		//end of experiments
		std::cout << "Aleshin's experiment 3D simulation is completed";
		std::cout << std::endl;
	};
}



#endif