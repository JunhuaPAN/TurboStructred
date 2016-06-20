#ifndef TurboStructured_Tests_1DTests_RarefactionTest2
#define TurboStructured_Tests_1DTests_RarefactionTest2

#include <iostream>
#include <vector>
#include "kernel\kernel.h"
#include "Methods\ExplicitRungeKuttaFVM.h"
#include "TestsUtility.h"

//	The base of this test is Test #3 from E.Toro book (p. 129)					//
//	In that Riemann problem we have one rarefaction wave folowing to the left	//
//	We consider only inner part of the rarefaction for order estimation			//

namespace RarefactionTest2
{
	// special variable for rarefaction segment length
	double subgrid_length;

	// Parameters structure
	struct ShockTubeParameters {
		double gamma;
		double roL;
		double PL;
		double uL;
		double roR;
		double PR;
		double uR;
		double x0;
	};

	//Types of nonlinear waves
	enum class WaveType {
		Shock,
		Rarefaction
	};

	//Find star region quantites
	struct StarVariables {
		double pStar;
		double uStar;
		double roStarL;
		double roStarR;
		WaveType leftWave;
		WaveType rightWave;

		//Propagation speeds
		double SL;
		double SHL;
		double STL;

		double SR;
		double SHR;
		double STR;

		double MaxSpeed;
	} starValues;

	// Default settings and default state
	KernelConfiguration DefaultSettings() {
		KernelConfiguration conf;

		// Values by default
		conf.nDims = 1;
		conf.nX = 200;
		conf.LX = 1.0;
		conf.isPeriodicX = false;
		conf.isUniformAlongX = true;
		conf.qx = 1.00;
		conf.Gamma = TestsUtility::gamma1 + 1;

		conf.SolutionMethod = KernelConfiguration::Method::ExplicitRungeKuttaFVM;

		conf.DummyLayerSize = 1;
		conf.methodConfiguration.CFL = 0.45;
		conf.methodConfiguration.RungeKuttaOrder = 1;
		conf.methodConfiguration.Eps = 0.05;
		conf.methodConfiguration.ReconstructionType = Reconstruction::PiecewiseConstant;
		conf.methodConfiguration.RiemannProblemSolver = RPSolver::RoePikeSolver;

		conf.xLeftBoundary.BCType = BoundaryConditionType::Natural;
		conf.xLeftBoundary.Gamma = 1.4;
		conf.xRightBoundary.BCType = BoundaryConditionType::Natural;
		conf.xRightBoundary.Gamma = 1.4;

		conf.MaxTime = 0.25;
		conf.MaxIteration = 100000;
		conf.SaveSolutionTime = 0;
		conf.SaveSolutionIterations = 0;
		conf.ResidualOutputIterations = 100;

		return conf;
	};

	ShockTubeParameters DefaultState() {
		ShockTubeParameters params;
		params.gamma = TestsUtility::gamma1 + 1;
		params.roL = 1.0;
		params.roR = 0.125;
		params.PL = 1.0;
		params.uL = 0.0;
		params.uR = 0.0;
		params.PR = 0.1;
		params.x0 = 0.5;		// initial discontinuity position
		return params;
	};

	// Create a subgrid that contain cells with rarefaction fan
	Grid GetRarefactionCells(StarVariables& starV, Grid& g, double time, double x0, double Lx) {
		Grid res;
		double cut_factor = 0.12;	// to cut both end of rarefaction fan
											// define head and tail speeds
		double velL = starV.SHL;
		double velR = starV.STL;

		// define head and tail coordinates
		double xl = x0 + velL * time;
		double xr = x0 + velR * time;

		// define indexes of ends of rarefaction cell set
		double deltax = xr - xl;
		xl -= cut_factor * deltax;
		xr += cut_factor * deltax;
		int iMin = (int)(xl * g.nX / Lx) + g.dummyCellLayersX;
		int iMax = (int)(xr * g.nX / Lx) + g.dummyCellLayersX;

		// store segment length
		subgrid_length = g.CoordinateX[iMax] - g.CoordinateX[iMin] + 0.5 * (g.hx[iMax] + g.hx[iMin]);

		// create subgrid for each processors
		res = CreateSubGrid(iMin, iMax, g.jMin, g.jMax, g.kMin, g.kMax, g);
		return std::move(res);
	};

	// Choose star values for N test
	void SetStarValues(ShockTubeParameters& params) {

		// Write approximate solution
		starValues.leftWave = WaveType::Rarefaction;
		starValues.rightWave = WaveType::Shock;
		starValues.pStar = 0.30313;
		starValues.uStar = 0.92745;
		starValues.roStarL = 0.42632;
		starValues.roStarR = 0.26557;

		// compute speeds of three waves
		double gamma = params.gamma;

		//Left side of contact
		double pRatioL = starValues.pStar / params.PL;
		double aL = sqrt(gamma * params.PL / params.roL);

		//Left rarefaction head
		starValues.SHL = params.uL - aL;

		//Determine rarefaction tail propagation speed		
		double aStarL = aL * pow(pRatioL, (gamma - 1) / (2 * gamma));
		starValues.STL = starValues.uStar - aStarL;

		//Right side of contact
		double pRatioR = starValues.pStar / params.PR;
		double aR = sqrt(gamma * params.PR / params.roR);

		//Determine shock propagation speed			
		starValues.SR = sqrt((gamma + 1) * pRatioR / (2 * gamma) + (gamma - 1) / (2 * gamma));
		starValues.SR = params.uR + aR * starValues.SR;
	};

	// Compute exact solution in Cell
	std::vector<double> ComputeExactSolutionInCell(ShockTubeParameters& pars, double x, double t) {
		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t
		double S = x / t;
		double ro;
		double u;
		double p;

		if (starValues.uStar >= S) {
			//Left side of contact
			//Shock wave
			if (starValues.leftWave == WaveType::Shock) {
				if (starValues.SL >= S) {
					// Case a1
					// Left of the shock
					ro = pars.roL;
					u = pars.uL;
					p = pars.PL;
				}
				else {
					// Case a2
					// Right of the shock shock
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.leftWave == WaveType::Rarefaction) {
				if (starValues.SHL > S) {
					//Left region
					ro = pars.roL;
					u = pars.uL;
					p = pars.PL;
				}
				else if (S > starValues.STL) {
					//Star region
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				}
				else {
					//Rarefaction fan region
					double aL = sqrt(pars.gamma * pars.PL / pars.roL);
					double CL = 2.0 / (pars.gamma + 1) + (pars.gamma - 1) * (pars.uL - S) / ((pars.gamma + 1) * aL);
					double CLRo = pow(CL, 2.0 / (pars.gamma - 1));
					double CLP = pow(CL, 2.0 * pars.gamma / (pars.gamma - 1));

					//Density
					ro = CLRo * pars.roL;

					//Velocity
					u = aL + 0.5 * (pars.gamma - 1) * pars.uL + S;
					u *= 2 / (pars.gamma + 1);

					//Pressure					
					p = CLP * pars.PL;
				};
			};

		}
		else {
			//Right side of contact

			//Shock wave
			if (starValues.rightWave == WaveType::Shock) {
				if (starValues.SR <= S) {
					//Case a1
					//Supersonic shock
					ro = pars.roR;
					u = pars.uR;
					p = pars.PR;
				}
				else {
					//Case a2
					//Subsonic shock
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.rightWave == WaveType::Rarefaction) {
				if (starValues.SHR < S) {
					//Right region
					ro = pars.roR;
					u = pars.uR;
					p = pars.PR;
				}
				else if (S < starValues.STR) {
					//Star region
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				}
				else {
					//Rarefaction fan region
					double aR = sqrt(pars.gamma * pars.PR / pars.roR);
					double CR = 2.0 / (pars.gamma + 1) - (pars.gamma - 1) * (pars.uR - S) / ((pars.gamma + 1) * aR);
					double CRRo = pow(CR, 2.0 / (pars.gamma - 1));
					double CRP = pow(CR, 2.0 * pars.gamma / (pars.gamma - 1));

					//Density
					ro = CRRo * pars.roR;

					//Velocity
					u = -aR + (pars.gamma - 1) * pars.uR / 2.0 + S;
					u *= 2 / (pars.gamma + 1);

					//Pressure					
					p = CRP * pars.PR;
				};
			};
		};

		std::vector<double> result;
		result.push_back(ro);
		result.push_back(u);
		result.push_back(p);
		return result;
	};

	// Compute exact solution for test N
	std::valarray<double> ComputeExactSolution(ShockTubeParameters& pars, Grid& g, double t) {
		// init result valarray
		std::valarray<double> res(g.nlocalX * TestsUtility::nVar);

		// Compute star state from Toro book p. 133
		SetStarValues(pars);

		// Compute exact solution in every cell and write the result valarray
		for (int i = 0; i < g.nlocalX; i++) {
			double x = g.CoordinateX[g.iMin + i] - pars.x0;

			// second order for integrall
			double xl, xr;
			const double d2 = 0.5773502691896257;		// Legandr delta
			double delta = d2 * 0.5 * g.hx[g.iMin + i];

			// convert primitive variables to conservative ones
			auto convert = [](const std::vector<double>& arg) {
				double rho = arg[0];
				double u = arg[1];
				double p = arg[2];
				double e = p / (TestsUtility::gamma1 * rho);
				return std::vector<double> { rho, rho * u, 0, 0, rho * e + 0.5 * rho * u * u };
			};

			// compute conservative variables in Legandr points
			std::vector<double> vars1 = convert(ComputeExactSolutionInCell(pars, x - delta, t));
			std::vector<double> vars2 = convert(ComputeExactSolutionInCell(pars, x + delta, t));

			// cell averaged conservative variables ( 2nd order )
			for (int j = 0; j < TestsUtility::nVar; j++) res[TestsUtility::nVar * i + j] = 0.5 * (1.0 * vars1[j] + 1.0 * vars2[j]);
		};
		return res;
	};

	// Collision of two media
	std::vector<double> RunSingleExperiment(int argc, char *argv[], int Nx, Reconstruction RecType, RPSolver _RPsolver) {
		// Result
		std::vector<double> errors;

		// Use default settings
		KernelConfiguration conf = DefaultSettings();
		ShockTubeParameters params = DefaultState();
		conf.nX = Nx;
		conf.methodConfiguration.ReconstructionType = RecType;
		conf.methodConfiguration.RiemannProblemSolver = _RPsolver;

		// Solution file
		std::ostringstream fname;

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
		kernel->Init(conf);

		// IC
		auto initD = [&conf, &params](Vector r) {
			double p = params.PL;
			double u = params.uL;
			double ro = params.roL;
			if (r.x > params.x0) {
				p = params.PR;
				u = params.uR;
				ro = params.roR;
			};

			std::vector<double> res(5);
			res[0] = ro;
			res[1] = ro * u;
			res[2] = 0;
			res[3] = 0;
			res[4] = p / TestsUtility::gamma1 + 0.5 * ro * u * u;
			return res;
		};
		kernel->SetInitialConditions(initD);

		// Run computation
		kernel->Run();

		// Compute exact solution in rarefaction part
		SetStarValues(params);
		Grid rarefGrid = GetRarefactionCells(starValues, kernel->grid, kernel->stepInfo.Time, params.x0, conf.LX);
		TestsUtility::exact_solution = ComputeExactSolution(params, rarefGrid, kernel->stepInfo.Time);

		// Compute accuracy
		int local_iMin = kernel->grid.dummyCellLayersX + rarefGrid.iMin - kernel->grid.iMin;
		std::valarray<double> inner_values(&(kernel->values[TestsUtility::nVar * local_iMin]), rarefGrid.nlocalX * kernel->nVariables);
		std::vector<double> L2 = TestsUtility::ComputeL2Error(inner_values, rarefGrid, *(kernel->pManager), subgrid_length);
		std::vector<double> L1 = TestsUtility::ComputeL1Error(inner_values, rarefGrid, *(kernel->pManager), subgrid_length);

		// Write in errors
		errors = L2;
		for (int i = 0; i < TestsUtility::nVar; i++) errors.push_back(L1[i]);

		// Show the errors
		if (kernel->pManager->IsMaster()) {
			std::cout << "Rarefaction test has finished. ";
			std::cout << std::endl;
			std::cout << "ReconstructionType: " << fname.str();
			std::cout << ", Nx = " << Nx << std::endl;
			std::cout << "L2_rho = " << L2[0];
			std::cout << ", L1_rho = " << L1[0];
			std::cout << std::endl;
		};

		// Save both solution in TecPlot
		fname << ", Nx=";
		fname << Nx;
		fname << ", t=";
		fname << conf.MaxTime;
		fname << ".dat";
		std::string filename = fname.str();
		kernel->SaveSolution(filename);
		TestsUtility::exact_solution = ComputeExactSolution(params, kernel->grid, kernel->stepInfo.Time);
		TestsUtility::SaveExactSolution(filename, kernel->grid, *(kernel->pManager));

		// Save the errors
		std::stringstream fname2;
		fname2 << "Rarefaction_test2.dat";
		TestsUtility::WriteErrors(fname2.str(), errors, kernel->grid, *(kernel->pManager));

		// Finalize kernel
		kernel->Finalize();

		return errors;
	};

	void RunExperiment(int argc, char *argv[]) {
		int Nx = 800;

		// Reconstruction type
		// Reconstruction RecType{ Reconstruction::PiecewiseConstant };
		Reconstruction RecType{ Reconstruction::ENO2PointsStencil };
		
		// RP solver
		RPSolver rSolver{ RPSolver::GodunovSolver };
		// RPSolver rSolver{};

		// collect all errors in file
		std::vector<double> err;
		err = RunSingleExperiment(argc, argv, Nx, RecType, rSolver);

		return;
	};

}	//end of namespace area


#endif
