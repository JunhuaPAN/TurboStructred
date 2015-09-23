#ifndef TurboStructured_Tests_1DTests_ContactDisTest
#define TurboStructured_Tests_1DTests_ContactDisTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h"
#include "Methods\ExplicitRungeKuttaFVM.h"

// single contact discontinuity test
// roL = 1.0 || roR = 0.125
// uL = 1.0  || uR = 1.0
// pL = 1.0  || pR = 1.0

namespace ContactDisTest
{
	const int nVar = 5;
	const double gamma1 = 0.4;
	std::valarray<double> exact_solution;

	// Error norm functions
	std::vector<double> ComputeL2Error(std::valarray<double>& comp_val, Grid& g, ParallelManager& par) {
		std::vector<double> res_temp(nVar, 0);
		std::vector<double> res(nVar, 0);
		for (int i = 0; i < comp_val.size(); i++) {
			double err = comp_val[i] - exact_solution[i];
			res_temp[i % nVar] += (err * err) * g.hx[g.iMin + i / nVar];
		};
		par.Barrier();
		for (int i = 0; i < nVar; i++) res[i] = par.Sum(res_temp[i]);
		for (int i = 0; i < nVar; i++) res[i] = sqrt(res[i]);
		return res;
	}

	std::vector<double> ComputeLinfError(std::valarray<double>& comp_val, Grid& g, ParallelManager& par) {
		std::vector<double> res_temp(nVar, 0);
		std::vector<double> res(nVar, 0);
		for (int i = 0; i < comp_val.size(); i++) {
			double err = comp_val[i] - exact_solution[i];
			if (res_temp[i % nVar] < abs(err)) res_temp[i % nVar] = abs(err);
		};
		par.Barrier();
		for (int i = 0; i < nVar; i++) res[i] = par.Max(res_temp[i]);

		return res;
	}

	// Save both solutions
	void SaveExactSolution(std::string fname, Grid& g, ParallelManager& pManager) {
		//Tecplot version    
		int rank = pManager.getRank();
		
		// Open the file
		if (!pManager.IsFirstNode()) pManager.Wait(rank - 1);

		//Reopen file for writing
		std::ofstream ofs(fname, std::ios_base::app);
		ofs << R"(ZONE T="Exact solution")";
		ofs << std::endl;
		ofs << std::scientific;

		// Exact solution
		for (int i = 0; i < g.nlocalX; i++) {
			//Obtain cell data
			double rho = exact_solution[i * nVar];
			double u = exact_solution[i * nVar + 1] / rho;
			double v = 0;
			double w = 0;
			double e = exact_solution[i * nVar + 4] / rho - 0.5 * u * u;
			double P = gamma1 * rho * e;

			//Write to file
			ofs << g.CoordinateX[g.iMin + i] << " ";
			ofs << rho << " ";
			ofs << u << " ";
			ofs << v << " ";
			ofs << w << " ";
			ofs << P << " ";
			ofs << e << " ";
			ofs << std::endl;
		};	//	end cycle

			//Signal to next process to begin writing
		if (!pManager.IsLastNode()) {
			pManager.Signal(rank + 1);
		};

		//Syncronize
		pManager.Barrier();
		return;
	};

	// IC structure
	struct ShockTubeParameters {
		double roL;
		double roR;
		double P;
		double uL;
		double uR;
		double x0;
		double ref_velocity;
	};

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
		conf.Gamma = gamma1 + 1;

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

		conf.MaxTime = 0.2;
		conf.MaxIteration = 1000000;
		conf.SaveSolutionSnapshotTime = 0;
		conf.SaveSolutionSnapshotIterations = 0;
		conf.ResidualOutputIterations = 100;

		return conf;
	};

	ShockTubeParameters DefaultStates() {
		ShockTubeParameters params;
		params.roL = 1.0;
		params.roR = 0.125;
		params.P = 1.0;
		params.uL = 1.0;
		params.uR = 1.0;
		params.x0 = 0.2;		// initial discontinuity position
		params.ref_velocity = 0;

		return params;
	};

	std::valarray<double> ComputeExactSolution(ShockTubeParameters& pars, Grid& g, double time) {
		double x_disc = pars.x0 + 1.0 * time;		// uL = uR = 1
		std::valarray<double> res(g.nlocalX * nVar);
		double u = 1.0;
		double p = pars.P;
		for (int i = 0; i < g.nlocalX; i++) {
			double rho = pars.roR;
			if (g.CoordinateX[g.iMin + i] <= x_disc) rho = pars.roL;
			double e = p / (gamma1 * rho);
			res[i * nVar] = rho;
			res[i * nVar + 1] = rho*u;
			res[i * nVar + 2] = 0;
			res[i * nVar + 3] = 0;
			res[i * nVar + 4] = rho * e + 0.5 * rho * u * u;
		};
		return res;
	};

	// Collision of two media
	std::vector<double> RunSingleExperiment(int argc, char *argv[], int Nx, double MaxTime, KernelConfiguration::Reconstruction RecType) {
		// Result
		std::vector<double> errors;

		// Use default settings
		KernelConfiguration conf = DefaultSettings();
		ShockTubeParameters params = DefaultStates();
		conf.nX = Nx;
		conf.MaxTime = MaxTime;
		conf.ReconstructionType = RecType;

		// Solution file
		std::ostringstream fname;

		// Init kernel
		std::unique_ptr<Kernel> kernel;
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::PiecewiseConstant) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
			fname << "PWConstant";
		};
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::ENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
			fname << "ENO2";
		};
		if (conf.ReconstructionType == KernelConfiguration::Reconstruction::WENO2PointsStencil) {
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<WENO2PointsStencil>(&argc, &argv));
			fname << "WENO2";
		};
		kernel->Init(conf);

		// IC
		auto initD = [&conf, &params](Vector r) {
			double p = params.P;
			double u = 0;
			double ro = 0;
			if (r.x < params.x0) {
				u = params.uL;
				ro = params.roL;
			}
			else {
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
		kernel->SaveSolution("init.dat");

		// Run computation
		kernel->Run();

		// Compute exact solution
		exact_solution = ComputeExactSolution(params, kernel->g, conf.MaxTime);

		// Compute accuracy
		std::valarray<double> inner_values(&(kernel->values[nVar * kernel->g.dummyCellLayersX]), kernel->g.nlocalX * kernel->nVariables);
		std::vector<double> L2 = ComputeL2Error(inner_values, kernel->g, *(kernel->pManager));
		std::vector<double> Linf = ComputeLinfError(inner_values, kernel->g, *(kernel->pManager));

		// Write in errors
		errors = L2;
		for (int i = 0; i < nVar; i++) errors.push_back(Linf[i]);

		// Show the errors
		if (kernel->pManager->IsMaster()) {
			std::cout << "Contact discontinuity test was finished. ";
			std::cout << std::endl;
			std::cout << "ReconstructionType: " << fname.str();
			std::cout << ", Nx = " << Nx << std::endl;
			std::cout << "L2_rho = " << L2[0];
			std::cout << ", Linf_rho = " << Linf[0];
			std::cout << std::endl;
		};
		
		// Save both solution in TecPlot
		fname << ", Nx=";
		fname << Nx;
		fname << ", t=";
		fname << MaxTime;
		fname << ".dat";
		std::string filename = fname.str();
		kernel->SaveSolution(filename);
		SaveExactSolution(filename, kernel->g, *(kernel->pManager));

		// Finalize kernel
		kernel->Finalize();

		return errors;
	};

	void RunExperiment(int argc, char *argv[]) {
		int Nx = 1600;
		double MaxTime = 0.2;
		//KernelConfiguration::Reconstruction RecType = KernelConfiguration::Reconstruction::PiecewiseConstant;
		//KernelConfiguration::Reconstruction RecType = KernelConfiguration::Reconstruction::ENO2PointsStencil;
		KernelConfiguration::Reconstruction RecType = KernelConfiguration::Reconstruction::WENO2PointsStencil;

		// collect all errors in file
		std::vector<double> err;
		//ofs << R"(VARIABLES = "Nx" "L2_rho" "L2_rhou" "L2_rhoE" "Linf_rho" "Linf_rhou" "Linf_rhoE")" << std::endl;
		err = RunSingleExperiment(argc, argv, Nx, MaxTime, RecType);

		// write result in file
		std::string fname("Test-ContactDiscontinuity.dat");
		std::ofstream ofs(fname, std::ios_base::app);
		ofs << Nx << ' ';
		ofs << err[0] << ' ';
		ofs << err[1] << ' ';
		ofs << err[4] << ' ';
		ofs << err[nVar] << ' ';
		ofs << err[nVar + 1] << ' ';
		ofs << err[nVar + 4] << ' ';
		ofs << std::endl;
		ofs.close();
	};

}	//end of namespace area


#endif
