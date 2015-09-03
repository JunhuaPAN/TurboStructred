#ifndef TurboStructured_Tests_1DTests_ToroTest
#define TurboStructured_Tests_1DTests_ToroTest

#include <iostream>
#include <vector>
#include "kernel\kernel.h";
#include "Methods\ExplicitRungeKuttaFVM.h"

namespace ToroTests
{

// some usefull structures
struct ShockTubeParameters {
	double gamma;
	double roL;
	double PL;
	double uL;
	double roR;
	double PR;
	double uR;
};

// Error norm functions
double ComputeL2Error(std::valarray<double>& comp_val, std::valarray<double>& exac_val, int nV) {
	return 0;
}

double ComputeLinfError(std::valarray<double>& comp_val, std::valarray<double>& exac_val, int nV) {
	return 0;
}


// class for test frob Toro book
class ToroTest
{
public:
	double x0;
	ShockTubeParameters par;

	// Initial states function depending on test number
	void SetStates(int TestNumber, double gamma) {
		assert(TestNumber != 2, "Toro test #2 isn't implemented");
		if (TestNumber == 1) {
			par.roL = 1.0;
			par.uL = 0.0;
			par.PL = 1.0;
			par.gamma = gamma;
			par.roR = 0.125;
			par.uR = 0.0;
			par.PR = 0.1;
		};

		return;
	};

	// Compute exact solution depending on test number
	/*std::vector<double> ComputeExactSolutionInCell(int TestNumber, double x, double t) {
		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t
		double S = x / t;
		double ro;
		double u;
		double p;

		if (starValues.uStar >= S) {
			//Left side of contact
			//Shock wave
			if (starValues.leftWave == Godunov3DSolverPerfectGas::Shock) {
				if (starValues.SL >= S) {
					//Case a1
					//Supersonic shock
					ro = roL;
					u = uL;
					p = pL;
				}
				else {
					//Case a2
					//Subsonic shock
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.leftWave == Godunov3DSolverPerfectGas::Rarefaction) {
				if (starValues.SHL > S) {
					//Left region
					ro = roL;
					u = uL;
					p = pL;
				}
				else if (S > starValues.STL) {
					//Star region
					ro = starValues.roStarL;
					u = starValues.uStar;
					p = starValues.pStar;
				}
				else {
					//Rarefaction fan region
					double aL = sqrt(gammaL * pL / roL);
					double CL = 2 / (gammaL + 1) + (gammaL - 1) * (uL - S) / ((gammaL + 1) * aL);
					CL = pow(CL, 2 / (gammaL - 1));

					//Density
					ro = CL * roL;

					//Velocity
					u = aL + (gammaL - 1) * uL / 2.0 + S;
					u *= 2 / (gammaL + 1);

					//Pressure					
					p = CL * pL;
				};
			};

		}
		else {
			//Right side of contact

			//Shock wave
			if (starValues.rightWave == Godunov3DSolverPerfectGas::Shock) {
				if (starValues.SR <= S) {
					//Case a1
					//Supersonic shock
					ro = roR;
					u = uR;
					p = pR;
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
			if (starValues.rightWave == Godunov3DSolverPerfectGas::Rarefaction) {
				if (starValues.SHR < S) {
					//Right region
					ro = roR;
					u = uR;
					p = pR;
				}
				else if (S < starValues.STR) {
					//Star region
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				}
				else {
					//Rarefaction fan region
					double aR = sqrt(gammaR * pR / roR);
					double CR = 2 / (gammaR + 1) - (gammaR - 1) * (uR - S) / ((gammaR + 1) * aR);
					CR = pow(CR, 2 / (gammaR - 1));

					//Density
					ro = CR * roR;

					//Velocity
					u = -aR + (gammaR - 1) * uR / 2.0 + S;
					u *= 2 / (gammaR + 1);

					//Pressure					
					p = CR * pR;
				};
			};
		};

		std::vector<double> result;
		result.push_back(ro);
		result.push_back(u);
		result.push_back(p);
		return result;
	};*/

	std::pair<double, double> Run(int argc, char *argv[]) {
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
			//kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
			kernel = std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
		};
		kernel->Init(conf);

		// initial conditions
		SetStates(1, conf.Gamma);

		//save solution
		kernel->SaveSolution("init.dat");

		//run computation
		kernel->Run();

		//finalize kernel
		kernel->Finalize();

		return{};
	};
};

}	//end of namespace area


#endif
