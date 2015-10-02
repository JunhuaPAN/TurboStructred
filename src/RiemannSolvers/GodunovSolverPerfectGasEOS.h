#ifndef TurboStructured_RiemannSolvers_GodunovSolverPerfectGasEOS
#define TurboStructured_RiemannSolvers_GodunovSolverPerfectGasEOS

#include "RiemannSolver.h"
#include "alglibinternal.h"
#include "optimization.h"


class GodunovSolverPerfectGasEOS : public RiemannSolver {
	// Number of variables
	static const int nVariables = 5;

	// Required info
	double _gamma;	//	Specific heat ratio (inherited from gas model)
	double _operatingPressure; //	Operating pressure (optional, 0 by default)

	// Get Pressure function
	double GetPressure(const std::valarray<double>& celldata) {
		double ro = celldata[0];
		double vx = celldata[1] / ro;
		double vy = celldata[2] / ro;
		double vz = celldata[3] / ro;
		double E = celldata[4] / ro;
		double P = (_gamma - 1.0) * ro * (E - (vx * vx + vy * vy + vz * vz) / 2.0);
		return P;
	};

	//	Numerical flux
	std::vector<double> F(std::valarray<double> &U, Vector n)
	{
		std::vector<double> res(nVariables, 0);
		double ro = U[0];
		double vx = U[1] / ro;
		double vy = U[2] / ro;
		double vz = U[3] / ro;
		double roE = U[4];	//ro*e
		double p = (_gamma - 1.0) * (roE - ro * (vx * vx + vy * vy + vz * vz) / 2.0) - _operatingPressure;
		double vn = vx * n.x + vy * n.y + vz * n.z;

		res[0] = ro * vn;
		res[1] = ro * vn * vx + n.x * p;
		res[2] = ro * vn * vy + n.y * p;
		res[3] = ro * vn * vz + n.z * p;
		res[4] = vn * (roE + p);
		return res;
	};
	
	// usefull structure
	struct fParameters {
		double roL;
		double pL;
		double uL;
		double roR;
		double pR;
		double uR;
	} params;

	// fL part of algebraic equation for pressure in exact RP solver		(proposition 4.2.1 from Toro)
	double f1(double roL, double uL, double pL, double pStar)
	{
		double res = 0;

		if (pStar > pL) {
			//Left shock
			double AL = 2 / ((_gamma + 1.0) * roL);
			double BL = pL * (_gamma - 1.0) / (_gamma + 1.0);
			res = (pStar - pL) * sqrt(AL / (pStar + BL));
		}
		else {
			//Left rarefaction
			double aL = sqrt(_gamma * pL / roL);
			res = pow(pStar / pL, (_gamma - 1.0) / (2 * _gamma)) - 1.0;
			res *= 2 * aL / (_gamma - 1.0);
		};

		return res;
	};

	// fR part
	double f2(double roR, double uR, double pR, double pStar)
	{
		double res = 0;

		if (pStar > pR) {
			//Right shock
			double AR = 2 / ((_gamma + 1.0) * roR);
			double BR = pR * (_gamma - 1.0) / (_gamma + 1.0);
			res = (pStar - pR) * sqrt(AR / (pStar + BR));
		}
		else {
			//Right rarefaction
			double aR = sqrt(_gamma * pR / roR);
			res = pow(pStar / pR, (_gamma - 1.0) / (2 * _gamma)) - 1.0;
			res *= 2 * aR / (_gamma - 1.0);
		};

		return res;
	};

	// Target function for pressure equation in exact RP solver
	static void fFunction(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *object)
	{
		GodunovSolverPerfectGasEOS* solver = (GodunovSolverPerfectGasEOS*) object;
		//
		// this callback calculates		
		// f(pStar) = f1(pStar, WL) + f2(pStar, WR) + deltaU
		//
		fParameters par = solver->params;
		double pStar = x[0];
		double deltaU = par.uR - par.uL;
		double res = solver->f1(par.roL, par.uL, par.pL, pStar) + solver->f2(par.roR, par.uR, par.pR, pStar) + deltaU;

		fi[0] = res;
	};

	// Types of nonlinear waves
	enum class WaveType {
		Shock,
		Rarefaction
	};

	// Find star region quantites
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
	};

	// Compute solution configuration
	StarVariables ComputeStarVariables(double roL, double uL, double pL, double roR, double uR, double pR) {
		StarVariables result;

		//	Left state
		params.roL = roL;
		params.uL = uL;
		params.pL = pL;

		//	Right state
		params.roR = roR;
		params.uR = uR;
		params.pR = pR;

		//	Compute star region pressure
		alglib::real_1d_array x;
		alglib::real_1d_array bl;
		alglib::real_1d_array bu;
		x.setlength(1);
		bl.setlength(1);
		bu.setlength(1);

		//Set boundary constraints (non negative pressure)
		bl[0] = 0.0;
		bu[0] = std::numeric_limits<double>::max();

		//Initial guess (possible othe choises, the simplest for now)
		double pStar = 0.5*(params.pL + params.pR);
		x[0] = pStar;

		//Iterative scheme parameters
		double epsg = 1e-10;
		double epsf = 0;
		double epsx = 0;
		double diffstep = 1e-6;
		alglib::ae_int_t maxits = 0;
		alglib::minlmstate state;
		alglib::minlmreport rep;

		alglib::minlmcreatev(1, x, diffstep, state);
		alglib::minlmsetbc(state, bl, bu);
		alglib::minlmsetcond(state, epsg, epsf, epsx, maxits);
		alglib::minlmoptimize(state, fFunction, NULL, (void*)this);
		alglib::minlmresults(state, x, rep);

		//Check if solution converged (TO DO)		
		result.pStar = x[0];

		//Compute star region velocity
		result.uStar = 0.5*(uL + uR) + 0.5*(f2(roR, uR, pR, result.pStar) - f1(roL, uL, pL, result.pStar));

		//Determine nonlinear wave type and properties
		double C = (_gamma - 1.0) / (_gamma + 1.0);
		result.MaxSpeed = 0;

		//Left side of contact
		double pRatioL = result.pStar / pL;
		double aL = sqrt(_gamma * pL / roL);
		if (result.pStar > pL) {
			//Left shock
			result.leftWave = WaveType::Shock;

			//Determine density in star region			
			result.roStarL = pRatioL + C;
			result.roStarL /= C * pRatioL + 1.0;
			result.roStarL *= roL;

			//Determine shock propagation speed			
			result.SL = sqrt((_gamma + 1) * pRatioL / (2 * _gamma) + (_gamma - 1) / (2 * _gamma));
			result.SL = uL - aL * result.SL;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SL));
		}
		else {
			//Left rarefaction
			result.leftWave = WaveType::Rarefaction;

			//Determine density in star region
			result.roStarL = roL * pow(pRatioL, 1.0 / _gamma);

			//Determine rarefaction head propagation speed		
			result.SHL = uL - aL;

			//Determine rarefaction tail propagation speed		
			double aStarL = aL * pow(pRatioL, (_gamma - 1) / (2 * _gamma));
			result.STL = result.uStar - aStarL;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SHL));
		};

		//Right side of contact
		double pRatioR = result.pStar / pR;
		double aR = sqrt(_gamma * pR / roR);
		if (result.pStar > pR) {
			//Right shock
			result.rightWave = WaveType::Shock;

			//Determine density in star region			
			result.roStarR = pRatioR + C;
			result.roStarR /= C * pRatioR + 1.0;
			result.roStarR *= roR;

			//Determine shock propagation speed			
			result.SR = sqrt((_gamma + 1) * pRatioR / (2 * _gamma) + (_gamma - 1) / (2 * _gamma));
			result.SR = uR + aR * result.SR;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SR));
		}
		else {
			//Left rarefaction
			result.rightWave = WaveType::Rarefaction;

			//Determine density in star region
			result.roStarR = roR * pow(pRatioR, 1.0 / _gamma);

			//Determine rarefaction head propagation speed		
			result.SHR = uR + aR;

			//Determine rarefaction tail propagation speed		
			double aStarR = aR * pow(pRatioR, (_gamma - 1) / (2 * _gamma));
			result.STR = result.uStar + aStarR;

			result.MaxSpeed = max(result.MaxSpeed, std::abs(result.SHR));
		};

		return result;
	};
	
public:
	//Constructor
	GodunovSolverPerfectGasEOS(double gamma, double eps, double opPressure) : _gamma(gamma), _operatingPressure(opPressure) {
	};

	//Solve riemann problem
	virtual RiemannProblemSolutionResult ComputeFlux(std::valarray<double> &UL, std::valarray<double> &UR, Vector fn) override {
		RiemannProblemSolutionResult result;
		result.Fluxes.resize(nVariables);

		// Left and right states
		double roL = UL[0];
		double roR = UR[0];
		double pL = GetPressure(UL);
		double pR = GetPressure(UR);

		// velocity vertors
		Vector velocityL = Vector(UL[1], UL[2], UL[3]) / UL[0];
		Vector velocityR = Vector(UR[1], UR[2], UR[3]) / UR[0];
		
		// Compute normal velocity
		double uL = velocityL * fn;
		double uR = velocityR * fn;

		// Star region variables
		StarVariables starValues = ComputeStarVariables(roL, uL, pL, roR, uR, pR);
		
		//Compute flux (Toro p. 219) 
		//Sample exact solution at S = x/t
		double S = 0;
		double ro;
		double u;
		double p;
		double gamma{ _gamma };

		if (starValues.uStar >= S) {
			//Left side of contact
			//Shock wave
			if (starValues.leftWave == WaveType::Shock) {
				if (starValues.SL >= S) {
					// Case a1
					// Left of the shock
					ro = roL;
					u = uL;
					p = pL;
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
					double aL = sqrt(gamma * pL / roL);
					double CL = 2.0 / (gamma + 1) + (gamma - 1) * (uL - S) / ((gamma + 1) * aL);
					double CLRo = pow(CL, 2.0 / (gamma - 1));
					double CLP = pow(CL, 2.0 * gamma / (gamma - 1));

					//Density
					ro = CLRo * roL;

					//Velocity
					u = aL + 0.5 * (gamma - 1) * uL + S;
					u *= 2 / (gamma + 1);

					//Pressure					
					p = CLP * pL;
				};
			};

		} else {
			//Right side of contact

			//Shock wave
			if (starValues.rightWave == WaveType::Shock) {
				if (starValues.SR <= S) {
					//Case a1
					ro = roR;
					u = uR;
					p = pR;
				}
				else {
					//Case a2
					ro = starValues.roStarR;
					u = starValues.uStar;
					p = starValues.pStar;
				};
			};

			//Rarefaction wave
			if (starValues.rightWave == WaveType::Rarefaction) {
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
					double aR = sqrt(gamma * pR / roR);
					double CR = 2.0 / (gamma + 1) - (gamma - 1) * (uR - S) / ((gamma + 1) * aR);
					double CRRo = pow(CR, 2.0 / (gamma - 1));
					double CRP = pow(CR, 2.0 * gamma / (gamma - 1));

					//Density
					ro = CRRo * roR;

					//Velocity
					u = -aR + 0.5 * (gamma - 1) * uR + S;
					u *= 2 / (gamma + 1);

					//Pressure					
					p = CRP * pR;
				};
			};
		};

		// Create conservative variables vector
		Vector Vface = u * fn;
		std::valarray<double> Uface(nVariables);
		Uface[0] = ro;
		Uface[1] = ro * Vface.x;
		Uface[2] = ro * Vface.y;
		Uface[3] = ro * Vface.z;
		Uface[4] = ro * (p / (ro * (gamma - 1.0)) + 0.5 * u * u);

		// write result
		result.Fluxes = F(Uface, fn);
		result.MaxEigenvalue = starValues.MaxSpeed;
		result.Pressure = p;
		result.Velocity = Vface;

		return result;
	};

};

#endif