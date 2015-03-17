#ifndef TurboStructured_RiemannSolvers_RoeSolverPerfectGasEOS
#define TurboStructured_RiemannSolvers_RoeSolverPerfectGasEOS

#include "RiemannSolver.h"

class RoeSolverPerfectGasEOS : public RiemannSolver {
	//Number of variables
	static const int nVariables = 5;

	//Required info
	double _gamma; //Specific heat ratio (inherited from gas model)
	double _eps;	//Harten entropy correction coefficient (optional, 0 by default)
	double _operatingPressure; //Operating pressure (optional, 0 by default)

	//Numerical flux
	std::vector<double> F(double* U, Vector n)
	{		
		std::vector<double> res(nVariables,0);		
		double ro = U[0];
		double vx = U[1]/ro;
		double vy = U[2]/ro;
		double vz = U[3]/ro;
		double roE = U[4];	//ro*e
		double p = (_gamma-1.0)*(roE-ro*(vx*vx+vy*vy+vz*vz)/2.0) - _operatingPressure;		
		double vn = vx*n.x + vy*n.y + vz*n.z;

		res[0] = ro*vn;
		res[1] = ro*vn*vx + n.x*p;
		res[2] = ro*vn*vy + n.y*p;
		res[3] = ro*vn*vz + n.z*p;
		res[4] = vn*(roE + p);
		return res;
	};

	//fabs() and Harten's entropy correction procedure
	double Harten(double z, double eps) 
	{
		z = fabs(z);
		if (z<eps) z = ((z*z)/eps + eps)*0.5;
		return z;
	};

public:	
	//Constructor
	RoeSolverPerfectGasEOS(double gamma, double eps, double opPressure) : _gamma(gamma), _eps(eps), _operatingPressure(opPressure) {
	};

	//Get required width of dummy cell layer
	virtual int GetDummyCellLayerSize() override {
		return 1;
	};

	//Solve riemann problem
	virtual RiemannProblemSolutionResult ComputeFlux(std::vector<double *> values, Vector fn) override {
		RiemannProblemSolutionResult result;
		result.Fluxes.resize(5);

		// Obtain values
		double *UL = values[0];
		double *UR = values[1];

		// Calculate symmetric flux part
		std::vector<double> FL = F(UL, fn);
		std::vector<double> FR = F(UR, fn);
		for (int i = 0; i < nVariables; i++) result.Fluxes[i] = FL[i] + FR[i];

		// Calculates stabilization term which is a part of numerical
		// flux vector i.e. |A|(Q{R}-Q{L})
		// Roe type averaging procedure first
		double ro_l = UL[0];
		double ro_r = UR[0];
		Vector velocity_l, velocity_r;
		velocity_l.x = UL[1]/ro_l;
		velocity_l.y = UL[2]/ro_l;
		velocity_l.z = UL[3]/ro_l;
		velocity_r.x = UR[1]/ro_r;
		velocity_r.y = UR[2]/ro_r;
		velocity_r.z = UR[3]/ro_r;
		double e_l = UL[4]/ro_l;
		double e_r = UR[4]/ro_r;
		double k;
		k = 0.5*(velocity_l*velocity_l);		//kinetik energy
		double h_l = (e_l-k)*_gamma + k;	//enthalpy
		k = 0.5*(velocity_r*velocity_r);		//kinetik energy
		double h_r = (e_r-k)*_gamma + k;	//enthalpy
		double ro = sqrt(ro_l*ro_r);          // (Roe averaged) density		
		double ql  = sqrt(ro_l)/(sqrt(ro_l)+sqrt(ro_r));
		double qr  = sqrt(ro_r)/(sqrt(ro_l)+sqrt(ro_r));
		Vector velocity = ql*velocity_l + qr*velocity_r;	// (Roe averaged) velocity	
		double h  = ql*h_l + qr*h_r;  // (Roe averaged) total enthalpy
		//Proceed to solution
		double phi2 = 0.5*(_gamma - 1)*(velocity*velocity);
		double dn = 1.0;
		double c = sqrt((_gamma - 1)*h - phi2);	//acoustic velocity
		assert(c>0);
		//Debug	
		double uw = velocity * fn;
		double eig_max = fabs(uw)+c*dn;	
		double AA1 = Harten(uw, _eps*eig_max);       // AA1, AA3, AA1 -
		double AA3 = Harten(uw+c*dn, _eps*eig_max);  // eigenvalues of a flux vector
		double AA4 = Harten(uw-c*dn, _eps*eig_max);  // Jacobian matrix
		double Eig1= AA1;
		double Eiga=(AA3 - AA4)*0.5/(dn*c);
		double Eigb=(AA3 + AA4)*0.5 - AA1;
		//parametrs vectors Qa and Qb (i guess)
		std::vector<double> Qa(5, 0);
		std::vector<double> Qb(5, 0);
		Qa[0]=0;
		Qa[1]=fn.x;
		Qa[2]=fn.y;
		Qa[3]=fn.z;
		Qa[4]=uw;    
		Qb[0]=1;
		Qb[1]=velocity.x;
		Qb[2]=velocity.y;
		Qb[3]=velocity.z;
		Qb[4]=h;
		//Calculate solution
		//Some quotients
		double R1 =phi2*ro_r-(_gamma-1)*ro_r*(velocity*velocity_r-e_r);	//PR =R1
		double D1 =ro_r*fn*(velocity_r - velocity);				//DR =D1

		double R2 =phi2*ro_l-(_gamma-1)*ro_l*(velocity*velocity_l-e_l);	//PR =R1
		double D2 =ro_l*fn*(velocity_l - velocity);		

		double C2 = 1.0/(c*c);
		double DN2= 1.0/(dn*dn);

		double *ul = UL;
		double *ur = UR;
		for(int i=0; i<5; i++){
				  result.Fluxes[i] +=(Eig1*ul[i] + Eiga*(R2*Qa[i]     +D2*Qb[i])
									+ Eigb*(R2*Qb[i]*C2  +D2*Qa[i]*DN2)
						-Eig1*ur[i] - Eiga*(R1*Qa[i]     +D1*Qb[i])
									- Eigb*(R1*Qb[i]*C2  +D1*Qa[i]*DN2));
				  result.Fluxes[i] *= 0.5;
		};	

		//Estimate for maximum eigenvalue
		result.MaxEigenvalue = eig_max;	

		//Estimate for velocity
		result.Velocity = velocity;

		//Estimate for pressure
		result.Pressure = _operatingPressure;

		return result;
	};

};

#endif