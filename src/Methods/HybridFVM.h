#ifndef TurboStructured_Methods_HybridFVM
#define TurboStructured_Methods_HybridFVM

#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "kernel.h"

//Base class for Hybrid methods (works with one phase ideal gas model)
class HybridFVM : public Kernel {
public:
	//splitting over directions types
	enum direction {
		X_direction,
		Y_direction,
		Z_direction
	};

	//for splitting over x, y, z directions
	inline void ExchangeVelocityComponents(std::vector<double> &U, direction dir) {
		if(dir == X_direction) return;
		
		//	(u, v, w) -> (v, u, w)
		if(dir == Y_direction) {
			double temp = U[1];
			U[1] = U[2];
			U[2] = temp;
			return;
		};

		//	(u, v, w) -> (w, v, u)
		if(dir == Z_direction) {
			double temp = U[1];
			U[1] = U[3];
			U[3] = temp;
			return;
		};
	};

	//Prepare conservative variables for appropriate directions
	//(u v w) components of velocity have feature that u is a normal velocity component everytime and v and w are considered as passive scalar
	//only for x_direction now
	inline std::vector<double> PrepareConservativeVariables(int i, int j, int k, direction dir) {
		std::vector<double> U(nVariables, 0);
		double* Value = getCellValues(i, j, k);
		for(int nv = 0; nv < nVariables; nv++) U[nv] = Value[nv];
		ExchangeVelocityComponents(U, dir);
		return U;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	virtual void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) {
		//Grunaisen
		double gr = gamma - 1.0;

		//Left cell
		double rol = UL[0];
		double ul  = UL[1]/rol;
		double vl = UL[2]/rol;
		double wl  = UL[3]/rol;
		double El  = UL[4]/rol;
		double el = El - 0.5*(ul*ul + vl*vl + wl*wl);	//specific internal energy
		
		//EOS of ideal gas
		double pl  = gr*rol*el;	
		double hl = El + pl/rol;	//total entalpy

		//Right cell
		double ror = UR[0];
		double ur  = UR[1]/ror;
		double vr = UR[2]/ror;
		double wr  = UR[3]/ror;
		double Er  = UR[4]/ror;
		double er = Er - 0.5*(ur*ur + vr*vr + wr*wr);	//specific internal energy

		//EOS of ideal gas
		double pr  = gr*ror*er;	
		double hr = Er + pr/ror;	//total entalpy

		//avarage values on faces
		double roa = sqrt(rol*ror);
		double sl  = sqrt(rol)/(sqrt(rol) + sqrt(ror));			//density based weight coefficients
		double sr  = sqrt(ror)/(sqrt(rol) + sqrt(ror));
		double ua = ul*sl+ur*sr;
		double va = vl*sl+vr*sr;
		double wa = wl*sl+wr*sr;
		double ha  = hl*sl+hr*sr;
		double q2  = ua*ua + va*va + wa*wa;
		double q2p = q2/2;							//kinetic energy
		double c2 = gr*(ha-q2p);
		double c  = sqrt(c2);						//averaged sound speed

		//write eigenvalues
		eigenvals[0] = ua - c;
		eigenvals[1] = ua;
		eigenvals[2] = ua;
		eigenvals[3] = ua;
		eigenvals[4] = ua + c;

		//write right eigenvectors
		R.element[0][0] = 1;
		R.element[1][0] = ua - c;
		R.element[2][0] = va;
		R.element[3][0] = wa;
		R.element[4][0] = ha - c*ua;

		R.element[0][1] = 1;
		R.element[1][1] = ua;
		R.element[2][1] = va;
		R.element[3][1] = wa;
		R.element[4][1] = ha - c2/gr;
		
		R.element[0][2] = 0;
		R.element[1][2] = 0;
		R.element[2][2] = 1;
		R.element[3][2] = 0;
		R.element[4][2] = va;

	    R.element[0][3] = 0;
		R.element[1][3] = 0;
		R.element[2][3] = 0;
		R.element[3][3] = 1;
		R.element[4][3] = wa;
		
		R.element[0][4] = 1;
		R.element[1][4] = ua + c;
		R.element[2][4] = va;
		R.element[3][4] = wa;
		R.element[4][4] = ha + c*ua;
		
		double dc2 = 2*c2;

		//write inverse matrix (left eigenvectors are rows)
		Rinv.element[0][0] = (gr*q2p + ua*c)/dc2;
		Rinv.element[0][1] = (-1.0)*(c + ua*gr)/dc2;
		Rinv.element[0][2] = -va*gr/dc2;
		Rinv.element[0][3] = -wa*gr/dc2;
		Rinv.element[0][4] = gr/dc2;

		Rinv.element[1][0] = 1.0 - gr*q2p/c2;
		Rinv.element[1][1] = ua*gr/c2;
		Rinv.element[1][2] = va*gr/c2;
		Rinv.element[1][3] = wa*gr/c2;
		Rinv.element[1][4] = -gr/c2;

		Rinv.element[2][0] = -va;
		Rinv.element[2][1] = 0;
		Rinv.element[2][2] = 1;
		Rinv.element[2][3] = 0;
		Rinv.element[2][4] = 0;

		Rinv.element[3][0] = -wa;
		Rinv.element[3][1] = 0;
		Rinv.element[3][2] = 0;
		Rinv.element[3][3] = 1;
		Rinv.element[3][4] = 0;

		Rinv.element[4][0] = (gr*q2p - ua*c)/dc2;
		Rinv.element[4][1] = (c - ua*gr)/dc2;
		Rinv.element[4][2] = va*gr/dc2;
		Rinv.element[4][3] = wa*gr/dc2;
		Rinv.element[4][4] = gr/dc2;
	};
	 
	//Compute fluxes in X direction and write values_new array
	void Xfluxes(double dt) {
		direction dir = X_direction;
		for (double& r : residual) r = 0;	//Nullify residuals

		//initialize data containers
		std::vector<double> Uim1, Ui, Uip1, Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;	//vectors of characteristic variables computed as Gjj = Rinv(i+1/2)*Ujj,  j = im1, im, ip1, ip2
		double alf, bet, gam, del;					//coefficients of monotone scheme
		std::vector<double> fr(nVariables), fl(nVariables);					//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(nVariables);
		Matrix R(nVariables, nVariables);
		Matrix Rinv(nVariables, nVariables);

		//Compute face square
		double fS = 0;
		if (nDims == 1) fS = 1.0;
		if (nDims == 2) fS = hy;
		if (nDims == 3) fS = hy * hz;

		double dx = dt/hx;			//need for local courant number
		for (int iz = kMin; iz <= kMax; iz++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				Uim1 = PrepareConservativeVariables(iMin-2, iy, iz, dir);
				Ui = PrepareConservativeVariables(iMin-1, iy, iz, dir);
				Uip1 = PrepareConservativeVariables(iMin, iy, iz, dir);

				for (int ix = iMin; ix <= iMax + 1; ix++)
				{
					Uip2 = PrepareConservativeVariables(ix+1, iy, iz, dir);
					PrepareEigenMatrix(Ui, Uip1, R, Rinv, eigenvals);

					//Characteristic values
					Gim1 = Rinv*Uim1;
					Gi = Rinv*Ui;
					Gip1 = Rinv*Uip1;
					Gip2 = Rinv*Uip2;

					//Compute characteristic flux
					for (int k = 0; k < nVariables; k++)
					{
						//Check hybridization condition
						double dG = Gip1[k] - Gi[k];
						double dGlr = (Gip2[k] - Gip1[k]) - (Gi[k] - Gim1[k]);
						double a = eigenvals[k];
						double Gk = a*dG*dGlr;		//the condition for stencil switching
						double c = a*dx;			//local courant number
						
						//Compute stencil coefficients
						if(a >= 0)
						{
							if(Gk >= 0)
							{
								alf = -0.5*(1-c);
								bet = 1-alf;
								gam = 0;
								del = 0;
							} else {
								alf = 0;
								gam = 0.5*(1-c);
								bet = 1-gam;
								del = 0;
							};
						} else {
							if(Gk >= 0)
							{
								alf = 0;
								del = -0.5*(1-c);
								gam = 1-del;
								bet = 0;
							} else {
								alf = 0;
								bet = 0.5*(1-c);
								gam = 1-bet;
								del = 0;
							};
						};
            
						//Compute characteristic flux
						fr[k] = Gim1[k]*alf + Gi[k]*bet + Gip1[k]*gam + Gip2[k]*del;
						fr[k] *= a;
					};

					//Return to the conservative variables
					fr = R*fr;
					
					//for all inner cells
					if(ix > iMin)
					{
						//Fluxes difference equals residual 
						int idx = getSerialIndexLocal(ix - 1, iy, iz);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] = - (fr[nv] - fl[nv]) * fS;
					};
          
					fl = fr;
					Uim1 = Ui;
					Ui = Uip1;
					Uip1 = Uip2;
				 };
			};
		};

		RewriteDuringSolution(dt);
	};

	//Compute fluxes in Y direction and write values_new array
	void Yfluxes(double dt) {
		direction dir = Y_direction;
		for (double& r : residual) r = 0;	//Nullify residuals

		//initialize data containers
		std::vector<double> Uim1, Ui, Uip1, Uip2;						//vectors of conservative variables in j-1  j  j+1  j+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;						//vectors of characteristic variables computed as Gjj = Rinv(j+1/2)*Ujj,  jj = jm1, jm, jp1, jp2
		double alf, bet, gam, del;										//coefficients of monotone scheme
		std::vector<double> fr(nVariables), fl(nVariables);				//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(nVariables);
		Matrix R(nVariables, nVariables);
		Matrix Rinv(nVariables, nVariables);

		//Compute face square
		double fS = 0;
		if (nDims == 2) fS = hx;
		if (nDims == 3) fS = hx * hz;

		double dy = dt/hy;			//need for local courant number
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iz = kMin; iz <= kMax; iz++)
			{
				Uim1 = PrepareConservativeVariables(ix, jMin-2, iz, dir);
				Ui = PrepareConservativeVariables(ix, jMin-1, iz, dir);
				Uip1 = PrepareConservativeVariables(ix, jMin, iz, dir);

				for (int iy = jMin; iy <= jMax + 1; iy++)
				{
					Uip2 = PrepareConservativeVariables(ix, iy+1, iz, dir);
					PrepareEigenMatrix(Ui, Uip1, R, Rinv, eigenvals);

					//Characteristic values
					Gim1 = Rinv*Uim1;
					Gi = Rinv*Ui;
					Gip1 = Rinv*Uip1;
					Gip2 = Rinv*Uip2;

					//Compute characteristic flux
					for (int k = 0; k < nVariables; k++)
					{
						//Check hybridization condition
						double dG = Gip1[k] - Gi[k];
						double dGlr = (Gip2[k] - Gip1[k]) - (Gi[k] - Gim1[k]);
						double a = eigenvals[k];
						double Gk = a*dG*dGlr;		//the condition for stencil switching
						double c = a*dy;			//local courant number
						
						//Compute stencil coefficients
						if(a >= 0)
						{
							if(Gk >= 0)
							{
								alf = -0.5*(1-c);
								bet = 1-alf;
								gam = 0;
								del = 0;
							} else {
								alf = 0;
								gam = 0.5*(1-c);
								bet = 1-gam;
								del = 0;
							};
						} else {
							if(Gk >= 0)
							{
								alf = 0;
								del = -0.5*(1-c);
								gam = 1-del;
								bet = 0;
							} else {
								alf = 0;
								bet = 0.5*(1-c);
								gam = 1-bet;
								del = 0;
							};
						};
            
						//Compute characteristic flux
						fr[k] = Gim1[k]*alf + Gi[k]*bet + Gip1[k]*gam + Gip2[k]*del;
						fr[k] *= a;
					};

					//Return to the conservative variables
					fr = R*fr;
					
					//for all inner cells
					if(iy > jMin)
					{
						//compute delta of fluxes
						int idx = getSerialIndexLocal(ix, iy - 1, iz);
						std::vector<double> deltaFlux(nVariables);
						for(int nv = 0; nv < nVariables; nv++) deltaFlux[nv] = (fr[nv] - fl[nv]) * fS;
						ExchangeVelocityComponents(deltaFlux, dir);

						//Fluxes difference equals residual 
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] = - deltaFlux[nv];
					};
          
					fl = fr;
					Uim1 = Ui;
					Ui = Uip1;
					Uip1 = Uip2;
				 };
			};
		};
		RewriteDuringSolution(dt);
	};

	//Compute fluxes in Z direction and write values_new array			START FROM HERE
	void Zfluxes(double dt) {
		direction dir = Z_direction;
		for (double& r : residual) r = 0;	//Nullify residuals

		//initialize data containers
		std::vector<double> Uim1, Ui, Uip1, Uip2;											//vectors of conservative variables in k-1  k  k+1  k+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;									//vectors of characteristic variables computed as Gjj = Rinv(k+1/2)*Ujj,  jj = km1, km, kp1, kp2
		double alf, bet, gam, del;													//coefficients of monotone scheme
		std::vector<double> fr(nVariables), fl(nVariables);							//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(nVariables);
		Matrix R(nVariables, nVariables);
		Matrix Rinv(nVariables, nVariables);

		//Compute face square
		double fS = hx * hy;

		double dz = dt/hz;			//need for local courant number
		for (int iy = jMin; iy <= jMax; iy++)
		{
			for (int ix = iMin; ix <= iMax; ix++)
			{
				Uim1 = PrepareConservativeVariables(ix, iy, kMin-2, dir);
				Ui = PrepareConservativeVariables(ix, iy, kMin-1, dir);
				Uip1 = PrepareConservativeVariables(ix, iy, kMin, dir);

				for (int iz = kMin; iz <= kMax + 1; iz++)
				{
					Uip2 = PrepareConservativeVariables(ix, iy, iz+1, dir);
					PrepareEigenMatrix(Ui, Uip1, R, Rinv, eigenvals);

					//Characteristic values
					Gim1 = Rinv*Uim1;
					Gi = Rinv*Ui;
					Gip1 = Rinv*Uip1;
					Gip2 = Rinv*Uip2;

					//Compute characteristic flux
					for (int k = 0; k < nVariables; k++)
					{
						//Check hybridization condition
						double dG = Gip1[k] - Gi[k];
						double dGlr = (Gip2[k] - Gip1[k]) - (Gi[k] - Gim1[k]);
						double a = eigenvals[k];
						double Gk = a*dG*dGlr;		//the condition for stencil switching
						double c = a*dz;			//local courant number
						
						//Compute stencil coefficients
						if(a >= 0)
						{
							if(Gk >= 0)
							{
								alf = -0.5*(1-c);
								bet = 1-alf;
								gam = 0;
								del = 0;
							} else {
								alf = 0;
								gam = 0.5*(1-c);
								bet = 1-gam;
								del = 0;
							};
						} else {
							if(Gk >= 0)
							{
								alf = 0;
								del = -0.5*(1-c);
								gam = 1-del;
								bet = 0;
							} else {
								alf = 0;
								bet = 0.5*(1-c);
								gam = 1-bet;
								del = 0;
							};
						};
            
						//Compute characteristic flux
						fr[k] = Gim1[k]*alf + Gi[k]*bet + Gip1[k]*gam + Gip2[k]*del;
						fr[k] *= a;
					};

					//Return to the conservative variables
					fr = R*fr;
					
					//for all inner cells
					if(iz > kMin)
					{
						////compute delta of fluxes
						int idx = getSerialIndexLocal(ix, iy, iz - 1);
						std::vector<double> deltaFlux(nVariables);
						for(int nv = 0; nv < nVariables; nv++) deltaFlux[nv] = (fr[nv] - fl[nv]) * fS;
						ExchangeVelocityComponents(deltaFlux, Z_direction);

						//Fluxes difference equals residual 
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] = deltaFlux[nv];
					};
          
					fl = fr;
					Uim1 = Ui;
					Ui = Uip1;
					Uip1 = Uip2;
				 };
			};
		};
		RewriteDuringSolution(dt);
	};

	//compute full tine step
	virtual void ComputeTimeStep() {
		//variables to collect maximum values
		double dmax  = 0;
		double ccmax = 0;	//maximum of speed of sound (not used yet)

		//conservative variable temporal container		
	
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				for (int iz = kMin; iz <= kMax; iz++)
				{
					int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
					std::valarray<double>&& U = values[std::slice(sBegin, nVariables, 1)];
					double ro = U[0];
					double rou = U[1];
					double rov = U[2];
					double row = U[3];
					double roE = U[4];
					double uu = rou/ro;
					double vv = rov/ro;
					double ww = row/ro;
					double E = roE/ro;
					double ek = (uu*uu + vv*vv + ww*ww)/2;
					double e = E - ek;
					
					//Compute sound speed for ideal gas
					double gr = gamma - 1.0;
					double c = sqrt(gr*gamma*e);

					double um = fabs(uu) + c;
					ccmax = std::max(ccmax, c);
					um /= hx;
					if(um >= dmax) dmax = um;
          
					if (nDims > 1) {
						um = fabs(vv) + c;
						um /= hy;
						if(um >= dmax) dmax = um;
					};

					if (nDims > 2) {
						um = fabs(ww) + c;
						um /= hz;
						if(um >= dmax) dmax = um; 
					};
				};
			};
		};

		double dt = CFL/dmax;
		if(stepInfo.NextSnapshotTime>0) dt = std::min(stepInfo.NextSnapshotTime - stepInfo.Time, dt);
		stepInfo.TimeStep = dt;
	};

	//Update solution without reseting of residuals (for intermidiate step)
	void RewriteDuringSolution(double dt) {
		double volume = hx * hy * hz;
		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);

					//Update cell values
					for(int nv = 0; nv < nVariables; nv++) {
						values[idx * nVariables + nv] += residual[idx * nVariables + nv] * dt / volume;
						//Compute total residual
						stepInfo.Residual[nv] += abs(residual[idx * nVariables + nv]);
					};
				};
			};
		};
	};

	//Constuctor inherited
	HybridFVM(int* argc, char **argv[]) : Kernel(argc, argv) {};

	//Method parameters
	double CFL;

	//Initizalization
	virtual void Init(KernelConfiguration& config) override {
		Kernel::Init(config);
		CFL = config.methodConfiguration.CFL;
	};

	//Number of required dummy layers
	virtual int GetDummyCellLayerSize() override {
		return 2;
	};

	//one iteration for spatial splitting scheme accuracy by time (Update solution and write residuals)
	void CompleteOneIteration() {
		stepInfo.Residual.resize(nVariables);
		//Determine order of spacial directions
		double dt = 0;
		switch (nDims) {
		case 1 :
			dt = stepInfo.TimeStep;
			Xfluxes(dt);
			// source terms
			if(isExternalForce == true) {
				for (double r : residual) r = 0;
				ProcessExternalForces();
				RewriteDuringSolution(dt);
			};
			if(isExternalAccelaration == true) {
				for (double r : residual) r = 0;
				ProcessExternalAcceleration();
				RewriteDuringSolution(dt);
			};
			break;
		case 2 :
			dt = 0.5*stepInfo.TimeStep;
			if ( stepInfo.Iteration % 2 == 0 ) {
				//half of a Step
				Xfluxes(dt);
				ComputeDummyCellValues();
				Yfluxes(dt);
				ComputeDummyCellValues();

				//second half
				Yfluxes(dt);
				ComputeDummyCellValues();
				Xfluxes(dt);
				
				// source terms
				if(isExternalForce == true) {
					for (double r : residual) r = 0;
					ProcessExternalForces();
					RewriteDuringSolution(2*dt);
				};
				if(isExternalAccelaration == true) {
					for (double r : residual) r = 0;
					ProcessExternalAcceleration();
					RewriteDuringSolution(2*dt);
				};
			} else {
				//half of Step
				Yfluxes(dt);
				ComputeDummyCellValues();
				Xfluxes(dt);
				ComputeDummyCellValues();

				//second half
				Xfluxes(dt);
				ComputeDummyCellValues();
				Yfluxes(dt);

				// source terms
				if(isExternalForce == true) {
					for (double r : residual) r = 0;
					ProcessExternalForces();
					RewriteDuringSolution(2*dt);
				};
				if(isExternalAccelaration = true) {
					for (double r : residual) r = 0;
					ProcessExternalAcceleration();
					RewriteDuringSolution(2*dt);
				};
			};
			break;
		};
	};

	//Explicit time step
	virtual void IterationStep() override {
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Iteration started\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Determine timestep
		ComputeTimeStep();

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Timestep calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		CompleteOneIteration();

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Residual calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Update conservative variables
		UpdateResiduals();
		
		// Update time
		stepInfo.Time += stepInfo.TimeStep;
		stepInfo.Iteration++;

		pManager->Barrier();
	};	

};

#endif