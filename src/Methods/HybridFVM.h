#ifndef TurboStructured_Methods_HybridFVM
#define TurboStructured_Methods_HybridFVM

#include "Methods\Method.h"
#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "kernel.h"

//Base class for all solution methods that desribe iterations process in detail
class HybridFVM : public Kernel {
public:
	//splitting over directions types
	enum direction {
		X_direction,
		Y_direction,
		Z_direction
	};

	//for splitting over x, y, z directions
	void ExchangeVelocityComponents(double* U, direction dir) {
		if(dir == X_direction) return;
		
		//	(u, v, w) -> (v, w, u)
		if(dir == Y_direction) {
			double temp = U[1];
			U[1] = U[2];
			U[2] = U[3];
			U[3] = temp;
			return;
		};

		//	(u, v, w) -> (w, u, v)
		if(dir == Z_direction) {
			double temp = U[1];
			U[1] = U[3];
			U[3] = U[2];
			U[2] = temp;
			return;
		};
	};

	//Prepare conservative variables for appropriate directions
	//(u v w) components of velocity have feature that u is a normal velocity component everytime and v and w are considered as passive scalar
	//only for x_direction now
	double* PrepareConservativeVariables(int i,int j, int k, direction dir) {
		//EKR need to implement MPI request part and dummy cells
		double* U = getCellValues(i, j, k);
		ExchangeVelocityComponents(U, dir);
		return U;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	void PrepareEigenMatrix(double* &UL, double* &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) {
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

		Matrix E = R*Rinv;
		int test;
	};

	//Compute fluxes in X direction and write values_new array
	void Xfluxes() {
		direction dir = X_direction;

		//initialize data containers
		double* Uim1, *Ui, *Uip1, *Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;	//vectors of characteristic variables computed as Gj = Rinv(i+1/2)*Uj,  j = im1, im, ip1, ip2
		double alf, bet, gam, del;					//coefficients of monotone scheme
		std::vector<double> fr(nVariables), fl(nVariables);					//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(5);
		Matrix R(5, 5);
		Matrix Rinv(5, 5);

		double dx = stepInfo.TimeStep/hx;
		for (int iz = kMin; iz <= kMax; iz++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				//Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = hy;
				if (nDims == 3) fS = hy * hz;

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
						//compute new conservative values
						//std::vector<double> res(nVariables);
						//for(int k = 0; k < nVariables; k++) res[k] = Ui[k] - (fr[k] - fl[k])*dx;
						//ExchangeVelocityComponents(&res.front(), dir);

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
	};

	//Compute fluxes in Y direction and write values_new array
	void Yfluxes() {
		direction dir = Y_direction;

		//initialize data containers
		double* Uim1, *Ui, *Uip1, *Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;	//vectors of characteristic variables computed as Gj = Rinv(i+1/2)*Uj,  j = im1, im, ip1, ip2
		double alf, bet, gam, del;					//coefficients of monotone scheme
		std::vector<double> fr, fl;					//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(5);
		Matrix R(5, 5);
		Matrix Rinv(5, 5);

		double dy = stepInfo.TimeStep/hy;
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
						//compute new conservative values
						std::vector<double> res(nVariables);
						for(int k = 0; k < nVariables; k++) res[k] = Ui[k] - (fr[k]-fl[k])*dy;
						ExchangeVelocityComponents(&res.front(), dir);

						//write them
						int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
						//for (int k = 0; k < nVariables; k++) values_new[sBegin + k] = res[k];
					};
          
					fl = fr;
					Uim1 = Ui;
					Ui = Uip1;
					Uip1 = Uip2;
				 };
			};
		};
	};

	//Compute fluxes in Z direction and write values_new array
	void Zfluxes() {
		direction dir = Z_direction;

		//initialize data containers
		double* Uim1, *Ui, *Uip1, *Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;	//vectors of characteristic variables computed as Gj = Rinv(i+1/2)*Uj,  j = im1, im, ip1, ip2
		double alf, bet, gam, del;					//coefficients of monotone scheme
		std::vector<double> fr, fl;					//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(5);
		Matrix R(5, 5);
		Matrix Rinv(5, 5);

		double dz = stepInfo.TimeStep/hz;
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
						//compute new conservative values
						std::vector<double> res(nVariables);
						for(int k = 0; k < nVariables; k++) res[k] = Ui[k] - (fr[k] - fl[k])*dz;
						ExchangeVelocityComponents(&res.front(), dir);

						//write them
						int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
						//for (int k = 0; k < nVariables; k++) values_new[sBegin + k] = res[k];
					};
          
					fl = fr;
					Uim1 = Ui;
					Ui = Uip1;
					Uip1 = Uip2;
				 };
			};
		};
	};

	void ComputeTimeStep() {
		//variables to collect maximum values
		double dmax  = 0;
		double ccmax = 0;

		//conservative variable temporal container
		std::vector<double> U;
	
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				for (int iz = kMin; iz <= kMax; iz++)
				{
					int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
					U =  std::vector<double>(values.begin() + sBegin, values.begin() + sBegin + nVariables);
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

	
public:
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

	//Explicit time step
	virtual void IterationStep() override {			
		//Determine timestep
		ComputeTimeStep();		

		//Determine order of spacial directions
		Xfluxes();

		//Compute residuals and update solution
		stepInfo.Residual.resize(5, 0);
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				for (int iz = kMin; iz <= kMax; iz++)
				{
					int idx = getSerialIndexLocal(ix, iy, iz);

					//Compute total residual
					stepInfo.Residual[0] += abs(residual[idx * nVariables]);			//ro
					stepInfo.Residual[1] += abs(residual[idx * nVariables + 1]);		//rou
					stepInfo.Residual[2] += abs(residual[idx * nVariables + 2]);		//rov
					stepInfo.Residual[3] += abs(residual[idx * nVariables + 3]);		//row
					stepInfo.Residual[4] += abs(residual[idx * nVariables + 4]);		//roE
				};
			};
		};
		

		//Update conservative variables
		UpdateSolution(stepInfo.TimeStep);
		
		// Update time
		stepInfo.Time += stepInfo.TimeStep;
		stepInfo.Iteration++;
	};	

};

#endif