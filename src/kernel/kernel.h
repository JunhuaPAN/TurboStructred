#ifndef TurboStructured_kernel_kernel
#define TurboStructured_kernel_kernel

#include <array>
#include <vector>
#include <chrono>
#include <thread>
#include <functional>
#include <sstream>
#include <cassert>
#include <fstream>

#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"

//Step info
class StepInfo {
public:
	double Time;
	double TimeStep;	
	int Iteration;
	std::vector<double> Residual;
	bool NeedResidual;
};

//Calculation kernel
class Kernel {
public:
	//Parallel run information
	int rank; //Rank of process
	int rankCart[3]; //Position of process in cartesian mpi grid
	int dimsCart[3]; //Dimensions of cartesian mpi grid
	int iMin; 
	int iMax;
	int jMin;
	int jMax;
	int kMin;
	int kMax;

	//splitting over directions types
	enum direction {
		X_direction,
		Y_direction,
		Z_direction
	};

	//Current step information
	StepInfo stepInfo;

	//Grid information	
	int nDims; //Number of dimensions
	int nX; //Number of cells in x dimension
	int nY; //Number of cells in y dimension
	int nZ; //Number of cells in z dimension
	int nXAll; //Number of cells in x dimension including dummy layers
	int nYAll; //Number of cells in y dimension including dummy layers
	int nZAll; //Number of cells in z dimension including dummy layers
	//Same for local cells
	int nlocalX;
	int nlocalY;
	int nlocalZ;
	int nlocalXAll;
	int nlocalYAll;
	int nlocalZAll;
	double hx, hy, hz;	//cell sizes
	bool IsPeriodicX; //X periodicity
	bool IsPeriodicY; //Y periodicity
	bool IsPeriodicZ; //Z periodicity
	std::vector<double> CoordinateX; //Cell center coordinates
	std::vector<double> CoordinateY; //Cell center coordinates
	std::vector<double> CoordinateZ; //Cell center coordinates
	int dummyCellLayers; //number of dummy cell layers
	
	//Solution information
	//Gas model information
	int nVariables; // number of conservative variables
	double gamma; //
	std::vector<double> values; //array for storing values
	std::vector<double> values_new;	//array of new values

	//Get serial index for cell
	inline int getSerialIndexGlobal(int i, int j, int k) {
		int sI = (k * nXAll * nYAll + j * nXAll + i);
		return sI;
	};

	inline int getSerialIndexLocal(int i, int j, int k) {
		int sI =  (k - kMin + dummyCellLayers) * nlocalXAll * nlocalYAll + (j - jMin + dummyCellLayers) * nlocalXAll + (i - iMin + dummyCellLayers);
		return sI;
	};
	
	//Get cell values
	inline double* getCellValues(int i, int j, int k) {
		int sBegin = getSerialIndexLocal(i, j, k) * nVariables;		
		return &values[sBegin];
	};

	//Calculation parameters
	int RungeKuttaOrder;
	double CFL;	
	double MaxIteration;
	double MaxTime;
	double CurrentTime;
	double NextSnapshotTime;
	double SaveSolutionSnapshotTime;
	int SaveSolutionSnapshotIterations;

	//Constructor
	Kernel() {
	};

	//Set initial conditions
	void SetInitialConditions(std::function<std::vector<double>(Vector r)> initF) {		
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					int sBegin = getSerialIndexLocal(i, j, k) * nVariables;
					std::vector<double> U = initF(Vector(x,y,z));
					for (int varn = 0; varn < nVariables; varn++) values[sBegin + varn] = U[varn];					
				};
			};
		};
	};

	//Initialize kernel
	void Init(KernelConfiguration& config) {
		//Initialize MPI
		rank = 0;
		rankCart[0] = 0;
		rankCart[1] = 0;
		rankCart[2] = 0;
		dimsCart[0] = 1;
		dimsCart[1] = 1;
		dimsCart[2] = 1;		

		//Initialize local grid
		nDims = config.nDims;
		nX = config.nX;
		IsPeriodicX = config.isPeriodicX;
		nY = 1;
		IsPeriodicY = true;
		nZ = 1;
		IsPeriodicZ = true;
		if (nDims > 1) {
			nY = config.nY;
			IsPeriodicY = config.isPeriodicY;
		};
		if (nDims > 2) {
			nZ = config.nZ;
			IsPeriodicZ = config.isPeriodicZ;
		};
		dummyCellLayers = 2; //number of dummy cell layers
		nlocalX = nX / dimsCart[0];
		nlocalY = nY / dimsCart[1];
		nlocalZ = nZ / dimsCart[2];		
		iMin = rankCart[0] * nlocalX + dummyCellLayers;
		iMax = (rankCart[0]+1) * nlocalX + dummyCellLayers - 1;
		jMin = rankCart[1] * nlocalY + dummyCellLayers;
		jMax = (rankCart[1]+1) * nlocalY + dummyCellLayers - 1;
		kMin = rankCart[2] * nlocalZ + dummyCellLayers;
		kMax = (rankCart[2]+1) * nlocalZ + dummyCellLayers - 1;
		double Lx = config.LX;
		double Ly = config.LY;
		double Lz = config.LZ;		

		//Generate cells (uniform grid)
		nXAll = nX + 2 * dummyCellLayers;
		nYAll = nY + 2 * dummyCellLayers;
		nZAll = nZ + 2 * dummyCellLayers;
		nlocalXAll = nlocalX + 2 * dummyCellLayers;
		nlocalYAll = nlocalY + 2 * dummyCellLayers;
		nlocalZAll = nlocalZ + 2 * dummyCellLayers;
		CoordinateX.resize(nXAll);
		CoordinateY.resize(nYAll);
		CoordinateZ.resize(nZAll);
		double dx = Lx / nlocalX;
		double xMin = 0;
		xMin = xMin - (dummyCellLayers * dx) + 0.5 * dx;
		double xMax = Lx;
		xMax = xMax + (dummyCellLayers * dx) - 0.5 * dx;
		for (int i = iMin - dummyCellLayers; i <= iMax + dummyCellLayers; i++) {
			double x = xMin + (xMax - xMin) * 1.0 * i / nlocalX;							
			CoordinateX[i] = x;
		};		

		double dy = Ly / nlocalY;
		double yMin = 0;
		yMin = yMin - (dummyCellLayers * dy) + 0.5 * dy;
		double yMax = Ly;
		yMax = yMax + (dummyCellLayers * dy) - 0.5 * dy;
		for (int j = jMin - dummyCellLayers; j <= jMax + dummyCellLayers; j++) {
			double y = yMin + (yMax - yMin) * 1.0 * j / nlocalY;					
			CoordinateY[j] = y;
		};

		double dz = Lz / nlocalZ;
		double zMin = 0;
		zMin = zMin - (dummyCellLayers * dz) + 0.5 * dz;
		double zMax = Lz;
		zMax = zMax + (dummyCellLayers * dz) - 0.5 * dz;
		for (int k = kMin - dummyCellLayers; k <= kMax + dummyCellLayers; k++) {
			double z = zMin + (zMax - zMin) * 1.0 * k / nlocalZ;
			CoordinateZ[k] = z;
		};

		//Initialize gas model parameters
		gamma = 1.4;
		nVariables = 5;

		//Allocate memory for values
		values.resize(nVariables * nlocalXAll * nlocalYAll * nlocalZAll);		

		//Initialize calculation parameters
		CFL = 0.5;
		MaxTime = 0.2;
		MaxIteration = 100000;
		SaveSolutionSnapshotTime = 0;	
		SaveSolutionSnapshotIterations = 1;

		//Initialize boundary conditions
		//TO DO
	};

	//Exchange values between processors
	void ExchangeValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;
		
		//X direction exchange
		for (j = jMin; j <= jMax; j++) {
			for (k = kMin; k <= kMax; k++) {
				for (int layer = 1; layer <= dummyCellLayers; layer++) {
					// minus direction
					i = iMin - layer; // layer index					
					if ((rankCart[0] == 0) && (IsPeriodicX)) {
						int iMirror = nX - i + 1; // mirror cell index
						int sI = getSerialIndexLocal(i, j, k);
						int sIMirror = getSerialIndexLocal(iMirror, j, k);
						for (int nv = 0; nv < nVariables; nv++) {
							values[sI * nVariables + nv] = values[sIMirror * nVariables + nv]; //TO DO serial version for now	
						};
					};

					// plus direction		
					i = iMax + layer; // layer index					
					if ((rankCart[0] == dimsCart[0] - 1) && (IsPeriodicX)) {
						int iMirror =  dummyCellLayers - layer; // mirror cell index
						int sI = getSerialIndexLocal(i, j, k);
						int sIMirror = getSerialIndexLocal(iMirror, j, k);
						for (int nv = 0; nv < nVariables; nv++) {
							values[sI * nVariables + nv] = values[sIMirror * nVariables + nv]; //TO DO serial version for now	
						};
					};
				};
			};
		};
	};

	//Compute dummy cell values as result of boundary conditions and interprocessor exchange communication
	void ComputeDummyCellValues() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//Interprocessor exchange
		ExchangeValues();

		//X direction		
		for (j = jMin; j <= jMax; j++) {
			for (k = kMin; k <= kMax; k++) {
				for (int layer = 1; layer <= dummyCellLayers; layer++) {
					i = iMin - layer; // layer index
					if (!IsPeriodicX) {
						//Apply boundary conditions
					};
				};
			};
		};
	};

	//Save solution
	void SaveSolution(std::string fname) {
		std::ofstream ofs(fname);
		//Header
		ofs.open("history.dat", std::ofstream::out);
		ofs<<"VARIABLES = ";
		ofs<<"\""<<"X"<<"\" ";
		ofs<<"\""<<"ro"<<"\" ";
		ofs<<"\""<<"u"<<"\" ";
		ofs<<"\""<<"P"<<"\" ";
		ofs<<"\""<<"e"<<"\" ";
		ofs<<std::endl;

		//Solution
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					double* U = getCellValues(i,j,k);
					double ro = U[0];
					double u = U[1] / ro;
					double v = U[2] / ro;
					double w = U[3] / ro;
					double e = U[4] - ro * (u*u+v*v+w*w) / 2.0;
					double P = (gamma - 1.0) * ro * e;

					//Write to file
					ofs<<x<<" ";
					ofs<<ro<<" ";
					ofs<<u<<" ";
					ofs<<P<<" ";
					ofs<<e<<" ";
					ofs<<std::endl;
				};
			};
		};

		ofs.close();
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
	void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) {
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
		R.element[0][1] = ua - c;
		R.element[0][2] = va;
		R.element[0][3] = wa;
		R.element[0][4] = ha - c*ua;

		R.element[1][0] = 1;
		R.element[1][1] = ua;
		R.element[1][2] = va;
		R.element[1][3] = wa;
		R.element[1][4] = ha - c2/gr;
		
		R.element[2][0] = 0;
		R.element[2][1] = 0;
		R.element[2][2] = 1;
		R.element[2][3] = 0;
		R.element[2][4] = va;

	    R.element[3][0] = 0;
		R.element[3][1] = 0;
		R.element[3][2] = 0;
		R.element[3][3] = 1;
		R.element[3][4] = wa;
		
		R.element[4][0] = 1;
		R.element[4][1] = ua + c;
		R.element[4][2] = va;
		R.element[4][3] = wa;
		R.element[4][4] = ha + c*ua;
		
		double dc2 = 2*c2;

		//write left eigenvectors
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
	void Xfluxes() {
		direction dir = X_direction;

		//initialize data containers
		std::vector<double> Uim1, Ui, Uip1, Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
		std::vector<double> Gim1, Gi, Gip1, Gip2;	//vectors of characteristic variables computed as Gj = Rinv(i+1/2)*Uj,  j = im1, im, ip1, ip2
		double alf, bet, gam, del;					//coefficients of monotone scheme
		std::vector<double> fr, fl;					//left and right fluxes
		
		//Initialize eigenvalues and eigenvectors structures
		std::vector<double> eigenvals(5);
		Matrix R(5, 5);
		Matrix Rinv(5, 5);

		double dx = stepInfo.TimeStep/hx;
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
						//compute new conservative values
						std::vector<double> res(nVariables);
						for(int k = 0; k < nVariables; k++) res[k] = Ui[k] - (fr[k] - fl[k])*dx;
						ExchangeVelocityComponents(res, dir);

						//write them
						int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
						for (int k = 0; k < nVariables; k++) values_new[sBegin + k] = res[k];
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
		std::vector<double> Uim1, Ui, Uip1, Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
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
						ExchangeVelocityComponents(res, dir);

						//write them
						int sBegin = getSerialIndex(ix, iy, iz) * nVariables;
						for (int k = 0; k < nVariables; k++) values_new[sBegin + k] = res[k];
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
		std::vector<double> Uim1, Ui, Uip1, Uip2;	//vectors of conservative variables in i-1  i  i+1  i+2 cells
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
						ExchangeVelocityComponents(res, dir);

						//write them
						int sBegin = getSerialIndex(ix, iy, iz) * nVariables;
						for (int k = 0; k < nVariables; k++) values_new[sBegin + k] = res[k];
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
					int sBegin = getSerialIndex(ix, iy, iz) * nVariables;
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
          
					um = fabs(vv) + c;
					um /= hy;
					if(um >= dmax) dmax = um;

					um = fabs(ww) + c;
					um /= hz;
					if(um >= dmax) dmax = um; 
				};
			};
		};

		stepInfo.TimeStep = CFL/dmax;
	};

	//a half of full time step
	void HalfStep() {
		//Compute residuals if needed
		if(stepInfo.NeedResidual == true)
		{
			std::vector<double> res(values.size());
			for(int i = 0; i < res.size(); i++) {
				res[i] = values_new[i] - values[i];
			};

			stepInfo.Residual.resize(5, 0);
			for (int ix = iMin; ix <= iMax; ix++)
			{
				for (int iy = jMin; iy <= jMax; iy++)
				{
					for (int iz = kMin; iz <= kMax; iz++)
					{
						int idx = getSerialIndex(ix, iy, iz);
						stepInfo.Residual[0] = abs(res[idx]);			//ro
						stepInfo.Residual[1] = abs(res[idx + 1]);		//rou
						stepInfo.Residual[2] = abs(res[idx + 2]);		//rov
						stepInfo.Residual[3] = abs(res[idx + 3]);		//row
						stepInfo.Residual[4] = abs(res[idx + 4]);		//roE
					};
				};
			};
		};
	};

	//Explicit time step
	void ExplicitTimeStep() {		
		//Determine timestep
		ComputeTimeStep();

		//Determine order of spacial directions
		Xfluxes();

		//Update conservative variables
		values = values_new;

		//refresh dummy valuies
	};	

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = 1 - std::floor(std::log10(SaveSolutionSnapshotTime));
		};

		//Start timer		

		//Calc loop		
		if (rank == 0) {
			std::cout<<"Calculation started!\n";
		};

		for (stepInfo.Iteration = 0; stepInfo.Iteration <= MaxIteration; stepInfo.Iteration++) {
			//Calculate one time step
			ExplicitTimeStep();

			//Output step information						
			if (rank == 0) {
				std::cout<<"Iteration = "<<stepInfo.Iteration<<"; Total time = "<< stepInfo.Time << "; Time step = " <<stepInfo.TimeStep << "; RMSrou = "<<stepInfo.Residual[1]<<"\n";			
			};			

			//Solution snapshots
			//Every few iterations
			if ((SaveSolutionSnapshotIterations != 0) && (stepInfo.Iteration % SaveSolutionSnapshotIterations) == 0) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<"dataI"<<stepInfo.Iteration<<".cgns";						
				SaveSolution(snapshotFileName.str());
			};

			//Every fixed time interval
			if ((SaveSolutionSnapshotTime > 0) && (NextSnapshotTime == stepInfo.Time)) {
				//Save snapshot
				std::stringstream snapshotFileName;
				snapshotFileName.str(std::string());
				snapshotFileName<<std::fixed;
				snapshotFileName.precision(snapshotTimePrecision);								
				snapshotFileName<<"dataT"<<stepInfo.Time<<".cgns";							
				SaveSolution(snapshotFileName.str());

				//Adjust next snapshot time
				NextSnapshotTime += SaveSolutionSnapshotTime;
			};

			//Save history			

			//Convergence criteria
			if (stepInfo.Iteration == MaxIteration) {
				if (rank == 0) {
					std::cout<<"Maximal number of iterations reached.\n";
				};
				break;
			};

			if (stepInfo.Time >= MaxTime) {
				if (rank == 0) {
					std::cout<<"Maximal time reached.\n";
				};				
				break;
			};

			//Synchronize			
		};

		//Synchronize
		if (rank == 0) {
			std::cout<<"Calculation finished!\n";
		};

		//Stop timer
		double workTime = 0.0;
		std::cout<<"Total work time = " << workTime << " seconds.\n";			
	};

	//Finalize kernel
	void Finilaze() {
		//Finalize MPI
	};

};

#endif