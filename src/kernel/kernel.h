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

	//Function for shifting indexes
	inline int shiftIndex(int ind[3], int shift[3]) {
		return 0;
	};

	//Grid information	
	int nDims; //Number of dimensions
	int nX; //Number of cells in x dimension
	int nY; //Number of cells in y dimension
	int nZ; //Number of cells in z dimension
	double hx, hy, hz;	//cell sizes
	bool IsPeriodicX; //X periodicity
	bool IsPeriodicY; //Y periodicity
	bool IsPeriodicZ; //Z periodicity
	std::vector<double> CoordinateX; //Cell center coordinates
	std::vector<double> CoordinateY; //Cell center coordinates
	std::vector<double> CoordinateZ; //Cell center coordinates	
	
	//Solution information
	//Gas model information
	int nVariables; // number of conservative variables
	double gamma; //
	std::vector<double> values; //array for storing values
	std::vector<double> values_new;	//array of new values

	//Get serial index for cell
	inline int getSerialIndex(int i, int j, int k) {
		int sI = (i * nX * nY + j * nX + k);
		return sI;
	};
	
	//Get cell values
	inline std::vector<double> getCellValues(int i, int j, int k) {
		int sBegin = getSerialIndex(i, j, k) * nVariables;
		return std::vector<double>(values.begin() + sBegin, values.begin() + sBegin + nVariables);
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
		for (int i = iMin; i < iMax; i++) {
			for (int j = jMin; j < jMax; j++) {
				for (int k = kMin; k < kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					int sBegin = getSerialIndex(i, j, k) * nVariables;
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
		int nlocalX = nX / dimsCart[0];
		int nlocalY = nY / dimsCart[1];
		int nlocalZ = nZ / dimsCart[2];		
		iMin = rankCart[0] * nlocalX;
		iMax = (rankCart[0]+1) * nlocalX;
		jMin = rankCart[1] * nlocalY;
		jMax = (rankCart[1]+1) * nlocalY;
		kMin = rankCart[2] * nlocalZ;
		kMax = (rankCart[2]+1) * nlocalZ;

		//Initialize gas model parameters
		gamma = 1.4;

		//Initialize calculation parameters
		CFL = 0.5;
		MaxTime = 0.2;
		MaxIteration = 100000;
		SaveSolutionSnapshotTime = 0;	
		SaveSolutionSnapshotIterations = 1;

		//Initialize boundary conditions
		//TO DO
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
		for (int i = iMin; i < iMax; i++) {
			for (int j = jMin; j < jMax; j++) {
				for (int k = kMin; k < kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					std::vector<double> U = getCellValues(i,j,k);
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
	void ExchangeVelocityComponents(std::vector<double> &U, direction dir) {
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
	std::vector<double> PrepareConservativeVariables(int i,int j, int k, direction dir) {
		//EKR need to implement MPI request part and dummy cells
		std::vector<double> U = getCellValues(i, j, k);
		ExchangeVelocityComponents(U, dir);
		return U;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) {
		return;
	};

	//Compute fluxes in X direction and write values_new array
	void Xfluxes(double dt) {
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

		double dx = dt/hx;
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
						for(int k = 0; k < nVariables; k++) res[k] = Ui[k] - (fr[k]-fl[k])*dx;
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

	//Explicit time step
	void ExplicitTimeStep() {		
		//Determine timestep
		double dt = 0.01;

		//Update solution		
		for (int i = iMin; i < iMax; i++) {
			for (int j = jMin; j < jMax; j++) {
				for (int k = kMin; k < kMax; k++) {
					//Obtain cell data
					double x = CoordinateX[i];
					double y = CoordinateY[j];
					double z = CoordinateZ[k];
					std::vector<double> U = getCellValues(i,j,k);
					double ro = U[0];
					double u = U[1] / ro;
					double v = U[2] / ro;
					double w = U[3] / ro;
					double e = U[4] - ro * (u*u+v*v+w*w) / 2.0;
					double P = (gamma - 1.0) * ro * e;


				};
			};
		};
	};	

	//Run calculation
	void Run() {
		//Calculate snapshot times order of magnitude
		int snapshotTimePrecision = 0;
		if (SaveSolutionSnapshotTime > 0) {
			snapshotTimePrecision = 1 - std::floor(std::log10(SaveSolutionSnapshotTime));
		};

		//Start timer
		double workTime = 0.0;
		clock_t start, stop;
		/* Start timer */
		assert((start = clock())!=-1);

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

		/* Stop timer */
		stop = clock();
		workTime = (double) (stop-start)/CLOCKS_PER_SEC;		
		std::cout<<"Total work time = " << workTime << " seconds.\n";			
	};

	//Finalize kernel
	void Finilaze() {
		//Finalize MPI
	};

};

#endif