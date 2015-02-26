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
#include "Vector.h"

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