#ifndef TurboStructured_kernel_ParallelManager
#define TurboStructured_kernel_ParallelManager

#include "mpi.h"
#include <chrono>

struct StructuredGridInfo {
	int nDims;
	int nX;
	int nY;
	int nZ;
	bool IsPeriodicX;
	bool IsPeriodicY;
	bool IsPeriodicZ;
};

class ParallelManager {
	std::chrono::high_resolution_clock::duration _idleDuration; //total idle time
	bool isInitilized;	//was MPI initialization successful
	int _nProcessors;	//total number of processes
	int _rank;			//current processor rank in _comm
	MPI_Comm _comm;		//global communicator

public:	
	//Cartezian topology
	int rankCart[3];	//Position of process in cartesian mpi grid
	int dimsCart[3];	//Dimensions of cartesian mpi grid
	bool periodic[3];	//Periodic boundaries

	//Accessor functions
	inline MPI_Comm getComm() {
		return _comm;
	};

	inline int getRank() {
		return _rank;
	};

	inline int getProcessorNumber() {
		return _nProcessors;
	};

	inline std::chrono::high_resolution_clock::duration getIdleTime() {
		return _idleDuration;
	};

	//Constructor
	ParallelManager(int *argc, char **argv[]) {
		Init(argc, argv);
		isInitilized = false;
	};

	//Initialize MPI programm
	void Init(int *argc, char **argv[]) {
		MPI_Init(argc, argv);		
		_comm = MPI_COMM_WORLD;
		MPI_Comm_size(_comm, &_nProcessors);
		MPI_Comm_rank(_comm, &_rank);
		isInitilized = true;
		_idleDuration = std::chrono::high_resolution_clock::duration(0);
	};

	void InitCartesianTopology(StructuredGridInfo info) {
		//Simple one dimensional topology for now
		rankCart[0] = _rank;
		rankCart[1] = 0;
		rankCart[2] = 0;
		dimsCart[0] = _nProcessors;
		dimsCart[1] = 1;
		dimsCart[2] = 1;
		periodic[0] = info.IsPeriodicX;
		periodic[1] = info.IsPeriodicY;
		periodic[2] = info.IsPeriodicZ;
	};

	//Finilize MPI programm
	void Finalize() {
		MPI_Finalize();
		isInitilized = false;
	};

	//Is master node
	inline bool IsMaster() {
		return _rank == 0;
	};

	//Barrier syncronization
	void Barrier() {
		auto before = std::chrono::high_resolution_clock::now();
		MPI_Barrier(_comm);
		auto after = std::chrono::high_resolution_clock::now();
		_idleDuration += after - before; // update idle time duration
	};
};

#endif