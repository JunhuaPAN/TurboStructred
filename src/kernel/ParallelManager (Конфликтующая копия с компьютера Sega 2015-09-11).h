#ifndef TurboStructured_kernel_ParallelManager
#define TurboStructured_kernel_ParallelManager

#include "mpi.h"
#include <chrono>

enum class Direction {
	XDirection = 0,
	YDirection = 1,
	ZDirection = 2
};

enum class SubDirection {
	Left = -1,
	Right = 1
};

struct StructuredGridInfo {
	int nDims;
	int nX;
	int nY;
	int nZ;
	bool IsPeriodicX;
	bool IsPeriodicY;
	bool IsPeriodicZ;
	
	//compression coefficient of grid sizes towards the borders in X Y and Z directions (as geometric progression)
	double qx;
	double qy;
	double qz;
};

class ParallelManager {
	std::chrono::high_resolution_clock::duration _idleDuration; //total idle time
	bool isInitilized;	//was MPI initialization successful
	int _nProcessors;	//total number of processes
	int _rank;			//current processor rank in _comm
	MPI_Comm _comm;		//global communicator

public:	
	//Cartezian topology
	MPI_Comm _commCart; //cartesian communicator
	int _nDims;			//number of dimensions
	int _rankCart;		//current processor rank in _commCart
	int _nProcessorsCart;	//total number of processes
	int rankCart[3];	//Position of process in cartesian mpi grid
	int dimsCart[3];	//Dimensions of cartesian mpi grid
	bool periodic[3];	//Periodic boundaries

	//Get rank using cartesian index
	int GetRankByCartesianIndex(int i, int j, int k) {
		int rank = -1;
		int coords[3];
		coords[0] = i;
		coords[1] = j;
		coords[2] = k;
		int MPIResult = MPI_Cart_rank(_commCart, coords, &rank);
		if (MPIResult == MPI_PROC_NULL) rank = -1;
		return rank;
	};

	//Get rank using cartesian index
	int GetRankByCartesianIndexShift(int i, int j, int k) {
		int rank = -1;
		int coords[3];
		coords[0] = rankCart[0] + i;
		coords[1] = rankCart[1] + j;
		coords[2] = rankCart[2] + k;
		if (!periodic[0] && ((coords[0] < 0) || (coords[0] >= dimsCart[0]))) return rank;
		if (!periodic[1] && ((coords[1] < 0) || (coords[1] >= dimsCart[1]))) return rank;
		if (!periodic[2] && ((coords[2] < 0) || (coords[2] >= dimsCart[2]))) return rank;
		coords[0] = (coords[0] + dimsCart[0]) % dimsCart[0];
		coords[1] = (coords[1] + dimsCart[1]) % dimsCart[1];
		coords[2] = (coords[2] + dimsCart[2]) % dimsCart[2];
		int MPIResult = MPI_Cart_rank(_commCart, coords, &rank);
		if (MPIResult == MPI_PROC_NULL) rank = -1;
		return rank;
	};

	void InitCartesianTopology(StructuredGridInfo info) {
		//Simple one dimensional topology for now
		_nDims = info.nDims;				
		dimsCart[0] = 0;
		dimsCart[1] = 1;
		if (_nDims > 1) dimsCart[1] = 0;		
		dimsCart[2] = 1;
		if (_nDims > 2) dimsCart[2] = 0;

		rankCart[0] = 0;
		rankCart[1] = 0;
		rankCart[2] = 0;
		//if (_nDims == 2) dimsCart[1] = _nProcessors;
		//if (_nDims == 3) dimsCart[2] = _nProcessors;

		periodic[0] = info.IsPeriodicX;
		periodic[1] = info.IsPeriodicY;
		periodic[2] = info.IsPeriodicZ;	

		//Determine grid dimensions (could be optimized)
		MPI_Dims_create(_nProcessors, _nDims, dimsCart);

		//Create cartesian topology communicator
		int period[3];
		period[0] = (int)periodic[0];
		period[1] = (int)periodic[1];
		period[2] = (int)periodic[2];
		int MPIResult = MPI_Cart_create(_comm, _nDims, dimsCart, period, true, &_commCart);

		//Get rank for
		MPI_Comm_rank(_commCart, &_rankCart);

		//Get number of processors
		MPI_Comm_size(_commCart, &_nProcessorsCart);

		//Determine catezian coordinates
		MPI_Cart_coords(_commCart, _rankCart, _nDims, rankCart);

		//Sync
		Barrier();
	};

	//Wrappers for MPI functions

	//Sum of integers
	inline int Sum(int& x) {
		int res;
		MPI_Allreduce(&x, &res, 1, MPI_INT, MPI_SUM, _comm);
		return res;
	};

	//Sum of doubles
	inline double Sum(double& x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_SUM, _comm);
		return res;
	};

	//Minimum
	inline double Min(double& x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_MIN, _comm);
		return res;
	};

	//Maximum
	inline double Max(double& x) {
		double res;
		MPI_Allreduce(&x, &res, 1, MPI_LONG_DOUBLE, MPI_MAX, _comm);
		return res;
	};


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

	template <class PeriodType> 
	inline double getIdleTime() {
		//using durationType = std::chrono::duration<double, PeriodType>;
		//return std::chrono::duration_cast<durationType>(_idleDuration).count();
		return std::chrono::duration_cast<std::chrono::duration<double, PeriodType>>(_idleDuration).count();
	};

	//Is master node
	inline bool IsMaster() {
		return _rank == 0;
	};
	
	//Is first node
	inline bool IsFirstNode() {
		return _rank == 0;
	};

	//Is last node
	inline bool IsLastNode() {
		return _rank == _nProcessors - 1;
	};

	//Is master node in cartesian topology
	inline bool IsMasterCart() {
		return _rankCart == 0;
	};

	//Is first node in cartesian topology in given direction
	inline bool IsFirstNodeCart(Direction direction) {
		int dimIndex = (int)direction;
		bool result = (!periodic[dimIndex]) && (rankCart[dimIndex] == 0);
		return result;
	};

	//Is last node in cartesian topology in given direction
	inline bool IsLastNodeCart(Direction direction) {
		int dimIndex = (int)direction;
		bool result = (!periodic[dimIndex]) && (rankCart[dimIndex] == dimsCart[dimIndex] - 1);
		return result;
	};

	//Constructor
	ParallelManager(int *argc, char **argv[]) {
		Init(argc, argv);
		isInitilized = false;		//TO DO WHY false??
	};

	//Initialize MPI programm
	void Init(int *argc, char **argv[]) {
		MPI_Init(argc, argv);		
		_comm = MPI_COMM_WORLD;
		MPI_Comm_size(_comm, &_nProcessors);
		MPI_Comm_rank(_comm, &_rank);
		isInitilized = true;
		_idleDuration = std::chrono::high_resolution_clock::duration(0);

		//std::cout<<"rank = "<<_rank<<", Process started.\n";
		//Sync
		Barrier();
	};	

	//Finilize MPI programm
	void Finalize() {
		std::cout<<"rank = "<<_rank<<", Process finished\n";
		Barrier();
		MPI_Finalize();
		isInitilized = false;
	};	

	//Barrier syncronization
	void Barrier() {
		auto before = std::chrono::high_resolution_clock::now();
		MPI_Barrier(_comm);
		auto after = std::chrono::high_resolution_clock::now();
		_idleDuration += after - before; // update idle time duration
	};

	//Blocking wait
	void Wait(int rank) {
		MPI_Status status;
		auto before = std::chrono::high_resolution_clock::now();
		MPI_Recv(nullptr, 0, MPI_INT, rank, 0, _comm, &status);
		auto after = std::chrono::high_resolution_clock::now();
		_idleDuration += after - before; // update idle time duration
	};

	//Signal for waiting process
	void Signal(int rank) {
		MPI_Send(nullptr, 0, MPI_INT, rank, 0, _comm);
	};

	//Exchange routines
	void SendRecvDouble(MPI_Comm& comm, int destRank, int sourceRank, double *sendbuf, int nSend, double *recvbuf, int nRecv) {
		MPI_Status status;
		int tag = 0;
		int rank = -1;		
		MPI_Sendrecv(sendbuf, nSend, MPI_LONG_DOUBLE, destRank, tag, recvbuf, nRecv, MPI_LONG_DOUBLE, sourceRank, tag, comm, &status); 
	};

	//Gather one integer from each process
	void GatherCounts(int N, std::vector<int>& result) {		
		if (IsMaster()) {			
			result.resize(_nProcessors);
			MPI_Gather(&N, 1, MPI_INT, &result[0], 1, MPI_INT, 0, _comm);
		} else {
			MPI_Gather(&N, 1, MPI_INT, NULL, 1, MPI_INT, 0, _comm);
		};		
	};
};

#endif