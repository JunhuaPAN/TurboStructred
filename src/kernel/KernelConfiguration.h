#ifndef TurboStructured_kernel_KernelConfiguration
#define TurboStructured_kernel_KernelConfiguration

//Class that manages all configurable parameters
class KernelConfiguration {
public:
	//Grid sizes and configuration
	int nDims; //Number of dimensions
	int nX; //Number of cells in x dimension
	int nY; //Number of cells in y dimension
	int nZ; //Number of cells in z dimension	
	bool isPeriodicX; //X periodicity
	bool isPeriodicY; //Y periodicity
	bool isPeriodicZ; //Z periodicity
	double LX; //X size
	double LY; //Y size
	double LZ; //Z size

	//Solution method parameters

	//Calculation parameters
};

#endif