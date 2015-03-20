#ifndef TurboStructured_CUDA_CUDAManager
#define TurboStructured_CUDA_CUDAManager

#include <memory>
#include <iostream>

#include "KernelConfiguration.h"
#include "RiemannSolver.h"
#include "RoeSolverPerfectGasEOS.h"
#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_functions.hpp"

struct CUDAComputeFluxesParameters {
	int nVariables;
	int nDimensions;
	int nFaces;
	double gamma;
	double eps;
	double opPressure;
};

//Incapsulates all CUDA related data and functionality
class CUDAManager {
	static const int TILE_WIDTH = 10;

	//Constant memory for storing parameters
	CUDAComputeFluxesParameters _parameters;

	//Data storage, device side
	RiemannSolver ** _deviceRiemannSolver;
	double *deviceFluxes;
	double *deviceLeftStates;
	double *deviceRightStates;
	Vector *deviceFaceNormals;

	//Device properties
	cudaDeviceProp deviceProp;

	//Kernel launch parameters
	dim3 DimGrid;
	dim3 DimBlock;

	//Kernel and other device functions

	//Init kernel
	void ComputeFluxesInitKernel();

	//Clean kernel
	//__device__ void ComputeFluxesFinalizeKernel() {
	//	//Create RiemannSolver class in GPU memory
	//	delete (*_deviceRiemannSolver);
	//};


	////Compute convective part of flux
	//__device__ void ComputeConvectiveFlux(double *UL, double *UR, Vector* fn, double *Flux) {
	//	Flux[0] = UR[0] + UL[0]; 
	//};

	////Compute total flux
	//__device__ void ComputeFlux(double *UL, double *UR, Vector* fn, double *Flux) {
	//	ComputeConvectiveFlux(UL, UR, fn, Flux);
	//};

	////Computational kernel
	//__device__ void ComputeFluxesSharedKernel(double* LeftStates, double *RightStates, Vector *faceNormals, double *Fluxes, double* Velocities) {

	//	//Shared memory for tiled input
	//	//__shared__ double ULs[TILE_WIDTH];
	//	//__shared__ double URs[TILE_WIDTH];

	//	//Get index and positions of arguments and result buffers
	//	int bx = blockIdx.x;
	//	int tx = threadIdx.x;
	//	int index = bx * blockDim.x + tx;
	//	int nVariables = _parameters.nVariables;
	//	int indexValues = nVariables * index;

	//	//Compute flux
	//	ComputeFlux(&LeftStates[indexValues], &LeftStates[indexValues], &faceNormals[index], &Fluxes[indexValues]);
	//};


public:
	CUDAManager () {
	};

	__host__ void Init(CUDAConfiguration config);

	//void CUDAInitComputeFluxes(CUDAComputeFluxesParameters parameters);

	////Vector operation of computing numerical flux
	//void CUDAComputeFluxes(double* LeftStates, double *RightStates, Vector *faceNormals, double *Fluxes, double* Velocities);

	//void CUDAFinalizeComputeFluxes();
};

#endif
