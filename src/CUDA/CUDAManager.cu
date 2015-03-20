#include "CUDAManager.cu.h"
#include "cuda_runtime.h"
#include "device_functions.hpp"

//__global__ void CUDAManager::ComputeFluxesInitKernel()  {
//	//Create RiemannSolver class in GPU memory
//	(*_deviceRiemannSolver) = new RoeSolverPerfectGasEOS(_parameters.gamma, _parameters.eps, _parameters.opPressure);
//};
//
//__global__ void CUDAManager::CUDAInitComputeFluxes(CUDAComputeFluxesParameters parameters) {
//	//Compute optimal run parameters and hardware capability
//	cudaError_t cudaStatus;
//
//	//Initialize the grid and block dimensions here
//	int maxThreadsPerBlock = deviceProp.maxThreadsPerBlock;
//	int maxSharedMemory = deviceProp.sharedMemPerBlock;
//	int warpSize = deviceProp.warpSize;
//	int nFaces = parameters.nFaces;
//	int nVariables = parameters.nVariables;
//
//	//Determine shared memory threads limit
//	int deviceFaceMemory = 0;
//	deviceFaceMemory = 3 * nVariables * sizeof(double);
//	deviceFaceMemory += sizeof(Vector);
//	int maxThreadsMemory = 1 + (maxSharedMemory - 1) / deviceFaceMemory;
//
//	//Block dimensions
//	int nThreads = std::min(maxThreadsMemory, maxThreadsPerBlock);
//	DimBlock.x = nThreads;
//	DimBlock.y = 1;
//	DimBlock.z = 1;
//
//	//Grid dimensions
//	DimGrid.x = 1 + (nFaces - 1) / nThreads;
//	DimGrid.y = 1;
//	DimGrid.z = 1;
//		
//
//	//Fill device constant memory once
//	cudaStatus = cudaMemcpyToSymbol(&_parameters, &parameters, sizeof(CUDAComputeFluxesParameters));
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaDeviceSynchronize returned error code "<<cudaStatus<<" after launching matrixMultiplyShared!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	//Create RiemannSolver class in GPU memory
//	cudaMalloc(&_deviceRiemannSolver, sizeof(RoeSolverPerfectGasEOS*));
//
//	//Allocate device memory for data storage
//	int size = parameters.nVariables * nFaces;
//
//	cudaStatus = cudaMalloc((void**)&deviceLeftStates, size * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMalloc failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	cudaStatus = cudaMalloc((void**)&deviceRightStates, size * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMalloc failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	cudaStatus = cudaMalloc((void**)&deviceFluxes, size * sizeof(double));
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMalloc failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	cudaStatus = cudaMalloc((void**)&deviceFaceNormals, nFaces * sizeof(Vector));
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMalloc failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	//Make device side preparations
//	ComputeFluxesInitKernel<<<DimGrid, DimBlock>>>();
//};
//
//__global__ void CUDAManager::CUDAFinalizeComputeFluxes() {
//	//Clean device
//	ComputeFluxesFinalizeKernel<<<DimGrid, DimBlock>>>();
//
//	//Free device memory
//	cudaFree(deviceLeftStates);
//	cudaFree(deviceRightStates);
//	cudaFree(deviceFluxes);
//	cudaFree(deviceFaceNormals);
//	cudaFree(_deviceRiemannSolver);
//};
//
//__global__ void CUDAManager::CUDAComputeFluxes(double* LeftStates, double *RightStates, Vector *faceNormals, double *Fluxes, double* Velocities) {
//	cudaError_t cudaStatus;
//	int size = _parameters.nVariables * _parameters.nFaces;
//
//	//Copy left and right states to device memory
//	cudaStatus = cudaMemcpy(deviceLeftStates, LeftStates, size * sizeof(double), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMemcpy failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	cudaStatus = cudaMemcpy(deviceRightStates, RightStates, size * sizeof(double), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMemcpy failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	//Copy face normals to GPU buffer
//	cudaStatus = cudaMemcpy(deviceRightStates, RightStates, size * sizeof(double), cudaMemcpyHostToDevice);
//	if (cudaStatus != cudaSuccess) {
//		std::cerr<<"cudaMemcpy failed!"<<std::endl;
//		CUDAFinalizeComputeFluxes(); exit(0);
//	};
//
//	//Launch calculation kernel
//	ComputeFluxesSharedKernel<<<DimGrid, DimBlock>>>(deviceLeftStates, deviceRightStates, deviceFaceNormals, deviceFluxes);
//
//	//Copy fluxes to host memory
//};


__host__ void CUDAManager::Init(CUDAConfiguration config) {
	int deviceCount;
	cudaError_t cudaStatus = cudaGetDeviceCount(&deviceCount);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaGetDeviceCount failed!");
		return;
	}	

	for (int dev = 0; dev < deviceCount; dev++) {
		cudaGetDeviceProperties(&deviceProp, dev);

		if (dev == 0) {
			if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
				std::cout<<"No CUDA GPU has been detected"<<std::endl;
				return;
			} else if (deviceCount == 1) {
				std::cout<<"There is 1 device supporting CUDA"<<std::endl;
			} else {
				std::cout<<"There are "<<deviceCount<<" devices supporting CUDA"<<std::endl;
			}
		}

		std::cout<<"Device "<<dev<<" name: "<<deviceProp.name<<std::endl;
		std::cout<<" Computational Capabilities: "<< deviceProp.major<< "."<< deviceProp.minor<<std::endl;
		std::cout<<" Maximum global memory size: "<< deviceProp.totalGlobalMem<<std::endl;
		std::cout<<" Maximum constant memory size: "<< deviceProp.totalConstMem<<std::endl;
		std::cout<<" Maximum shared memory size per block: "<< deviceProp.sharedMemPerBlock<<std::endl;
		std::cout<<" Maximum block dimensions: "<< deviceProp.maxThreadsDim[0]<< " x "<<
													deviceProp.maxThreadsDim[1]<< " x "<<
													deviceProp.maxThreadsDim[2]<<std::endl;
		std::cout<<" Maximum grid dimensions: "<< deviceProp.maxGridSize[0]<< " x "<<
													deviceProp.maxGridSize[1]<< " x "<<
													deviceProp.maxGridSize[2]<<std::endl;
		std::cout<<" Warp size: "<<deviceProp.warpSize<<std::endl;
	};

	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);	
	if (cudaStatus != cudaSuccess) {
		std::cerr<<"cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";
		exit(0);
	};

	// Get device properties
	cudaGetDeviceProperties(&deviceProp, 0);
};