#ifndef TurboStructured_Methods_Method
#define TurboStructured_Methods_Method

#include <vector>
#include <memory>
#include "KernelConfiguration.h"


////Base class for all solution methods that desribe iterations process in detail
//class Method {	
//    //virtual ~Method() = default; // enable deletion of a Derived* through a Base* c++11
//public:
//	//Solution information
//	std::vector<double> values;		//array for storing values
//	std::vector<double> residual;	//array of residuals ( - d(U * Volume)  / dt)
//	StepInfo* stepInfo;
//	virtual ~Method() {};
//	virtual int GetDummyCellLayerSize() = 0; // request as much dummy cell layers as required
//	virtual void Init( KernelConfiguration& config) = 0;
//	virtual void IterationStep( StepInfo& _stepInfo ) = 0; // main function
//};

#endif