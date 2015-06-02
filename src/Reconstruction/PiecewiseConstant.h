#ifndef TurboStructured_Reconstruction_PiecewiseConstant
#define TurboStructured_Reconstruction_PiecewiseConstant

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by piecewise constant approximation
class PiecewiseConstant : public IReconstruction {
protected:	
	std::valarray<double> values_; 
public:
	//Piecewise reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		return values_;
	};

	//constructor
	PiecewiseConstant() { };
	PiecewiseConstant(std::valarray<double>& _values) : values_ (_values) { };

	//copy semantics
	//PiecewiseConstant(const PiecewiseConstant& element) : values_ (element.values_) {}; // copy constructor 	
	//PiecewiseConstant& operator=(const PiecewiseConstant& element) { //copy assignment operator
	//	values_ = element.values_;
	//	return *this;
	//};

	//move semantics
	//PiecewiseConstant& operator=(PiecewiseConstant&& element) { //move assignment operator
	//	return *element;
	//};	
};

PiecewiseConstant ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector point) {
	return std::move(PiecewiseConstant(value));
};


#endif