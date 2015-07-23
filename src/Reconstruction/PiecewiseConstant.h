#ifndef TurboStructured_Reconstruction_PiecewiseConstant
#define TurboStructured_Reconstruction_PiecewiseConstant

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by piecewise constant approximation
class PiecewiseConstant : public IReconstruction {
protected:
	std::valarray<double> values_{0};
public:
	//Piecewise reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		return values_;
	};

	//constructor
	PiecewiseConstant() { };
	PiecewiseConstant(std::valarray<double>& _values, int _nDimensions) : values_ (_values)
	{
		nDimensions = _nDimensions;
		nValues = _values.size();
	};

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) { return nV; };
	virtual std::valarray<double> Serialize() override {
		return values_;
	};
	virtual void Deserialize(const std::valarray<double>& _values) {
		values_ = _values;
	};

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

template<>
PiecewiseConstant ComputeReconstruction<PiecewiseConstant>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector point, int nDims) {
	return std::move(PiecewiseConstant(value, nDims));
};


#endif