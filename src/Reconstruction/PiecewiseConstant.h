#ifndef TurboStructured_Reconstruction_PiecewiseConstant
#define TurboStructured_Reconstruction_PiecewiseConstant

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>

//Reconstruction by piecewise constant approximation
class PiecewiseConstant : public IReconstruction {
protected:
	const std::valarray<double>& values_;
public:
	//Piecewise reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		return values_;
	};

	//constructor
	PiecewiseConstant(std::valarray<double> const& _values) : values_ (_values) { };
};

PiecewiseConstant ComputeReconstruction(std::vector<std::valarray<double>& > const& values, std::vector<Vector> const& points, std::valarray<double> const& value, Vector const& point) {
	return std::move(PiecewiseConstant(value));
};


#endif