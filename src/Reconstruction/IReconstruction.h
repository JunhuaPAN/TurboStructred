#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility\Vector.h"

//Basic class for all reconstruction classes
class IReconstruction {


public:
	virtual std::valarray<double> SampleSolution(Vector const& point) = 0;
};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector point, int nDim) {
	static_assert(false, "We dont have required function.")
};

#endif