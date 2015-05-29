#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility\Vector.h"

//Basic class for all reconstruction classes
class IReconstruction {


public:
	virtual std::valarray<double> SampleSolution(Vector const& point) = 0;
};

#endif