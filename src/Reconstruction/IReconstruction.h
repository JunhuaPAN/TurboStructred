#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility/Vector.h"

//Basic class for all reconstruction classes
class IReconstruction {
public:
	int nDimensions;
	int nValues;

	virtual std::valarray<double> SampleSolution(Vector const& point) = 0;

	// Default constructor
	IReconstruction() {
		nDimensions = 0;
		nValues = 0;
	};

	//! Return required
	static std::size_t GetBufferLenght(int nD, int nV) { return 0; };
	virtual std::valarray<double> Serialize() = 0;
	virtual void Deserialize(const std::valarray<double>& ) = 0;

	// Update reconstruction center
	virtual void RefrashPosition(Vector point) = 0;
};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	//static_assert(false, "We dont have required function.");
	throw std::runtime_error("We dont have required function.");
};

#endif
