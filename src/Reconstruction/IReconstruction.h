#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility/Vector.h"
#include "utility/VariablesTransition.h"

// Enumerate all tupes of reconstruction here
enum class Reconstruction {
	PiecewiseConstant,
	Linear2PointsStencil,
	ENO2PointsStencil,
};

// Basic class for all reconstruction classes
class IReconstruction {
public:
	// base class perameters
	int nValues;
	int nDims;
	std::valarray<double> vals;		// values at zero point (0, 0, 0)

	// Sample solution or get reconstruction at point corresponding to the zero point
	virtual inline std::valarray<double> SampleSolution(Vector const& point) = 0;
	
	// Default constructor
	IReconstruction() {};
	IReconstruction(std::valarray<double> _vals, int _nDims, int _nValues) : vals(_vals), nDims(_nDims), nValues(_nValues) { };

	// Aditional inner initialization for derived classes
	virtual void Init(int _nValues, int _nDims) {
		nValues = _nValues;
		nDims = _nDims;
	};

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return 0;
	};
	virtual std::valarray<double> Serialize() = 0;
	virtual void Deserialize(const std::valarray<double>& _values) = 0;
};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	throw std::runtime_error("We dont have required function.");
};

#endif
