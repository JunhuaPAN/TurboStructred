#ifndef TurboStructured_Reconstruction_PiecewiseConstant
#define TurboStructured_Reconstruction_PiecewiseConstant

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include "Reconstruction/IReconstruction.h"


//Reconstruction by piecewise constant approximation
class PiecewiseConstant : public IReconstruction {
public:

	//Piecewise reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) override { return vals; }

	// Constructors
	PiecewiseConstant() { };
	PiecewiseConstant(std::valarray<double>& _values, int _nValues, int _nDim) : IReconstruction(_values, _nDim, _nValues) {};

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV;
	};
	virtual std::valarray<double> Serialize() override {
		return vals;
	};
	virtual void Deserialize(const std::valarray<double>& _values) override {
		vals = _values;
	};
};


template<>
PiecewiseConstant ComputeReconstruction<PiecewiseConstant>(std::vector<std::valarray<double> >& values, std::vector<CellInfo>& cells, std::valarray<double>& value, CellInfo& cell, int nDim) {
	return std::move(PiecewiseConstant(value, value.size(), nDim));
};

#endif
