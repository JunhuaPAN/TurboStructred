#ifndef TurboStructured_Reconstruction_PiecewiseConstant
#define TurboStructured_Reconstruction_PiecewiseConstant

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include "Reconstruction/IReconstruction.h"


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
	PiecewiseConstant(std::valarray<double>& _values) : values_ (_values){};

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) { return nV; };
	virtual std::valarray<double> Serialize() override {
		return values_;
	};
	virtual void Deserialize(const std::valarray<double>& _values) override {
		values_ = _values;
	};

	virtual void RefreshPosition(Vector point) override {};
};


template<>
PiecewiseConstant ComputeReconstruction<PiecewiseConstant>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	return std::move(PiecewiseConstant(value));
};

#endif
