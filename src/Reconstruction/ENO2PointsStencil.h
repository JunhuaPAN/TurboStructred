#ifndef TurboStructured_Reconstruction_ENO2PointsStencil
#define TurboStructured_Reconstruction_ENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by ENO interpolation procedure for 2 points stencil
class ENO2PointsStencil : public IReconstruction {
protected:
	std::valarray<double> values_;
	Vector center_;
	std::valarray< Vector > gradients_;

public:
	//ENO second order reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		int&& size = gradients_.size();
		std::valarray<double> res(size);

		Vector dr = (point - center_);
		for (int i = 0; i < size; i++) res[i] = dr * gradients_[i];
		return std::move(values_ + res);
	};

	//constructor
	ENO2PointsStencil() { };
	ENO2PointsStencil(std::valarray<double>& _values, Vector& _center, std::valarray< Vector >& _gradients, int _nDimensions, int _nValues) :
		values_(_values),
		center_(_center),
		gradients_(_gradients)
		{
			nDimensions = _nDimensions;
			nValues = _nValues;
		};

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV + 3 + nV * 3;		// values.size + center.size + gradients.size
	};
	virtual std::valarray<double> Serialize() override {
		std::vector<double> res;
		for (auto& r : values_) res.push_back(r);
		for (auto& r : static_cast<std::valarray<double>> (center_)) res.push_back(r);
		for (auto& r : gradients_) {
			res.push_back(r.x);
			res.push_back(r.y);
			res.push_back(r.z);
		};
		return std::valarray<double> {res.data(), res.size()};
	};
	virtual void Deserialize(const std::valarray<double>& _values) {
		values_ = {&_values[0], (size_t)nValues};
		center_ = Vector(_values[nValues], _values[nValues + 1], _values[nValues + 2]);

		gradients_.resize(nValues);
		for (int i = 0; i < gradients_.size(); i++) {
			gradients_[i] = Vector(_values[nValues + 3 * (i + 1)], _values[nValues + 3 * (i + 1) + 1], _values[nValues + 3 * (i + 1) + 2]);
		};

		return;
	};
};

template<>
ENO2PointsStencil ComputeReconstruction<ENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	//compute Gradients
	int size = value.size();
	std::valarray< Vector > gradients(size);
	
	//1D case
	if (nDim == 1) {
		//std::valarray<double> grad_l(size), grad_r(size);
		double grad_l, grad_r;
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[0][i]) / (point.x - points[0].x);
			grad_r = (values[1][i] - value[i]) / (points[1].x - point.x);
			gradients[i] = Vector(std::min(std::abs(grad_l), std::abs(grad_r)), 0, 0);
		};

	};

	//2D case

	//3D case
	
	//ENO2PointsStencil res = ENO2PointsStencil(value, point, gradients);
	return std::move(ENO2PointsStencil(value, point, gradients, nDim, size));
};


#endif