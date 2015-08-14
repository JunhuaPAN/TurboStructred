#ifndef TurboStructured_Reconstruction_WENO2PointsStencil
#define TurboStructured_Reconstruction_WENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by convex combination of 2 points stencil reconstructions
class WENO2PointsStencil : public IReconstruction {
protected:
	std::valarray<double> values_;
	Vector center_;
	std::valarray< Vector > gradientsL_;	// left stencil derivatives   TO DO vector< valarray< double> > > gradients_
	std::valarray< Vector > gradientsR_;	// right ones
	std::vector< std::valarray < double > > weights_;	//weights

public:
	//WENO 3 order reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		//result array
		std::valarray<double> res(values_.size());

		// compute WENO reconstruction
		for (int i = 0; i < values_.size(); i++) {

			//consider two linear reconstructions
			double val0 = values_[i] + gradientsR_[i] * (point - center_);
			double val1 = values_[i] + gradientsL_[i] * (point - center_);

			//consider left and right cases
			if (point.x > center_.x) res[i] = weights_[0][i] * val0 + weights_[1][i] * val1;
			else res[i] = weights_[2][i] * val0 + weights_[3][i] * val1;
		};

		return std::move(res);
	};

	//constructor
	WENO2PointsStencil() {	};
	WENO2PointsStencil(std::valarray<double>& _values, Vector& _center, std::valarray< Vector >& _gradientsL, std::valarray< Vector >& _gradientsR,
		std::vector< std::valarray <double> >& _weights) :
		values_(_values),
		center_(_center),
		gradientsL_(_gradientsL),
		gradientsR_(_gradientsR),
		weights_(_weights) { };

	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV + nD + 2 * nV * nD + 4 * nD;		// values.size + center.size + gradients.size + weights.size
	};
	virtual std::valarray<double> Serialize() override {
		std::vector<double> res;
		for (auto& r : values_) res.push_back(r);
		for (auto& r : static_cast<std::valarray<double>> (center_)) res.push_back(r);
		for (auto& r : gradientsL_) {
			res.push_back(r.x);
			res.push_back(r.y);
			res.push_back(r.z);
		};
		for (auto& r : gradientsR_) {
			res.push_back(r.x);
			res.push_back(r.y);
			res.push_back(r.z);
		};
		for (auto& r : weights_)
			for (auto& rr : r) res.push_back(rr);

		return std::valarray<double> {res.data(), res.size()};
	};
	virtual void Deserialize(const std::valarray<double>& _values) override {
		values_ = { &_values[0], (size_t)nValues };
		center_ = Vector(_values[nValues], _values[nValues + 1], _values[nValues + 2]);

		gradientsL_.resize(nValues * 3);
		for (int i = 0; i < gradientsL_.size(); i++) {
			gradientsL_[i] = Vector(_values[nValues + 3 * (i + 1)], _values[nValues + 3 * (i + 1) + 1], _values[nValues + 3 * (i + 1) + 2]);
		};
		gradientsR_.resize(nValues * 3);
		for (int i = gradientsR_.size(); i < 2 * gradientsR_.size(); i++) {
			gradientsR_[i - gradientsR_.size()] = Vector(_values[nValues + 3 * (i + 1)], _values[nValues + 3 * (i + 1) + 1], _values[nValues + 3 * (i + 1) + 2]);
		};

		weights_.resize(0);
		int idx = nValues + 3 + 6 * gradientsR_.size();
		for (int i = idx; i < _values.size(); i += 2) {
			weights_.push_back(std::valarray<double> {_values[i], _values[i + 1]});
		};

		return;
	};

	// Update center of Reconstruction if needed
	virtual void RefrashPosition(Vector point) override {
		center_ = point;
	};
};

template<>
WENO2PointsStencil ComputeReconstruction<WENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	//compute Gradients
	int size = value.size();
	std::valarray< Vector > gradientsL(size), gradientsR(size);		//for 1D case only
	std::vector< std::valarray < double> > weights;

	// constant value
	const double eps = 1.0e-16;
	const double d0 = 2.0 / 3.0;		//constant values for UNIFORM grid 2/3
	const double d1 = 1 - d0;

	//1D case
	if (nDim == 1) {
		weights.resize(4);
		for (auto& r : weights) r.resize(size);

		for (int i = 0; i < size; i++) {
			// compute smothness coefficients
			double betta0 = 0.5 * (points[1].x - points[0].x) * (values[1][i] - value[i]) / (points[1].x - point.x);
			betta0 *= betta0;
			betta0 += eps;
			betta0 *= betta0;

			double betta1 = 0.5 * (points[1].x - points[0].x) * (values[0][i] - value[i]) / (points[0].x - point.x);
			betta1 *= betta1;
			betta1 += eps;
			betta1 *= betta1;

			// compute weights for two interpolations (for left and right)
			weights[0][i] = d0 / betta0;
			weights[1][i] = d1 / betta1;
			weights[2][i] = d1 / betta0;
			weights[3][i] = d0 / betta1;

			double a0 = weights[0][i] + weights[1][i];
			double a1 = weights[2][i] + weights[3][i];

			// norm our weights
			weights[0][i] /= a0;
			weights[1][i] /= a0;
			weights[2][i] /= a1;
			weights[3][i] /= a1;

			// compute gradients
			gradientsL[i] = Vector( (value[i] - values[0][i]) / (point.x - points[0].x), 0, 0);
			gradientsR[i] = Vector( (value[i] - values[1][i]) / (point.x - points[1].x), 0, 0);
		};
	};

	//2D case

	//3D case

	return std::move(WENO2PointsStencil(value, point, gradientsL, gradientsR, weights));
};


#endif