#ifndef TurboStructured_Reconstruction_ENO2PointsStencil
#define TurboStructured_Reconstruction_ENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by piecewise constant approximation
class ENO2PointsStencil : public IReconstruction {
protected:
	std::valarray<double> values_;
	Vector center_;
	std::valarray< Vector > gradients_;

public:
	//Piecewise reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		int&& size = gradients_.size();
		std::valarray<double> res(size);

		Vector dr = (point - center_);
		for (int i = 0; i < size; i++) res[i] = dr * gradients_[i];
		return std::move(values_ + res);
	};

	//constructor
	ENO2PointsStencil() { };
	ENO2PointsStencil(std::valarray<double>& _values, Vector& _center, std::valarray< Vector >& _gradients):
		values_(_values),
		center_(_center),
		gradients_(_gradients)
		{ };

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
ENO2PointsStencil ComputeReconstruction<ENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector point, int nDim) {
	//compute Gradients
	int&& size = value.size();
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

		//compute convex hull of gradients
	};

	//2D case

	//3D case
	
	//ENO2PointsStencil res = ENO2PointsStencil(value, point, gradients);
	return std::move(ENO2PointsStencil(value, point, gradients));
};


#endif