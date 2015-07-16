#ifndef TurboStructured_Reconstruction_WENO2PointsStencil
#define TurboStructured_Reconstruction_WENO2PointsStencil

#include <valarray>
#include <vector>
#include "utility\Vector.h"
#include <Reconstruction\IReconstruction.h>


//Reconstruction by convex combination of 2 points stencil reconstructions
class WENO2PointsStencil : public IReconstruction {
protected:
	//std::valarray<double> values_;
	Vector center_;
	std::valarray< Vector > gradientsL_;
	std::valarray< Vector > gradientsR_;
	std::valarray< double > right_reconst_;
	std::valarray< double > left_reconst_;
	//double dl0, dl1;		//coefficients of weights for left and right stencil
	//double betta0, betta1;	// smothness indicators for left and right stencil : betta = (eps + br)^2	

public:
	//WENO 3 order reconstruction
	virtual inline std::valarray<double> SampleSolution(Vector const& point) {
		//right face 
		if (point.x > center_.x) return right_reconst_;
		//left face
		return left_reconst_;
	};

	//constructor
	WENO2PointsStencil() { };
	WENO2PointsStencil(std::valarray<double>& _left_reconst, std::valarray<double>& _right_reconst, Vector& _center) :
		left_reconst_(_left_reconst),
		center_(_center),
		right_reconst_(_right_reconst)
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
WENO2PointsStencil ComputeReconstruction<WENO2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector point, int nDim) {
	//compute Gradients
	int size = value.size();
	std::valarray< Vector > gradientsL(size), gradientsR(size);

	//result right and left reconstructions
	std::valarray<double> left_reconst(size), right_reconst(size);

	//right and left WENO reconstructions for elements of value
	//double recR = 0;
	//double recL = 0;

	//1D case
	if (nDim == 1) {
		std::valarray<double> wr0, wr1, wl0, wl1;		//weights
		//std::valarray<double> grad_l(size), grad_r(size);
		double grad_l, grad_r;
		for (int i = 0; i < size; i++) {
			grad_l = (value[i] - values[0][i]) / (point.x - points[0].x);
			grad_r = (values[1][i] - value[i]) / (points[1].x - point.x);
			gradientsL[i] = Vector(grad_l, 0, 0);
			gradientsR[i] = Vector(grad_r, 0, 0);
		};
	};

	//2D case

	//3D case

	//ENO2PointsStencil res = ENO2PointsStencil(value, point, gradients);
	return std::move(WENO2PointsStencil(left_reconst, right_reconst, point));
};


#endif