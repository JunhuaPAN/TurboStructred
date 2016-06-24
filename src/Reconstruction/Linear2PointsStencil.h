#ifndef TurboStructured_Reconstruction_LinearPointsStencil
#define TurboStructured_Reconstruction_LinearPointsStencil

#include <valarray>
#include <vector>
#include "utility/Vector.h"
#include <Reconstruction/IReconstruction.h>


// Linear reconstruction by 2 x 2 x 2 points stencil (i+-1) (j+-1) (k+-1)
class Linear2PointsStencil : public IReconstruction {
public:
	std::vector<Vector> grads;		// gradients of conservative variables at point (0, 0, 0)

	//constructors
	Linear2PointsStencil() {};
	Linear2PointsStencil(std::valarray<double> _vals, int _nDims, int _nValues, std::vector<Vector> _grads) :
		IReconstruction(_vals, _nDims, _nValues), grads(_grads) {};

	// Initialization
	virtual void Init(int _nValues, int _nDims) {
		IReconstruction::Init(_nValues, _nDims);
		grads.resize(_nValues);
	};

	virtual std::valarray<double> SampleSolution(Vector const& point) override {
		std::valarray<double> res(nValues);
		for(int i = 0; i < nValues; i++) res[i] = vals[i] + grads[i] * point;
		return res;
	};

	// Serrialization part
	static std::size_t GetBufferLenght(int nD, int nV) {
		return (1 + nD) * nV ;
	};
	virtual std::valarray<double> Serialize() override {
		std::valarray<double> res(GetBufferLenght(nDims, nValues));
		res[std::slice(0, nValues, 1)] = vals;	// first write values

		// then gradients are writen
		auto ibuf = nValues;
		for (int i = 0; i < nValues; i++) res[ibuf + i] = grads[i].x;
		ibuf += nValues;
		if (nDims > 1) for (int i = 0; i < nValues; i++) res[ibuf + i] = grads[i].y;
		ibuf += nValues;
		if (nDims > 2) for (int i = 0; i < nValues; i++) res[ibuf + i] = grads[i].z;

		return res;
	};
	virtual void Deserialize(const std::valarray<double>& _values) override {
		vals = _values[std::slice(0, nValues, 1)]; // catch the values first

		// then read the gradients
		auto ibuf = nValues;
		for (int i = 0; i < nValues; i++) grads[i].x = _values[ibuf + i];
		ibuf += nValues;
		if (nDims > 1) for (int i = 0; i < nValues; i++) grads[i].y = _values[ibuf + i];
		ibuf += nValues;
		if (nDims > 2) for (int i = 0; i < nValues; i++) grads[i].z = _values[ibuf + i];

		return;
	};
};

template<>
Linear2PointsStencil ComputeReconstruction<Linear2PointsStencil>(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	size_t size = value.size();
	std::vector<Vector> gradients(size);

	// 1D case
	std::valarray<double> pds;	// partial derivatives of reconstructed variables
	pds = (values[1] - values[0]) / (points[1].x - points[0].x);
	for (int i = 0; i < size; i++) {
		gradients[i].x = pds[i];
	};

	//2D case
	if (nDim > 1) {
		auto dy = points[3].y - points[2].y;
		pds = (values[3] - values[2]) / dy;
		for (int i = 0; i < size; i++) {
			gradients[i].y = pds[i];
		};
	};

	//3D case
	if (nDim > 2) {
		auto dz = points[5].z - points[4].z;
		pds = (values[5] - values[4]) / dz;
		for (int i = 0; i < size; i++) {
			gradients[i].z = pds[i];
		};
	};

	// Create reconstruction
	return std::move(Linear2PointsStencil(value, nDim, size, gradients));
};


#endif
