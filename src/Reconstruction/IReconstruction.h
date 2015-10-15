#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility/Vector.h"

// to define what face the value reconstructed at
enum class CubeFaces {
	xR, 
	xL, 
	yR, 
	yL,
	zR,
	zL
};

enum class Reconstruction {
	PiecewiseConstant,
	ENO2PointsStencil,
	WENO2PointsStencil
};

//Basic class for all reconstruction classes
class IReconstruction {
public:
	int nValues;
	int nDims;
	std::valarray< double > xR, xL, yR, yL, zR, zL;		// reconstruction storage

	virtual inline std::valarray<double> SampleSolution(CubeFaces const& point) {
		switch (point) {
		case CubeFaces::xL:
			return xL;
			break;
		case CubeFaces::xR:
			return xR;
			break;
		case CubeFaces::yL:
			return yL;
			break;
		case CubeFaces::yR:
			return yR;
			break;
		case CubeFaces::zL:
			return zL;
			break;
		case CubeFaces::zR:
			return zR;
			break;
		};
		std::cerr << "can't recognize reconstruction direction";
		std::cout << std::endl;
		return{};
	};

	// Default constructor
	IReconstruction() {
		nValues = 0;
	};
	IReconstruction(std::valarray<double>& _xR, std::valarray<double>& _xL, std::valarray<double>& _yR, std::valarray<double>& _yL, std::valarray<double>& _zR, std::valarray<double>& _zL, int _nDims, int _nValues) :
		xR(_xR),
		xL(_xL),
		yR(_yR),
		yL(_yL),
		zR(_zR),
		zL(_zL),
		nDims(_nDims),
		nValues(_nValues) { };

	//! Return required
	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV * nD * 2;
	};
	virtual std::valarray<double> Serialize()  {
		std::vector<double> res(0);
		for (auto& r : xL) res.push_back(r);
		for (auto& r : xR) res.push_back(r);
		for (auto& r : yL) res.push_back(r);
		for (auto& r : yR) res.push_back(r);
		for (auto& r : zL) res.push_back(r);
		for (auto& r : zR) res.push_back(r);

		return std::valarray<double> {res.data(), res.size()};
	};
	virtual void Deserialize(const std::valarray<double>& _values)  {
		xR = { &_values[0], (size_t)nValues };
		xL = { &_values[nValues], (size_t)nValues };
		
		// 2D case
		if (nDims > 1) {
			yR = { &_values[2 * nValues], (size_t)nValues };
			yL = { &_values[3 * nValues], (size_t)nValues };
		};

		// 3D case
		if (nDims > 2) {
			zR = { &_values[2 * nValues], (size_t)nValues };
			zL = { &_values[3 * nValues], (size_t)nValues };
		};

		return;
	};

};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	//static_assert(false, "We dont have required function.");
	throw std::runtime_error("We dont have required function.");
};

#endif
