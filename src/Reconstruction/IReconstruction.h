#ifndef TurboStructured_Reconstruction_Reconstruction
#define TurboStructured_Reconstruction_Reconstruction

#include <valarray>
#include "utility/Vector.h"
#include "utility/VariablesTransition.h"

// to define what face the value reconstructed at
enum class CubeFaces {
	xR, 
	xL, 
	yR, 
	yL,
	zR,
	zL
};

struct CellReconstruction {
	std::valarray< double > xR, xL, yR, yL, zR, zL;		// reconstruction storage
};

enum class Reconstruction {
	PiecewiseConstant,
	Linear2PointsStencil,
	ENO2PointsStencil,
};

//Basic class for all reconstruction classes
class IReconstruction {
public:
	int nValues;
	int nDims;
	CellReconstruction recons;

	// Sample solution or get reconstruction at point
	virtual inline std::valarray<double> SampleSolution(Vector const& point) { return{}; }
	virtual inline std::valarray<double> SampleSolution(CubeFaces const& point) {
		switch (point) {
		case CubeFaces::xL:
			return recons.xL;
			break;
		case CubeFaces::xR:
			return recons.xR;
			break;
		case CubeFaces::yL:
			return recons.yL;
			break;
		case CubeFaces::yR:
			return recons.yR;
			break;
		case CubeFaces::zL:
			return recons.zL;
			break;
		case CubeFaces::zR:
			return recons.zR;
			break;
		};
		std::cerr << "can't recognize reconstruction direction";
		std::cout << std::endl;
		return {};
	};
	
	// Default constructor
	IReconstruction() {};
	IReconstruction(CellReconstruction& _recons, int _nDims, int _nValues) :
		recons(_recons),
		nDims(_nDims),
		nValues(_nValues) { };

	// Aditional inner initialization for derived classes
	virtual void Init(int _nValues, int _nDims) {
		nValues = _nValues;
		nDims = _nDims;
	};

	//! Return required
	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV * nD * 2;
	};

	virtual std::valarray<double> Serialize()  {
		std::vector<double> res(0);
		for (auto& r : recons.xL) res.push_back(r);
		for (auto& r : recons.xR) res.push_back(r);
		for (auto& r : recons.yL) res.push_back(r);
		for (auto& r : recons.yR) res.push_back(r);
		for (auto& r : recons.zL) res.push_back(r);
		for (auto& r : recons.zR) res.push_back(r);

		return std::valarray<double> {res.data(), res.size()};
	};
	virtual void Deserialize(const std::valarray<double>& _values)  {
		recons.xL = { &_values[0], (size_t)nValues };
		recons.xR = { &_values[nValues], (size_t)nValues };
		
		// 2D case
		if (nDims > 1) {
			recons.yL = { &_values[2 * nValues], (size_t)nValues };
			recons.yR = { &_values[3 * nValues], (size_t)nValues };
		};

		// 3D case
		if (nDims > 2) {
			recons.zL = { &_values[4 * nValues], (size_t)nValues };
			recons.zR = { &_values[5 * nValues], (size_t)nValues };
		};

		return;
	};
};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim) {
	throw std::runtime_error("We dont have required function.");
};

#endif
