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
	ENO2PointsStencil,
	WENO2PointsStencil,
	ENO2CharactVars
};

//Basic class for all reconstruction classes
class IReconstruction {
public:
	int nValues;
	int nDims;
	CellReconstruction recons;

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
		return{};
	};

	// Default constructor
	IReconstruction() {
		nValues = 0;
	};
	IReconstruction(CellReconstruction& _recons, int _nDims, int _nValues) :
		recons(_recons),
		nDims(_nDims),
		nValues(_nValues) { };

	//! Return required
	// Serrialization
	static std::size_t GetBufferLenght(int nD, int nV) {
		return nV * nD * 2 + 2;
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
		recons.xR = { &_values[0], (size_t)nValues };
		recons.xL = { &_values[nValues], (size_t)nValues };
		
		// 2D case
		if (nDims > 1) {
			recons.yR = { &_values[2 * nValues], (size_t)nValues };
			recons.yL = { &_values[3 * nValues], (size_t)nValues };
		};

		// 3D case
		if (nDims > 2) {
			recons.zR = { &_values[2 * nValues], (size_t)nValues };
			recons.zL = { &_values[3 * nValues], (size_t)nValues };
		};

		return;
	};
};

template<typename T>
T ComputeReconstruction(std::vector<std::valarray<double> > values, std::vector<Vector> points, std::valarray<double> value, Vector& point, int nDim, double gamma) {
	//static_assert(false, "We dont have required function.");
	throw std::runtime_error("We dont have required function.");
};

#endif
