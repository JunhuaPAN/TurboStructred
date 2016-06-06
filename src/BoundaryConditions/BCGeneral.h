#ifndef TurboStructured_BoundaryConditions_BCGeneral
#define TurboStructured_BoundaryConditions_BCGeneral

#include <map>
#include <vector>
#include "utility/Vector.h"
#include "GasProperties.h"

//Compressible Grid scheme
//
//	_|________|_________|_
//	 |        |         |
//	 |        |         |
//	 |        |         |
//	 |        |         |					Second layer of Inner cells
//	 |        |         |
//	_|________|_________|_
//	 |        |         |		)
//	 |        |         |		> H			First layer of Inner cells
//	_|________|_________|_		)	
//	______________________		<-- Border is Here
//	 |        |         |		)
//	 |        |         |		> H			First lauer of Fictious cells
//	_|________|_________|_		)
//	 |        |         |		)
//	 |        |         |		> H			Second layer of Fictious cells
//	_|________|_________|_		)
//	 |		  |			|

//! Boundary conditions information for one side of the full domain
class BCinfo {
private:
	bool isTrivial;		// true if the whole side corresponds to only one boundary condition
	int trivialKey;		// Marker value for trivial case

public:
	std::function<int(const Vector& r)> computeMarker;		// function that returns BC marker

	// Trivial case constructor
	BCinfo(int bcMarker) : trivialKey(bcMarker) {
		isTrivial = true;
	};

	// Complex case constructor
	BCinfo(std::function<int(const Vector& r)> _computeMarker) : computeMarker(_computeMarker) {
		isTrivial = false;
	};

	// Default constructor
	BCinfo() {};

	// Get marker
	inline int getMarker(const Vector& r) {
		if (isTrivial) return trivialKey;
		return computeMarker(r);
	};

};

namespace BoundaryConditions {

//! Basic class for all boundary conditions
class BCGeneral{
public:
	// number of variables
	const int nVar{ 5 };

	//Gas model parameters
	std::unique_ptr<GasProperties> gas_prop;

	//! Set gas properties
	void BindGasProperties(GasProperties& _gas_prop) {
		gas_prop = std::make_unique<GasProperties>(_gas_prop);
	};

	//! Load configuration if need
	virtual void loadConfiguration(const BoundaryConditionConfiguration& config) {};

	//! Get dummy cell values
	virtual std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		return std::valarray<double>{};
	};

	//! Get dummy reconstructions
	virtual std::valarray<double> getDummyReconstructions(double* values, Vector faceNormal) {
		return std::move(getDummyValues(values, faceNormal, Vector(0, 0, 0), Vector(0, 0, 0)));
	};

	//! Get conservative variables on the border
	virtual std::valarray<double> getFaceValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		auto inVal = std::valarray<double>(values, nVar);
		auto duVal = getDummyValues(values, faceNormal, faceCenter, cellCenter);
		return 0.5 * (inVal + duVal);
	};
};

//! Describe all boundary conditions

//! Subsonic inlet BC
class SubsonicInlet : public BCGeneral {
public:
	double Pout, Rhout;
	Vector Vdir;		// unity vector of velocity direction in dummy cell
	Vector V;			// vector of velosity

	//! Get dummy cell values
	virtual std::valarray<double> getDummyValues2(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		// Compute dummy values
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		double E = values[4] / values[0];
		double roe = ro * E - 0.5 * ro * (u * u + v * v + w * w);
		double P = (gas_prop->gamma - 1.0) * roe;

		// Compute R- invariant
		double c = sqrt(gas_prop->gamma * P / ro);
		Vector Vel{ u, v, w };
		double Rm = Vel * faceNormal - 2.0 * c / (gas_prop->gamma - 1);
		
		// Compute face state
		double Pf = Pout;
		double rof = Rhout;
		double cf = sqrt(gas_prop->gamma * Pf / rof);
		Vector Velf = (Rm + 2.0 * cf / (gas_prop->gamma - 1)) * Vdir;
		Velf /= (Vdir * faceNormal);

		// Extrapolate in dummy cell
		double PDummy = P;
		double roDummy = Rhout;
		double roeDummy = PDummy / (gas_prop->gamma - 1);
		double cDummy = sqrt(gas_prop->gamma * PDummy / roDummy);
		Rm += 2.0 * cDummy / (gas_prop->gamma - 1);
		double cos_a = Vdir * faceNormal;

		// Compute dummy velocity vector
		Vector VDummy = (Rm / cos_a) * Vdir;

		// return conservative form
		std::valarray<double> res(5);
		res[0] = roDummy;
		res[1] = roDummy * VDummy.x;
		res[2] = roDummy * VDummy.y;
		res[3] = roDummy * VDummy.z;
		res[4] = roeDummy + 0.5 * roDummy * (VDummy * VDummy);
		return res;
	};
	virtual std::valarray<double> getDummyValues1(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter)  {
		// Compute dummy values
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		double E = values[4] / values[0];
		double roe = ro * E - 0.5 * ro * (u * u + v * v + w * w);
		double P = (gas_prop->gamma - 1.0) * roe;

		// Compute R- invariant
		double c = sqrt(gas_prop->gamma * P / ro);
		Vector Vel{ u, v, w };
		double Rm = Vel * faceNormal - 2.0 * c / (gas_prop->gamma - 1);

		// Compute dummy state
		double cDummy = 0.5 * (V * faceNormal - Rm) * (gas_prop->gamma - 1);
		double PDummy = Pout;
		double roDummy = gas_prop->gamma * PDummy / (cDummy * cDummy);
		double roeDummy = PDummy / (gas_prop->gamma - 1);

		// return conservative form
		std::valarray<double> res(5);
		res[0] = roDummy;
		res[1] = roDummy * V.x;
		res[2] = roDummy * V.y;
		res[3] = roDummy * V.z;
		res[4] = roeDummy + 0.5 * roDummy * (V * V);
		return res;
	};
	virtual std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) override  {
		
		// return conservative form
		std::valarray<double> res(5);
		res[0] = Rhout;
		res[1] = Rhout * V.x;
		res[2] = Rhout * V.y;
		res[3] = Rhout * V.z;
		res[4] = Pout / (gas_prop->gamma - 1) + 0.5 * Rhout * (V * V);
		return res;
	};
};

//! Sunsonic outlet BC
class SubsonicOutlet : public BCGeneral {

};

//! Values in dummy cells and at the face are the same as internal
class Natural : public BCGeneral {
	std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		return std::valarray<double>(values, nVar);
	};
};

//! Symmetry BC
class SymmetryY : public BCGeneral {

};

// No slip condition (velocity is zero at the wall)
class Wall : public BCGeneral {
	// Dummy values for no slip condition
	std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		auto res = std::valarray<double>(values, nVar);
		res[1] *= (-1);
		res[2] *= (-1);
		res[3] *= (-1);
		return res;
	};
};

// No slip condition (velocity is constant at the wall)
class MovingWall : public BCGeneral {
private:
	Vector Vb;		// Valocity at the wall

public:
	void loadConfiguration(const BoundaryConditionConfiguration& config) override {
		Vb = config.Velocity;
	};

	// Dummy values for no slip condition
	std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		std::valarray<double> res(nVar);
		auto rho = values[0];
		auto Vin = Vector(values[1], values[2], values[3]) / rho;
		auto rhoe = values[4] - 0.5 * rho * (Vin * Vin);
		
		// Extrapolate velocity
		Vector Vd = 2.0 * Vb - Vin;

		// Create and return conservative variables
		res[0] = rho;
		res[1] = rho * Vd.x;
		res[2] = rho * Vd.y;
		res[3] = rho * Vd.z;
		res[4] = rhoe + 0.5 * rho * (Vd * Vd);
		return res;
	};
};

//! Function that creates pointer to the BCGeneral class
std::unique_ptr<BCGeneral> CreateBC(BoundaryConditionConfiguration& myCondition) {
	std::unique_ptr<BCGeneral> res;

	switch (myCondition.BCType) {
	case BoundaryConditionType::SubsonicInlet:
		res = std::make_unique<SubsonicInlet>();
		break;
	case BoundaryConditionType::Natural:
		res = std::make_unique<Natural>();
		break;
	case BoundaryConditionType::Wall:
		res = std::make_unique<Wall>();
		break;
	case BoundaryConditionType::MovingWall:
		res = std::make_unique<MovingWall>();
		break;
	default:
		std::cout << "Can't find appropriate Boundary Condition!" << std::endl;
		break;
	};
	return res;
};

};

#endif
