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
	//Gas model parameters
	std::unique_ptr<GasProperties> gas_prop;

	//! Load configuration
	virtual void loadConfiguration(const BoundaryConditionConfiguration& config) {};

	//! Get dummy cell values
	virtual std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {
		return std::valarray<double>{};
	};

	//! Get dummy reconstructions
	virtual std::valarray<double> getDummyReconstructions(double* values, Vector faceNormal) {
		return std::valarray<double>{};
	};
};

//! Describe all boundary conditions
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

class SubsonicOutlet : public BCGeneral {

};

class Natural : public BCGeneral {

};

class SymmetryY : public BCGeneral {

};

class Wall : public BCGeneral {

};

//! Function that creates pointer to the BCGeneral class
std::unique_ptr<BCGeneral> CreateBC(BoundaryConditionConfiguration& myCondition) {
	std::unique_ptr<BCGeneral> res;

	switch (myCondition.BCType) {
	case BoundaryConditionType::SubsonicInlet:
		res = std::make_unique<SubsonicInlet>();
		break;
	default:
		std::cout << "Can't find appropriate Boundary Condition!" << std::endl;
		break;
	};
	return res;
};

};

#endif
