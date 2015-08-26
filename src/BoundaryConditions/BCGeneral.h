#ifndef TurboStructured_BoundaryConditions_BCGeneral
#define TurboStructured_BoundaryConditions_BCGeneral

#include "IBoundaryCondition.h"
#include <map>
#include <vector>
#include "utility/Vector.h"

namespace BoundaryConditions {

enum class BoundaryVariableType {
	Pressure,
	Density,
	VelocityX,
	VelocityY,
	VelocityZ,
	VelocityNormal,
	VelocityTangential,
	InternalEnergy,
	VolumeFraction
};

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

class CompositeBoundaryConditionInfo {
public:
	double CValue;
	double CGradient;
	double Value;

	//Set values
	void SetValues(double cValue, double cGradient, double value) {
		if (cValue > 0 || cGradient > 0) {
			CValue = cValue;
			CGradient = cGradient;
			Value = value;
			return;
		};

		//throw new Exception("One of coefficients must be non-zero");
	};

	//Set dirichlet boundary condition
	void SetDirichletBoundary(double Value) {
		SetValues(1.0, 0.0, Value);
	};

	//Set neuman boundary condition
	void SetNeumanBoundary(double Gradient) {
		SetValues(0.0, 1.0, Gradient);
	};

	//Interpolate
	inline double GetDummyValue(const double Vin, const Vector& faceNormal, const Vector& faceCenter, const Vector& cellCenter) {
		double dn = (faceCenter - cellCenter) * faceNormal;
		double Vout = Vin * (CGradient - CValue) + 2.0 * Value * (CGradient * dn + CValue);
		return Vout;
	};

	//Interpolate
	Vector GetDummyGradient(const double inV, const Vector inGrad, const Vector& faceNormal, const Vector& faceCenter, const Vector& cellCenter) {
		Vector res = inGrad;
		Vector GradNormal = (inGrad * faceNormal) * faceNormal;		//normal part of value gradient in inner cell
		Vector GradTangential = res - GradNormal;					//tangential to border part of value gradient vector

		//if we know normal part of value (inV) gradient
		if(CGradient != 0)
		{
			double sign = (faceNormal * Vector(1, 1, 1));
			res = GradTangential + (2.0 * sign * Value - GradNormal.mod()) * faceNormal;
		};

		//if we know value on the border
		if(CValue != 0) {
			res = GradNormal - GradTangential;
		};

		return res;
	};
};

class BCGeneral : public IBoundaryCondition {
public:
	//Gas model parameters
	double gamma;

	//Boundary conditions in general form
	std::map<BoundaryVariableType, CompositeBoundaryConditionInfo> boundaryConditions;

	// Get dummy cell values
	virtual std::valarray<double> getDummyValues(double* values, Vector& faceNormal, Vector& faceCenter, Vector& cellCenter) {		
		//Compute dummy values
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		double E = values[4] / values[0];
		double roe = ro * E - ro * (u*u + v*v + w*w) / 2.0;

		//Get pressure and interpolate internal energy and other variables
		double P = (gamma - 1.0) * roe;
		double PDummy = boundaryConditions[BoundaryVariableType::Pressure].GetDummyValue(P, faceNormal, faceCenter, cellCenter);
		double roDummy = boundaryConditions[BoundaryVariableType::Density].GetDummyValue(ro, faceNormal, faceCenter, cellCenter); 
		double uDummy = boundaryConditions[BoundaryVariableType::VelocityX].GetDummyValue(u, faceNormal, faceCenter, cellCenter);
		double vDummy = boundaryConditions[BoundaryVariableType::VelocityY].GetDummyValue(v, faceNormal, faceCenter, cellCenter);
		double wDummy = boundaryConditions[BoundaryVariableType::VelocityZ].GetDummyValue(w, faceNormal, faceCenter, cellCenter);
		double roeDummy = PDummy / (gamma - 1);

		std::valarray<double> res(5);
		res[0] = roDummy;
		res[1] = roDummy * uDummy;
		res[2] = roDummy * vDummy;
		res[3] = roDummy * wDummy;
		res[4] = roeDummy + roDummy * (uDummy*uDummy + vDummy*vDummy + wDummy*wDummy)/2.0;
		return res;
	};

	// Get dummy reconstructions
	virtual std::valarray<double> getDummyReconstructions(double* values) {
		return getDummyValues(values, Vector(0, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0));
	};
 
	virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {
		//Symmetry condition for left and right bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryX) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		//Symmetry condition for top and bottom bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryY) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		//Symmetry condition for front and back bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryZ) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::Natural) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::Wall) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::MovingWall) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(bcConfig.Velocity.x);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(bcConfig.Velocity.y);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(bcConfig.Velocity.z);
			boundaryConditions[BoundaryVariableType::Pressure].SetNeumanBoundary(0);

			gamma = bcConfig.Gamma;
			return;
		};

		throw 1;
	}; 

};

};

#endif