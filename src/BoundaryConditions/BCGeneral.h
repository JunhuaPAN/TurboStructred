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
	double Vreff;	// refference value

	//Set values
	inline void SetValues(double cValue, double cGradient, double value, double valueref) {
		CValue = cValue;
		CGradient = cGradient;
		Value = value;
		Vreff = valueref;
	};

	//Set dirichlet boundary condition
	void SetDirichletBoundary(double Value) {
		SetValues(1.0, 0.0, Value, 0);
	};

	//Set neuman boundary condition
	void SetNeumanBoundary(double Gradient) {
		SetValues(0.0, 1.0, Gradient, 0);
	};

	//Interpolate
	inline double GetDummyValue(const double Vin, const Vector& faceNormal, const Vector& faceCenter, const Vector& cellCenter) {
		double dn = (faceCenter - cellCenter) * faceNormal;
		double Vout = Vin * (CGradient - CValue) + 2.0 * Value * (CGradient * dn + CValue);
		return Vout + Vreff;
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

	//! Get dummy cell values
	virtual std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) {		
		// Compute dummy values
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		double E = values[4] / values[0];
		double roe = ro * E - 0.5 * ro * (u * u + v * v + w * w);

		// Get pressure and interpolate internal energy and other variables
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
		res[4] = roeDummy + 0.5 * roDummy * (uDummy * uDummy + vDummy * vDummy + wDummy * wDummy);
		return res;
	};

	//! Get dummy reconstructions
	virtual std::valarray<double> getDummyReconstructions(double* values) {
		return getDummyValues(values, Vector(0, 0, 0), Vector(0, 0, 0), Vector(0, 0, 0));
	};
 
	virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {
		gamma = bcConfig.Gamma;

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
						
			return;
		};

		throw 1;
	}; 

};

class SubsonicInletBC : public BCGeneral {
public:
	double Pout, Rhout;
	Vector Vdir;		// unity vector of velocity direction in dummy cell

	//! Get dummy cell values
	virtual std::valarray<double> getDummyValues(double* values, Vector faceNormal, Vector faceCenter, Vector cellCenter) override {
		// Compute dummy values
		double ro = values[0];
		double u = values[1] / values[0];
		double v = values[2] / values[0];
		double w = values[3] / values[0];
		double E = values[4] / values[0];
		double roe = ro * E - 0.5 * ro * (u * u + v * v + w * w);
		double P = (gamma - 1.0) * roe;

		// Compute R- invariant
		double c = sqrt(gamma * P / ro);
		Vector Vel{ u, v, w };
		double Rm = Vel * faceNormal - 2.0 * c / (gamma - 1);
		
		// Compute dummy state
		double PDummy = Pout;
		double roDummy = Rhout;
		double roeDummy = PDummy / (gamma - 1);
		double cDummy = sqrt(gamma * PDummy / roDummy);
		Rm += 2.0 * cDummy / (gamma - 1);
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

	virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {
		// Set output parameters
		Pout = bcConfig.Pstatic;
		Rhout = bcConfig.Density;
		Vdir = bcConfig.Vdirection;
		Vdir /= Vdir.mod();

		// create boundary conditions with zero velocity gradients
		boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
		boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);

		return;
	};

};

};

#endif
