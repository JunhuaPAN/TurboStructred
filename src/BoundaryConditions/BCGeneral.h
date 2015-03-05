#ifndef TurboStructured_BoundaryConditions_BCGeneral
#define TurboStructured_BoundaryConditions_BCGeneral

#include "BoundaryCondition.h"
#include <map>
#include <vector>
#include "utility/Vector.h"

enum class BoundaryVariableType {
	Pressure,
	Density,
	VelocityX,
	VelocityY,
	VelocityZ,
	InternalEnergy
};

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
	double GetDummyValue(const double inV, const Vector& faceNormal, const Vector& faceCenter, const Vector& cellCenter) {
		double dn = -(cellCenter - faceCenter) * faceNormal;
		double a = (Value - CValue * inV) / (CValue + CGradient / dn);
		double b = inV + a;
		double value = a + b;
		return value;
	};

	//Interpolate
	Vector GetDummyGradient(const double inV, const Vector inGrad, const Vector& faceNormal, const Vector& faceCenter, const Vector& cellCenter) {
		double dn = -(cellCenter - faceCenter) * faceNormal;
		double a = (Value - CValue * inV) / (CValue + CGradient / dn);
		double b = inV + a;
		double value = a + b;

		Vector outGrad = inGrad;
		double dudnIn = inGrad * faceNormal;
		double dudnOut = (2.0 * (inV - Value) / dn - dudnIn) * CValue + Value * CGradient;
		outGrad = inGrad + (dudnOut - dudnIn) * faceNormal;

		return outGrad;
	};
};

class BCGeneral : public BoundaryCondition {
public:
	//Gas model parameters
	double gamma;

	//Boundary conditions in general form
	std::map<BoundaryVariableType, CompositeBoundaryConditionInfo> boundaryConditions;

	//Get dummy cell values
	std::vector<double> getDummyValues(double* values, Vector& faceNormal, Vector& faceCenter, Vector& cellCenter) {		
		//Compute dummy values
		double ro = values[0];
		double u = values[1]/values[0];
		double v = values[2]/values[0];
		double w = values[3]/values[0];
		double E = values[4]/values[0];
		double roe = ro*E - ro*(u*u + v*v + w*w)/2.0;

		//Get pressure and interpolate internal energy and other variables
		double P = (gamma - 1.0) * roe;
		double PDummy = boundaryConditions[BoundaryVariableType::Pressure].GetDummyValue(P, faceNormal, faceCenter, cellCenter);
		double roDummy = boundaryConditions[BoundaryVariableType::Density].GetDummyValue(ro, faceNormal, faceCenter, cellCenter); 
		double uDummy = boundaryConditions[BoundaryVariableType::VelocityX].GetDummyValue(u, faceNormal, faceCenter, cellCenter);
		double vDummy = boundaryConditions[BoundaryVariableType::VelocityY].GetDummyValue(v, faceNormal, faceCenter, cellCenter);
		double wDummy = boundaryConditions[BoundaryVariableType::VelocityZ].GetDummyValue(w, faceNormal, faceCenter, cellCenter);
		double roeDummy = PDummy / (gamma - 1);

		std::vector<double> res(5);	
		res[0] = roDummy;
		res[1] = roDummy * uDummy;
		res[2] = roDummy * vDummy;
		res[3] = roDummy * wDummy;
		res[4] = roeDummy + roDummy*(uDummy*uDummy + vDummy*vDummy + wDummy*wDummy)/2.0;
		return res;
	};

	void loadConfiguration(BoundaryConditionConfiguration& bcConfig) {
		if (bcConfig.BCType == BoundaryConditionType::Symmetry) {				
			boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();
			boundaryConditions[BoundaryVariableType::Pressure] = CompositeBoundaryConditionInfo();

			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::Pressure].SetDirichletBoundary(0);

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

		throw 1;
	}; 

};

#endif