#ifndef TurboStructured_BoundaryConditions_BCGeneralBarotropic
#define TurboStructured_BoundaryConditions_BCGeneralBarotropic

#include "BCGeneral.h"
#include <map>
#include <vector>
#include "utility/Vector.h"

namespace BoundaryConditions {

class BCGeneralBarotropic : public BCGeneral {
public:
	//Get dummy cell values
	virtual std::vector<double> getDummyValues(double* values, Vector& faceNormal, Vector& faceCenter, Vector& cellCenter) override {		
		//Compute dummy values
		double ro = values[0];
		double u = values[1]/values[0];
		double v = values[2]/values[0];
		double w = values[3]/values[0];		

		//Get pressure and interpolate internal energy and other variables
		double roDummy = boundaryConditions[BoundaryVariableType::Density].GetDummyValue(ro, faceNormal, faceCenter, cellCenter); 
		double uDummy = boundaryConditions[BoundaryVariableType::VelocityX].GetDummyValue(u, faceNormal, faceCenter, cellCenter);
		double vDummy = boundaryConditions[BoundaryVariableType::VelocityY].GetDummyValue(v, faceNormal, faceCenter, cellCenter);
		double wDummy = boundaryConditions[BoundaryVariableType::VelocityZ].GetDummyValue(w, faceNormal, faceCenter, cellCenter);		

		std::vector<double> res(4);	
		res[0] = roDummy;
		res[1] = roDummy * uDummy;
		res[2] = roDummy * vDummy;
		res[3] = roDummy * wDummy;		
		return res;
	};
 
	virtual void loadConfiguration(BoundaryConditionConfiguration& bcConfig) override {
		boundaryConditions[BoundaryVariableType::Density] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityX] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityY] = CompositeBoundaryConditionInfo();
		boundaryConditions[BoundaryVariableType::VelocityZ] = CompositeBoundaryConditionInfo();

		//Symmetry condition for left and right bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryX) {						
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);						
			return;
		};

		//Symmetry condition for top and bottom bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryY) {				
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);						
			return;
		};

		//Symmetry condition for front and back bounds
		if (bcConfig.BCType == BoundaryConditionType::SymmetryZ) {							
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(0);						
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::Natural) {							
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetNeumanBoundary(0);			
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::Wall) {										
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(0);						
			return;
		};

		if (bcConfig.BCType == BoundaryConditionType::MovingWall) {							
			boundaryConditions[BoundaryVariableType::Density].SetNeumanBoundary(0);
			boundaryConditions[BoundaryVariableType::VelocityX].SetDirichletBoundary(bcConfig.Velocity.x);
			boundaryConditions[BoundaryVariableType::VelocityY].SetDirichletBoundary(bcConfig.Velocity.y);
			boundaryConditions[BoundaryVariableType::VelocityZ].SetDirichletBoundary(bcConfig.Velocity.z);			
			return;
		};

		throw 1;
	}; 

};

};

#endif