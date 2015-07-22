#ifndef TurboStructured_BoundaryConditionConfiguration
#define TurboStructured_BoundaryConditionConfiguration

#include "Vector.h"

enum class BoundaryConditionType {
	Wall,
	SymmetryX,
	SymmetryY,
	SymmetryZ,
	MovingWall,
	General,
	Natural
};

class BoundaryConditionConfiguration {
public:
	BoundaryConditionType BCType;
	Vector Velocity;
	double Gamma;
};

#endif
