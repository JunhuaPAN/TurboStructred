/*!
\file 
This file contains enumeration of boundary condition types and configuration common for them all 
*/

#ifndef TurboStructured_BoundaryConditionConfiguration
#define TurboStructured_BoundaryConditionConfiguration

#include "Vector.h"

/*! Some type of allowed boundary conditions */
enum class BoundaryConditionType {
	Wall,		///< no-slip condition
	SymmetryX,	///< symmetry condition for YZ plane border 
	SymmetryY,	///< symmetry condition for XZ plane border
	SymmetryZ,	///< symmetry condition for XY plane border
	MovingWall,	///< no-clip condition with some constant velocity
	General,	///< general type??
	Natural		///< subsonic outflow condition 
};

/*!	That class accumulates all info about applied boundary condition */
class BoundaryConditionConfiguration {
public:
	BoundaryConditionType BCType;	///<  One of allowed BC types 
	Vector Velocity;				///<  3D vector of Velocity (need for MovingWall type) 
	double Gamma;					///<  Medium property
};

#endif
