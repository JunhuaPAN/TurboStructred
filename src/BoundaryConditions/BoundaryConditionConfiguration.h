/*!
\file 
This file contains enumeration of boundary condition types and configuration common for them all 
*/

#ifndef TurboStructured_BoundaryConditionConfiguration
#define TurboStructured_BoundaryConditionConfiguration

#include "Vector.h"

/*! Some type of allowed boundary conditions */
enum class BoundaryConditionType {
	Wall,			///< no-slip condition
	SymmetryX,		///< symmetry condition for YZ plane border 
	SymmetryY,		///< symmetry condition for XZ plane border
	SymmetryZ,		///< symmetry condition for XY plane border
	MovingWall,		///< no-clip condition with some constant velocity
	SubsonicInlet,	///< inlet boundary condition
	Natural			///< subsonic outflow condition 
};

/*!	That class accumulates all info about applied boundary condition */
class BoundaryConditionConfiguration {
public:
	BoundaryConditionType BCType;	///<  One of allowed BC types 
	Vector Velocity;				///<  3D vector of Velocity (need for MovingWall type) 
	double Pstatic;					///<  Static pressure
	double Density;					///<  External density
	double Temperature;				///<  External temperature
	Vector Vdirection;				///<  Velocity direction

	// Constructor
	BoundaryConditionConfiguration(BoundaryConditionType _BCType) : BCType(_BCType) {};

	// Default constructor
	BoundaryConditionConfiguration() {};
};

/*!
That class describes if the side of the domain contains just one boundary condition
or several ones (complex case). In first case we set BCMarker at initialization,
and pass the function getMarker in second case.
*/
class DomainBCinfo {
private:
	int bcMarker;	// the number corresponding to the BC surface of the domain 

public:
	bool isComplex;		///< if we have side with several boundary conditions
	std::function<int(const Vector r)> getMarker;	///< function that return Marker value for a side of the domain

	//! Default constructor
	DomainBCinfo() {
		isComplex = false;
	};

	//! Set marker value for trivial case
	inline void SetMarker(int _bcMarker) {
		bcMarker = _bcMarker;
	};

	//! Get marker value for trivial case
	inline int GetMarker() {
		return bcMarker;
	};

};


#endif
