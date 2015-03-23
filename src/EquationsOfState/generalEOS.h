#ifndef TurboStructured_EquationsOfState_generalEOS
#define TurboStructured_EquationsOfState_generalEOS

#include <vector>

class GeneralEOS {

public:
	//Get pressure
	virtual double GetPressure(double ro, double e) = 0;

	//Get sound speed
	virtual double GetSoundSpeed(double ro, double e) = 0;

	//get internal energy
	virtual double GetInternalEnergy(double ro, double pressure) = 0;

	//get partial derivative of pressure by specific internal energy
	virtual double GetPressureEnergyDerivative(double ro, double e) = 0;

	//get partial derivative of pressure by density
	virtual double GetPressureDensityDerivative(double ro, double e) = 0;
};


#endif