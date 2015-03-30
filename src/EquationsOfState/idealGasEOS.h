#ifndef TurboStructured_EquationsOfState_idealGasEOS
#define TurboStructured_EquationsOfState_idealGasEOS

#include "generalEOS.h"


class IdealGasEOS : public GeneralEOS {
public:
	double gamma;

	//Constructor specifies material to use
	IdealGasEOS(double _gamma) {
		gamma = _gamma;
	};

	//set gamma
	inline void SetGamma(double _gamma) {
		gamma = _gamma;
	};

	//Get pressure
	inline double GetPressure(double ro, double e) {
		return (gamma - 1.0)*ro*e;
	};

	//Get sound speed
	inline double GetSoundSpeed(double ro, double e) {
		double pressure = (gamma - 1.0)*ro*e;
		return sqrt(gamma * pressure / ro);
	};

	//Get internal energy
	inline double GetInternalEnergy(double ro, double pressure) {
		return pressure / ((gamma - 1.0) * ro);
	};

	//get partial derivative of pressure by specific internal energy
	inline double GetPressureEnergyDerivative(double ro, double e) {
		return (gamma - 1.0) * ro;
	};

	//get partial derivative of pressure by density
	inline double GetPressureDensityDerivative(double ro, double e) {
		return (gamma - 1.0) * e;
	};
};

#endif