#ifndef TurboStructured_EquationsOfState_BaratropicEOS
#define TurboStructured_EquationsOfState_BaratropicEOS

#include "generalEOS.h"


class BaratropicEOS : public GeneralEOS {
public:
	double E;	//modul Yunga
	double P0;		//initial pressure
	double Ro0;		//initial dencity

	//Constructor specifies material to use
	BaratropicEOS(double _E, double _P0, double _Ro0) : E(_E), P0(_P0), Ro0(_Ro0) {};

	//set gamma
	inline void SetE(double _E) {
		E = _E;
	};

	//Get pressure
	inline double GetPressure(double ro, double e) {
		return P0 + E * ( ro - Ro0) / Ro0;
	};

	//Get sound speed
	inline double GetSoundSpeed(double ro, double e) {
		return sqrt(E/Ro0);
	};

	//Get internal energy
	inline double GetInternalEnergy(double ro, double pressure) {
		return 0;
	};

	//get partial derivative of pressure by specific internal energy
	inline double GetPressureEnergyDerivative(double ro, double e) {
		return 0;
	};

	//get partial derivative of pressure by density
	inline double GetPressureDensityDerivative(double ro, double e) {
		return E/Ro0;
	};
};

#endif