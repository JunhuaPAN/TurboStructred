#ifndef TurboStructured_EquationsOfState_VOFMixture_BaratropicTwoPhaseEOS
#define TurboStructured_EquationsOfState_VOFMixture_BaratropicTwoPhaseEOS

#include "EquationsOfState/generalEOS.h"

class BarotropicTwoPhaseEOS : public GeneralEOS {
private:
	double a1, a2, b1, b2;		// Pi = ai*Rho + bi

public:
	//Constructor specifies material to use
	BarotropicTwoPhaseEOS(double E1, double E2, double P01, double P02, double Ro1, double Ro2) {
		// squares of individual sound speeds
		a1 = E1 / Ro1;
		a2 = E2 / Ro2;
		// free member
		b1 = P01 - E1;
		b2 = P02 - E2;
	};

	//set gamma
	inline void SetE(double E1, double E2, double P01, double P02, double Ro1, double Ro2) {
		// squares of individual sound speeds
		a1 = E1 / Ro1;
		a2 = E2 / Ro2;
		// free member
		b1 = P01 - E1;
		b2 = P02 - E2;
	};

	//Get pressure from total density and volume fraction of first medium
	inline double GetPressure(double ro, double alpha) {
		double alpha2 = 1 - alpha;
		double a = a1 * a2 / (a2 * alpha + a1 * alpha2);	//common square of sound speed
		double b = ( b1 * a2 * alpha + b2 * a1 * alpha2 ) / (a2 * alpha + a1 * alpha2);
		return a * ro + b;
	};

	//Get sound speed
	inline double GetSoundSpeed(double ro, double alpha) {
		double a = a1 * a2 / (a2 * alpha + a1 * ( 1.0 - alpha ));
		return sqrt(a);
	};

	//get partial derivative of pressure by specific internal energy
	inline double GetPressureVolumeFractionDerivative(double ro, double alpha) {
		double a = a1 * a2 / (a2 * alpha + a1 * ( 1.0 - alpha ));
		a /= (a2 * alpha + a1 * ( 1.0 - alpha ));
		return a * ( (a1 - a2) * ro + (b1 - b2));
	};

	//get partial derivative of pressure by density
	inline double GetPressureDensityDerivative(double ro, double alpha) {
		return a1 * a2 / (a2 * alpha + a1 * ( 1.0 - alpha ));
	};

	//not usefull finctions
	//get internal energy
	virtual double GetInternalEnergy(double ro, double pressure) {
		return 0;
	};

	//get partial derivative of pressure by specific internal energy
	virtual double GetPressureEnergyDerivative(double ro, double e) {
		return 0;
	}
};

#endif
