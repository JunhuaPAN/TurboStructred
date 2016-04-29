#ifndef TurboStructured_GasProperties
#define TurboStructured_GasProperties


// Collect all properties of the medium in that structure
struct GasProperties {
public:
	double gamma;
	double molarMass;
	double thermalConductivity;
	double viscosity;
	double universalGasConstant;
};



#endif
