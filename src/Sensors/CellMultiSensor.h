#ifndef TurboStructured_Sensors_CellMultiSensor
#define TurboStructured_Sensors_CellMultiSensor

#include "Sensors/Sensor.h"
#include <vector>

class CellMultiSensor : public PointSensor {
private:
	std::vector< sen_func > senFunctions;		// target function

public:
	CellMultiSensor(std::string _filename, ParallelManager& _parM, Grid& _grid, std::vector< sen_func >& _senFunctions) :
		PointSensor(_filename, _parM, _grid), senFunctions(_senFunctions) {};

	virtual void Process(const std::valarray<double>& values) override {
		if (isActive == false) return;

		std::valarray<double> U = values[std::slice(idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
		ofs.open(filename, std::ios_base::app);
		ofs << iteration << ' ';
		for (auto f : senFunctions) ofs << f(U) << ' ';
		ofs << timer;
		std::endl(ofs);
		ofs.close();
	};
};


#endif