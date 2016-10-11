#ifndef TurboStructured_Sensors_CellSensor
#define TurboStructured_Sensors_CellSensor

#include "Sensors/Sensor.h"


class CellSensor : public PointSensor {
private:
	sen_func getValue;		// target function

public:
	CellSensor(std::string _filename, ParallelManager& _parM, Grid& _grid, std::function<double(const std::valarray<double>&)> _getValue) : PointSensor(_filename, _parM, _grid) {
		getValue = _getValue;
	};

	virtual void Process(const std::valarray<double>& values) override {
		if (isActive == false) return;

		std::valarray<double> U = values[std::slice(idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
		double val = getValue(U);
		ofs.open(filename, std::ios_base::app);
		ofs << iteration << ' ' << val << ' ' << timer;
		std::endl(ofs);
		ofs.close();
	};

};


#endif