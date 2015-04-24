#ifndef TurboStructured_Sensors_CellSensor
#define TurboStructured_Sensors_CellSensor

#include "Sensor.h"

class CellSensor : public Sensor {
private:
	int _idx;				//index of cell around sensor
public:
	CellSensor(std::string _filename, std::function<double(const std::valarray<double>&)> _getValue, int _nVariables) : Sensor(_filename, _getValue, _nVariables) {};

	inline void SetIndex(int idx) {
		_idx = idx;
	};

	virtual void Process(const std::valarray<double>& values) override {
		std::valarray<double>&& U = values[std::slice(_idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
		double val = getValue(U);
		ofs << iteration << ' ' << val << ' ' << timer;
		std::endl(ofs);
		ofs.flush();
	};
};


#endif