#ifndef TurboStructured_Sensors_CellSensor
#define TurboStructured_Sensors_CellSensor

#include "Sensors/Sensor.h"


class CellSensor : public Sensor {
private:
	int _idx;				// index of cell around sensor
	std::function<double(const std::valarray<double>&)> getValue;		// target function
	int nVariables;

public:
	CellSensor(std::string _filename, ParallelManager& _parM, Grid& _grid, std::function<double(const std::valarray<double>&)> _getValue) : Sensor(_filename, _parM, _grid) {
		getValue = _getValue;
	};

	inline void SetSensor(int i, int j, int k, int _nVariables) {
		// check active or not
		if ((i >= _grid.iMin) && (i <= _grid.iMax) && (j >= _grid.jMin) && (j <= _grid.jMax) && (k >= _grid.kMin) && (k <= _grid.kMax)) {
			_isActive = true;
			_idx = (k * _grid.nlocalXAll * _grid.nlocalYAll + j * _grid.nlocalXAll + i);
			nVariables = _nVariables;
		};
	};

	virtual void Process(const std::valarray<double>& values) override {
		if (_isActive == false) return;

		std::valarray<double> U = values[std::slice(_idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
		double val = getValue(U);
		ofs.open(filename, std::ios_base::app);
		ofs << iteration << ' ' << val << ' ' << timer;
		std::endl(ofs);
		ofs.close();
	};

};


#endif