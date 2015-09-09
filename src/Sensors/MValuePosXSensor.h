#ifndef TurboStructured_Sensors_MValuePosSensor
#define TurboStructured_Sensors_MValuePosSensor

#include "Sensor.h"

class MValuePosXSensor : public Sensor {
protected:
	int i0;							//first index of cell sequence
	int nVariables;					
	std::function<double(const std::valarray<double>&)> getValue;		// target function
public:
	MValuePosXSensor(std::string _filename, ParallelManager& _parM, Grid& _grid, std::function<double(const std::valarray<double>&)> _getValue) : Sensor(_filename, _parM, _grid) {
		getValue = _getValue;
	};

	inline void SetSensor(int j_slice, int k_slice, int _nVariables) {
		// check active or not
		if ((j_slice >= _grid.jMin) && (j_slice <= _grid.jMax) && (k_slice >= _grid.kMin) && (k_slice <= _grid.kMax)) {
			_isActive = true;
			i0 = k_slice * _grid.nXAll * _grid.nYAll + j_slice * _grid.nXAll + _grid.iMin;	// TO DO FIX
			nVariables = _nVariables;
		};
	};

	virtual void Process(const std::valarray<double>& values) override {
		double vMax = std::numeric_limits<double>::min();
		int iMax;
		if (_isActive == true) {
			// Find local maximum value of target function
			int s = _grid.nlocalX;
			for (auto i = 0; i < s; i++) {
				int idx = i0 + i;		// local index of cell
				std::valarray<double> U = values[std::slice(idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
				double v = getValue(U);
				if (v > vMax) {
					vMax = v;
					iMax = i;
				};
			};	// end of Find
		};
		double TotalMax = _parM.Max(vMax);	// choose the biggest one

		// choose the appropriate process and write the result
		if (vMax == TotalMax) {
			ofs.open(filename, std::ios_base::app);
			ofs << iteration << ' ' << _grid.coordsX[_grid.iMin + iMax] << ' ' << timer;
			std::endl(ofs);
			ofs.close();
		};
	};
};


#endif