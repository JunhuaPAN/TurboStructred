#ifndef TurboStructured_Sensors_MValuePosSensor
#define TurboStructured_Sensors_MValuePosSensor

#include "KernelConfiguration.h"
#include "Sensor.h"
#include <vector>

class CellSensor : public Sensor {
private:
	int _idy;					// indexes of cell's layer cell sensor
	int _idz;
	int _Nx, _Ny, _Nz;				// cell numbers
	std::vector<int> cell_list;	// list of cell indexes
public:
	CellSensor(std::string _filename, std::function<double(const std::valarray<double>&)> _getValue, int _nVariables) : Sensor(_filename, _getValue, _nVariables) {};

	inline void SetSensor(int idy, int idz, KernelConfiguration &conf) {
		_idy = idy;
		_idz = idz;
		_Nx = conf.nX;
		_Ny = conf.nY;
		_Nz = conf.nZ;

		cell_list.resize(_Nx);
		for (auto i = 0; i < _Nx; i++) {
			cell_id = 
		};
	};

	virtual void Process(const std::valarray<double>& values) override {
		std::valarray<double> U = values[std::slice(_idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
		double val = getValue(U);
		ofs << iteration << ' ' << val << ' ' << timer;
		std::endl(ofs);
		ofs.flush();
	};
};


#endif