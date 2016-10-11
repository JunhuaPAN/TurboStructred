#ifndef TurboStructured_Sensors_MValuePosSensor
#define TurboStructured_Sensors_MValuePosSensor

#include "Sensor.h"

class MValuePosXSensor : public Sensor {
protected:
	int i0;							// local index of first inner cell in slice
	int nVariables;					
	std::function<double(const std::valarray<double>&)> getValue;		// target function
public:
	MValuePosXSensor(std::string _filename, ParallelManager& _parM, Grid& _grid, std::function<double(const std::valarray<double>&)> _getValue) : Sensor(_filename, _parM, _grid) {
		getValue = _getValue;
	};

	inline void SetSensor(int j_slice, int k_slice, int _nVariables) {
		// check active or not
		if (	(j_slice >= grid.jMin - grid.dummyCellLayersY)
			&&	(j_slice <= grid.jMax - grid.dummyCellLayersY)
			&&	(k_slice >= grid.kMin - grid.dummyCellLayersZ)
			&&	(k_slice <= grid.kMax - grid.dummyCellLayersZ) )
		{
			isActive = true;
			i0 = (k_slice - grid.kMin + grid.dummyCellLayersZ) * grid.nlocalXAll * grid.nlocalYAll + (j_slice - grid.jMin + grid.dummyCellLayersY) * grid.nlocalXAll + grid.dummyCellLayersX;
			nVariables = _nVariables;
		};
	};

	virtual void Process(const std::valarray<double>& values) override {
		double uMax = -std::numeric_limits<double>::max();
		int iMax;		// local index of cell with maximum value
		if (isActive == true) {
			// Find local maximum value of target function
			int s = grid.nlocalX;
			for (auto i = 0; i < s; i++) {
				int idx = i0 + i;				// local cell index
				std::valarray<double> U = values[std::slice(idx * nVariables, nVariables, 1)];	//slice appropriate part of values array
				double u = getValue(U);
				if (u > uMax) {
					uMax = u;
					iMax = i;
				};
			};	// end of local maximum value search
		};
		parM.Barrier();
		double TotalMax = parM.Max(uMax);	// choose the biggest one

		// choose the appropriate process and write the result
		if (uMax == TotalMax) {
			ofs.open(filename, std::ios_base::app);
			ofs << iteration << ' ' << grid.CoordinateX[grid.iMin + iMax] << ' ' << timer;
			std::endl(ofs);
			ofs.close();
		};
	};
};


#endif