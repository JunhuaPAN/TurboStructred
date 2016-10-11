#ifndef TurboStructured_Sensors_MValuePosSensor2
#define TurboStructured_Sensors_MValuePosSensor2

#include "MValuePosXSensor.h"

// more accurate version of previous sensor //
// The position is computed by 2nd order polinom reconstruction for nonuniform grid //
class MValuePosXSensor2 : public MValuePosXSensor {
public:
	// constructor
	MValuePosXSensor2(std::string _filename, ParallelManager& _parM, Grid& _grid, std::function<double(const std::valarray<double>&)> _getValue) : MValuePosXSensor(_filename, _parM, _grid, _getValue) {};

	virtual void Process(const std::valarray<double>& values) override {
		double uMax = -std::numeric_limits<double>::max();
		int idxMax;			// local index of cell with maximum value
		if (isActive == true) {
			// Find local maximum value of target function
			int s = grid.nlocalX;
			for (auto i = 0; i < s; i++) {
				int idx = i0 + i;		// local index of cell
				std::valarray<double> U = values[std::slice(idx * nVariables, nVariables, 1)];	//	slice appropriate part of values array	
				double u = getValue(U);
				if (u > uMax) {
					uMax = u;
					idxMax = i;
				};
			};	// end of Find
		};
		parM.Barrier();
		double TotalMax = parM.Max(uMax);	// choose the biggest one

		// choose the appropriate process. Compute maximum position and write it in the file
		if (uMax == TotalMax) {
			double xmax;
			// left position case
			if (idxMax == 0) {
				xmax = grid.CoordinateX[grid.iMin] - 0.5 * grid.hx[grid.iMin];
			} // right position
			else if (idxMax == grid.nlocalX - 1) {
				xmax = grid.CoordinateX[grid.iMax] + 0.5 * grid.hx[grid.iMax];
			} // ordinary case
			else {		
				int i = grid.iMin + idxMax;			// i global
				double xi = grid.CoordinateX[i];		// position of central cell 
				double xip = grid.CoordinateX[i + 1];	// position of right cell
				double xim = grid.CoordinateX[i - 1];	// position of left cell
				double dxp = xip - xi;					// right cell scale
				double dxm = xi - xim;					// left one
				int idx = i0 + idxMax;					// local index of central cell

				// second order reconstruction
				// U = a (x - xi) (x - xi) + b (x - xi) + (u2 - u1)

				// find (a * dxm) and (b/a)
				std::valarray<double> U = values[std::slice((idx - 1) * nVariables, nVariables, 1)];
				double ul = getValue(U);
				U = values[std::slice((idx + 1) * nVariables, nVariables, 1)];
				double ur = getValue(U);
				double v2 = uMax - ul;		// v1 = ul - ul = 0
				double v3 = ur - ul;
				double adxm = v3 * dxm / ((dxm + dxp) * dxp) - v2 / dxp;	//  a * dxm
				double bda = dxm + v2 / adxm;								//  b / a
				xmax = xi - 0.5 * bda;
			}; // end of maximum position computing

			ofs.open(filename, std::ios_base::app);
			ofs << iteration << ' ' << xmax << ' ' << timer;
			std::endl(ofs);
			ofs.close();
		};
	};
};


#endif