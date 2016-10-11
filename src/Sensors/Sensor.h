#ifndef TurboStructured_Sensors_Sensor
#define TurboStructured_Sensors_Sensor

#include <string>
#include <valarray>
#include <functional>
#include <fstream>
#include "ParallelManager.h"

// sensor function
using sen_func = std::function<double(const std::valarray<double>&)>;

class Sensor
{
protected:
	bool isActive;			// flag which show the sensor is active or not
	ParallelManager& parM;	// parallel manager and grid
	Grid& grid;

public:
	int iteration;			//number of iteration
	double timer;			//time of working
	std::string filename;	//file to write data
	std::ofstream ofs;		//stream to write data

	Sensor(std::string _filename, ParallelManager& _parM, Grid& _grid) : filename(_filename), parM(_parM), grid(_grid) {
		isActive = false;
		iteration = 0;
		timer = 0;
	};

	virtual void Process(const std::valarray<double>& values) = 0;

	//update iteration number and work time of sensor
	inline void NewRecord() {
		iteration++;
	};
	inline void UpdateIteration(int iter) {
		iteration = iter;
	};
	inline void UpdateTimer(double& time) {
		timer = time;
	};
};

class PointSensor : public Sensor {
protected:
	int idx;				// index of cell that contains this sensor
	int nVariables;
	Vector pos;

public:
	PointSensor(std::string _filename, ParallelManager& _parM, Grid& _grid) : Sensor(_filename, _parM, _grid) {
	}

	// Set position of the sensor
	void SetSensor(int i, int j, int k, int _nVariables) {
		// check active or not
		if ((i >= grid.iMin) && (i <= grid.iMax) && (j >= grid.jMin) && (j <= grid.jMax) && (k >= grid.kMin) && (k <= grid.kMax)) {
			isActive = true;
			idx = grid.getSerialIndexLocal(i, j, k);
			nVariables = _nVariables;
			
			// write coordinates
			pos = Vector(grid.CoordinateX[i], grid.CoordinateY[j], grid.CoordinateZ[k]);

			// open file and write censor positions
			ofs.open(filename);
			ofs << "TITLE = \"X=" << pos.x;
			ofs << ", Y=" << pos.y;
			ofs << ", Z=" << pos.z << '\"';
			std::endl(ofs);
			ofs.close();
		};
	};
	void SetSensor(CellIdx cell, int _nVariables) {
		// check active or not
		this->SetSensor(cell.i, cell.j, cell.k, _nVariables);
	};

	virtual void Process(const std::valarray<double>& values) = 0;
};

#endif