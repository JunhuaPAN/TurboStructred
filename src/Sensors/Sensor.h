#ifndef TurboStructured_Sensors_Sensor
#define TurboStructured_Sensors_Sensor

#include <string>
#include <valarray>
#include <functional>
#include <fstream>
#include "ParallelManager.h"

class Sensor
{
protected:
	bool _isActive;			// flag which show the sensor is active or not
	ParallelManager& _parM;	// parallel manager and grid
	Grid& _grid;

public:
	int iteration;			//number of iteration
	double timer;			//time of working
	std::string filename;	//file to write data

	std::ofstream ofs;

	Sensor(std::string _filename, ParallelManager& parM, Grid& grid) : filename(_filename), _parM(parM), _grid(grid) {
		_isActive = false;
		if (_parM.getRank() == 0) {
			ofs.open(filename);
			ofs.close();
		};
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

#endif