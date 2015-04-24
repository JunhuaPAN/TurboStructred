#ifndef TurboStructured_Sensors_Sensor
#define TurboStructured_Sensors_Sensor

#include <string>
#include <valarray>
#include <functional>
#include <fstream>

class Sensor
{
protected:
	int nVariables;

public:
	int iteration;			//number of iteration
	double timer;			//time of working
	std::string filename;	//file to write data
	std::function<double(const std::valarray<double>&)> getValue;
	std::ofstream ofs;

	Sensor(std::string _filename, std::function<double(const std::valarray<double>&)> _getValue, int _nVariables) : filename(_filename), getValue(_getValue), nVariables(_nVariables) {
		//ofs.open(filename, std::ifstream::app);
		ofs.open(filename);
		iteration = 0;
		timer = 0;
	};

	~Sensor() {
		ofs.close();
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