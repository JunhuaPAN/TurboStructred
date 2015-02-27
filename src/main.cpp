#include <iostream>
#include <vector>
#include "kernel\kernel.h";

std::vector<double> SODinitialDistribution(Vector r) {
	//Left

	return std::vector<double>();
};

int main(int argc, char *argv[])
{
	//main function
	KernelConfiguration conf;
	conf.nDims = 2;
	conf.nX = 100;
	conf.nY = 100;
	conf.LX = 1.0;
	conf.LY = 1.0;
	conf.isPeriodicX = true;
	conf.isPeriodicY = true;

	//init kernel
	Kernel kernel;
	kernel.Init(conf);
	
	//run computation
	kernel.Run();

	//kernel.SetInitialConditions();

	//finalize kernel
	kernel.Finilaze();

	return 0;
};