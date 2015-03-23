#include "Tests\tests.h"


int main(int argc, char *argv[])
{
	//main function
	//RunSODTestRoe1D(argc, argv);
	//RunSODTestRoe2DX(argc, argv);
	//RunFluxesTest2D(argc, argv);
	//RunPoiseuille2DFVM(argc, argv);
	RunSODTestHybrid1DY(argc, argv);
	//RunDemchenkoTest2D(argc, argv);
	//RunSODTestHybrid2D(argc, argv);
	//RunDemchenkoTest3D(argc, argv);
	//RunSODTestHybrid1D(argc, argv);
	//RunPoiseuille3D(argc, argv);
	//RunShearFlow2D(argc, argv);

	return 0;
};