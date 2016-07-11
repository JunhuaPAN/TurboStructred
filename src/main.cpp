#include "Tests/tests.h"
#include <conio.h>



// Run the test
int main(int argc, char *argv[])
{
	RunSODTestRoe1D(argc, argv);
	//ToroTests::RunExperiment(argc, argv);
	//AleshinExp::RunExperiment(argc, argv);
	//BlasiusFlowTest::RunExperiment(argc, argv);
	//ExactEulerSolution::RunExperiment(argc, argv);
	
	_getch();
	return 0;
};