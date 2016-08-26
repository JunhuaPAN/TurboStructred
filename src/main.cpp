#include "Tests/tests.h"
#include <conio.h>



// Run the test
int main(int argc, char *argv[])
{
	SODxyzTest::RunExperiment(argc, argv, Direction::ZDirection);
	//RunSODTestRoe1D(argc, argv);
	//RunSODXTest(argc, argv);
	//BlasiusFlowTest::RunExperiment(argc, argv);
	//ContactDisTest::RunExperiment(argc, argv);

	//RunContactDisconTest1D(argc, argv);
	//ToroTests::RunExperiment(argc, argv);
	//AleshinExp::RunExperiment(argc, argv);
	//ExactEulerSolution::RunExperiment(argc, argv);
	
	system("pause");
	return 0;
};