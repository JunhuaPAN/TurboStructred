#include "Tests/tests.h"
#include <conio.h>



// Run the test
int main(int argc, char *argv[])
{
	BlasiusFlowTest::RunExperiment(argc, argv);
	//BlasiusFlowTest_SLtest::RunExperiment(argc, argv);		// ! that test detects a bug while use 8 cores
	//ExactEulerSolution::RunExperiment(argc, argv);
	//SODxyzTest::RunExperiment(argc, argv, Direction::ZDirection);
	//RunSODTestRoe1D(argc, argv);
	//RunSODXTest(argc, argv);

	//ContactDisTest::RunExperiment(argc, argv);

	//RunContactDisconTest1D(argc, argv);
	//ToroTests::RunExperiment(argc, argv);
	//AleshinExp::RunExperiment(argc, argv);
	
	system("pause");
	return 0;
};