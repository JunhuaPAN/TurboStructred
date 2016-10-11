#include "Tests/tests.h"
#include <conio.h>



// Run the test
int main(int argc, char *argv[])
{
	//PoiseuilleFlowTemp::RunExperiment(argc, argv);
	//RunPoiseuille1D(argc, argv);
	PoiseuilleChanel::RunExperiment(argc, argv);
	//Price2008KHI::RunExperiment(argc, argv);
	//AleshinExp::RunExperiment(argc, argv);
	//DrivenCavity::RunExperiment(argc, argv);
	//BlasiusFlow::RunExperiment(argc, argv);
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