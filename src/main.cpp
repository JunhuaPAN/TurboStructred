#include "Tests/tests.h"
#include <conio.h>



// Run the test
int main(int argc, char *argv[])
{
	RunShearFlow3DZ(argc, argv);
	//RunPoiseuille3D(argc, argv);
	//RunShearFlow2D(argc, argv);
	//RunPoiseuille2D(argc, argv);
	//RunShearFlow3D(argc, argv);
	
	_getch();
	return 0;
};