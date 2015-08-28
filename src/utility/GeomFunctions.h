#ifndef TurboStructured_Utility_GeomFunctions
#define TurboStructured_Utility_GeomFunctions

//Functions operating geometrical objects

#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Vector.h"
#include <Eigen/Dense>

using namespace Eigen;

//Compute gradient using least squares
Vector ComputeGradientByPoints(int ndims, Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values) {
  Vector grad;

  //Input	
  int n{ ndims }; //Number of dimensions (unknowns)	
  int m{ static_cast<int>(points.size()) }; //Number of equations
  int nrhs{ 1 }; //Number of right hand side
  MatrixXd A(n, m);
  VectorXd rhs(nrhs*m);

	//Validity check TO DO
	
	
	//Compose matrices 				
	for (int i = 0; i<points.size(); i++) {			
    Vector dr{ point - points[i] };
		A(i,0) = dr.x;
		if (n > 1) A(i,1) = dr.y;
		if (n > 2) A(i,2) = dr.z;
	};
	
	for (int i = 0; i<points.size(); i++) {						
    double dU{ value - values[i] };
		rhs(i) = dU;
	};	

	//Solve problem	
	JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);		
	Vector3d x = svd.solve(rhs); 

	//Check rank deficiency case
	if (svd.rank() < ndims) throw std::runtime_error{ "Could not solve for gradient due to rank deficiency"};

	//Output
	grad = Vector(x[0], 0, 0);
	if (ndims > 1) grad.y = x[1];
	if (ndims > 2) grad.z = x[2];

	return grad;
};


#endif
