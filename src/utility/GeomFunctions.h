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
#include "mkl.h"
#include "mkl_lapack.h"

//Compute gradient using least squares
Vector ComputeGradientByPoints(int ndims, Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values) {
	Vector grad;		
		
	//Input
	//ndims = 3;
	int n = ndims; //Number of dimensions (unknowns)
	int nPoints = points.size(); 
	int m = nPoints; //Number of equations
	int nrhs = 1; //Number of right hand side
	std::vector<double> a(n*m, 0);
	std::vector<double> b(nrhs*m, 0);
	std::vector<double> work(1 , 0);
	int lda = std::max(1,m); //lda = n;
	int ldb = std::max(lda, n); //ldb = nrhs;
	std::vector<int> jpvt(n, 0);
	double rcond = 0.01;
	int lwork = -1;

	//Output
	int _info;
	int rank;				
				
	//matrix.setlength(points.size(), 3);
	for (int i = 0; i<points.size(); i++) {			
		Vector dr = point - points[i];
		/*matrix[i][0] = dr.x;
		matrix[i][1] = dr.y;
		matrix[i][2] = dr.z;*/
		a[i + 0*lda] = dr.x;
		if (ndims > 1) a[i + 1*lda] = dr.y;
		if (ndims > 2) a[i + 2*lda] = dr.z;
	};

	//rhs.setlength(points.size());		
	for (int i = 0; i<points.size(); i++) {						
		double dU = value - values[i];
		//rhs[i] = dU;
		b[0 * ldb + i] = dU;
	};

	//Workspace querry
	dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
	lwork = (int)work[0];
	work.resize(lwork, 0);

	//Solve problem
	dgelsy(&m, &n, &nrhs, &a[0], &lda, &b[0], &ldb, &jpvt[0], &rcond, &rank, &work[0], &lwork, &_info);
	if (_info != 0) throw std::exception("Could not solve for gradient");
	if (rank < ndims) throw std::exception("Could not solve for gradient due to rank deficiency");
	grad = Vector(b[0], 0, 0);
	if (ndims > 1) grad.y = b[1];
	if (ndims > 2) grad.z = b[2];
	//std::cout<<"grad = "<<grad.x<<" "<<grad.y<<" "<<grad.z<<"\n";

	/*x.setlength(3);
	alglib::rmatrixsolvels(matrix, points.size(), 3, rhs, 0.0, info, rep, x);
	if (info != 1) throw Exception("Could not solve for gradient");
	grad = Vector(x[0], x[1], x[2]);	*/
	//std::cout<<"grad_old = "<<grad.x<<" "<<grad.y<<" "<<grad.z<<"\n";

	//Output for check
/*	std::cout<<"Matrix\n";
	for (int i = 0; i<points.size(); i++) {			
		Vector dr = point - points[i];
		std::cout<<dr.x<<" "<<dr.y<<" "<<dr.z<<"\n";			
	};
	std::cout<<"RHS\n";		
	for (int i = 0; i<points.size(); i++) {						
		double dU = value - values[i];
		std::cout<<dU<<"\n";	
	};*/

	return grad;
};


#endif