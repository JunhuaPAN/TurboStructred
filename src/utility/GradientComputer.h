#ifndef TurboStructured_Utility_GradientComputer
#define TurboStructured_Utility_GradientComputer


#include <vector>
#include "Vector.h"
#include <Eigen\Dense>

class GradientComputer {
private:
	const int nDims;
public:

	// constructors
	GradientComputer() : nDims(1) {};
	GradientComputer(int _nDims) : nDims(_nDims) {};

	// compute gradient using least squares
	Vector ExecuteByLeastSquares(Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values) {
		Vector grad;

		//Input	
		int m = points.size(); //Number of equations
		int nrhs = 1; //Number of right hand side
		Eigen::MatrixXd A(m, nDims);
		Eigen::VectorXd rhs(nrhs*m);

		//Compose matrices 				
		for (int i = 0; i<points.size(); i++) {
			Vector dr = point - points[i];
			A(i, 0) = dr.x;
			if (nDims > 1) A(i, 1) = dr.y;
			if (nDims > 2) A(i, 2) = dr.z;
		};

		for (int i = 0; i<points.size(); i++) {
			double dU = value - values[i];
			rhs(i) = dU;
		};

		//Solve problem	
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::VectorXd x = svd.solve(rhs);

		//Check rank deficiency case
		if (svd.rank() < nDims) throw std::exception("Could not solve for gradient due to rank deficiency");

		//Output
		grad = Vector(x(0), 0, 0);
		if (nDims > 1) grad.y = x(1);
		if (nDims > 2) grad.z = x(2);

		return grad;
	};

};
 

#endif