#ifndef TurboStructured_utility_Matrix
#define TurboStructured_utility_Matrix

#include <vector>
#include <iostream>

//Implementation of matrix product operation
//Matrix type
struct Matrix {
	int Nrows;	//number of rows
	int Mcolumns;	//number of columns

	//Elements of matrix
	typedef std::vector<double> DenseMatrixRow;
	std::vector<DenseMatrixRow > element;

	//constructor by default
	Matrix(): Nrows(0), Mcolumns(0)	
	{
		element.resize(0);
	};

	//constructor
	Matrix(int _Nrows, int _Mcolumns) : Nrows(_Nrows),	Mcolumns(_Mcolumns) {
		element.resize(Nrows);
		for(int i = 0; i < Nrows; i++)
			for(int j = 0; j < Mcolumns; j++)	{
				element[i].push_back(0);
			};
	};

	//operator []
	DenseMatrixRow& operator[](int i) {
		return element[i];
	};

};

inline Matrix operator*(const Matrix& A, const Matrix& B)
{
	int N = A.Nrows;
	int M = B.Mcolumns;
	int K = A.Mcolumns;
	assert(K == B.Nrows);
	
	Matrix res(N, M);
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < M; m++) {
			for(int k = 0; k < K; k++) {
				res.element[n][m] += A.element[n][k]*B.element[k][m];
			};
		};
	};

	return res;
};

inline std::vector<double> operator*(const Matrix& A, const Matrix::DenseMatrixRow& B)
{
	int N = A.Nrows;
	int M = A.Mcolumns;
	assert(M == B.size());
	
	std::vector<double> res(N, 0);
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < M; m++) {
				res[n] += A.element[n][m]*B[m];
		};
	};

	return res;
};

inline std::vector<double> operator*(const Matrix::DenseMatrixRow& x, const Matrix& A)
{
	int N = A.Nrows;
	int M = A.Mcolumns;
	assert(N == x.size());
	
	std::vector<double> res(M, 0);
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < M; m++) {
				res[m] += A.element[n][m]*x[m];
		};
	};

	return res;
};

inline double operator*(const Matrix::DenseMatrixRow& a, const Vector& b)
{	
	double res = 0;
	assert(a.size() == 3);
	res += a[0] * b.x;
	res += a[1] * b.y;
	res += a[2] * b.z;
	return res;
};

inline std::vector<double> operator*(const Matrix& A, const Vector& x)
{
	int N = A.Nrows;
	int M = A.Mcolumns; // == 3
	assert(M == 3);
	
	std::vector<double> res(N, 0);
	for(int n = 0; n < N; n++) {
		res[n] += A.element[n][0]*x.x;
		res[n] += A.element[n][1]*x.y;
		res[n] += A.element[n][2]*x.z;
	};

	return res;
};

//inline std::vector<double> operator*(const Matrix& A, const std::vector<double>& B)
//{
//	int N = A.Nrows;
//	int M = A.Mcolumns;
//
//	if(M != B.size()) {
//		std::cout << "Error: Can't perform matrix/vector production - incorrect sizes of multipliers\n";
//		return std::vector<double>();
//	};
//	
//	std::vector<double> res(N, 0);
//	for(int n = 0; n < N; n++) {
//		for(int m = 0; m < M; m++) {
//				res[n] += A.element[n][m]*B[m];
//		};
//	};
//
//	return res;
//};

inline std::vector<double> operator*(const Matrix& A, const double *B)
{
	int N = A.Nrows;
	int M = A.Mcolumns;

	//Unsafe 
	
	std::vector<double> res(N, 0);
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < M; m++) {
				res[n] += A.element[n][m]*B[m];
		};
	};

	return res;
};

#endif