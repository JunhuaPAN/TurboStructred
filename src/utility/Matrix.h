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
	std::vector<std::vector<double>> element;

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

};

inline Matrix operator*(const Matrix& A, const Matrix& B)
{
	int N = A.Nrows;
	int M = B.Mcolumns;
	int K = A.Mcolumns;

	if(K != B.Nrows) {
		std::cout << "Error: Can't perform matrix production - incorrect sizes of multipliers\n";
		return Matrix();
	};
	
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

inline std::vector<double> operator*(const Matrix& A, const std::vector<double>& B)
{
	int N = A.Nrows;
	int M = A.Mcolumns;

	if(M != B.size()) {
		std::cout << "Error: Can't perform matrix/vector production - incorrect sizes of multipliers\n";
		return std::vector<double>();
	};
	
	std::vector<double> res(N, 0);
	for(int n = 0; n < N; n++) {
		for(int m = 0; m < M; m++) {
				res[n] += A.element[n][m]*B[m];
		};
	};

	return res;
};

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