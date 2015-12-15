#ifndef TurboStructured_Utility_Slices
#define TurboStructured_Utility_Slices

#include <cmath>
#include <iostream>
#include <vector>

class Slice {
public:
	int nDims;
	int i;
	int j;
	int k;

	Slice(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {
		nDims = 0;
		if (i == -1) ++nDims;
		if (j == -1) ++nDims;
		if (k == -1) ++nDims;
	};
};

class SliceInfo {
	std::vector<Slice> slices;


};

#endif
