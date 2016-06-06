#ifndef TurboStructured_Utility_Stencil
#define TurboStructured_Utility_Stencil


#include <vector>
#include "Vector.h"
#include "grid.h"

struct CellIdx {
	int i, j, k;

	CellIdx(int _i, int _j, int _k) : i(_i), j(_j), k(_k) {};
};


class Stencil {
private:
	const int nDims;

public:
	const Grid* grid;

	// constructors
	Stencil() : nDims(1) {};
	Stencil(int _nDims, Grid& _grid) : nDims(_nDims), grid(&_grid) {};

	// create some type stencils
	auto InternalBorderStencil(Direction dir, CellIdx& cell) {
		std::vector<CellIdx> stn;
		stn.push_back(cell);	// first add the input cell

		// add other cells to create stencil 
		switch (dir) {
		case Direction::XDirection:
			if (nDims > 1) {
				stn.push_back(CellIdx(cell.i, cell.j + 1, cell.k));
				stn.push_back(CellIdx(cell.i, cell.j - 1, cell.k));
			};
			if (nDims > 2) {
				stn.push_back(CellIdx(cell.i, cell.j, cell.k + 1));
				stn.push_back(CellIdx(cell.i, cell.j, cell.k - 1));
			};
			break;
		case Direction::YDirection:
			stn.push_back(CellIdx(cell.i + 1, cell.j, cell.k));
			stn.push_back(CellIdx(cell.i - 1, cell.j, cell.k));
			if (nDims > 2) {
				stn.push_back(CellIdx(cell.i, cell.j, cell.k + 1));
				stn.push_back(CellIdx(cell.i, cell.j, cell.k - 1));
			};
			break;
		case Direction::ZDirection:
			stn.push_back(CellIdx(cell.i + 1, cell.j, cell.k));
			stn.push_back(CellIdx(cell.i - 1, cell.j, cell.k));
			stn.push_back(CellIdx(cell.i, cell.j + 1, cell.k));
			stn.push_back(CellIdx(cell.i, cell.j - 1, cell.k));
			break;
		default:
			break;
		};

		return stn;
	}


};


#endif
