#ifndef TurboStructured_kernel_grid
#define TurboStructured_kernel_grid

#include "KernelConfiguration.h"
#include <cassert>

// Info about cell
struct CellInfo {
	//int i, j, k;		// indexes of cell
	double hx, hy, hz;	// geometric sizes
	double x, y, z;		// center position
};

// Calculation kernel
class Grid {
public:
	int nDims;
	int iMin;		// grid par
	int iMax;		// grid par
	int jMin;		// grid par
	int jMax;		// grid par
	int kMin;		// grid par
	int kMax;		// grid par

	//Grid information
	int nX; //Number of cells in x dimension		// grid par
	int nY; //Number of cells in y dimension		// grid par
	int nZ; //Number of cells in z dimension		// grid par
	int nXAll; //Number of cells in x dimension including dummy layers		// grid par
	int nYAll; //Number of cells in y dimension including dummy layers		// grid par
	int nZAll; //Number of cells in z dimension including dummy layers		// grid par
	bool IsPeriodicX;
	bool IsPeriodicY;
	bool IsPeriodicZ;

	//compression coefficient of grid sizes towards the borders in X Y and Z directions (as geometric progression)
	double qx;
	double qy;
	double qz;

	//Same for local cells
	int nlocalX;		// grid par
	int nlocalY;		// grid par
	int nlocalZ;		// grid par
	int nlocalXAll;		// grid par
	int nlocalYAll;		// grid par
	int nlocalZAll;		// grid par
	int nCellsLocal;	// cells number
	int nCellsLocalAll;
	std::valarray<double> CoordinateX; //Cell center coordinates			// grid par
	std::valarray<double> CoordinateY; //Cell center coordinates			// grid par
	std::valarray<double> CoordinateZ; //Cell center coordinates			// grid par
	std::valarray<double> hx;
	std::valarray<double> hy;
	std::valarray<double> hz;

	// local grid volumes in order of local serial indexes of cells
	std::valarray<double> volumes;

	// number of dummy layers
	int dummyCellLayersX;
	int dummyCellLayersY;
	int dummyCellLayersZ;

	//Constructor
	Grid() {};

	// Get local index for cell
	inline int getSerialIndexLocal(int i, int j, int k) {
		int sI = (k - kMin + dummyCellLayersZ) * nlocalXAll * nlocalYAll + (j - jMin + dummyCellLayersY) * nlocalXAll + (i - iMin + dummyCellLayersX);
		return sI;
	};

	//Initialize Global Grid
	void InitGlobal(KernelConfiguration& config);

	// Initialize Local Grid
	void InitLocal(KernelConfiguration& config);

	// Create new grid as part of initial grid
	Grid CreateSubGrid(int iMin, int iMax, int jMin, int jMax, int kMin, int kMax);

	// Take information about the cell
	CellInfo CreateCell(int i, int j, int k);
};

//Initialize Global Grid
void Grid::InitGlobal(KernelConfiguration& config) {
	nDims = config.nDims;

	//Initialize local grid
	int dummyCellLayers = config.DummyLayerSize; //number of dummy cell layers				
	nX = config.nX;
	IsPeriodicX = config.isPeriodicX;
	dummyCellLayersX = dummyCellLayers;
	nY = 1;
	IsPeriodicY = true;
	dummyCellLayersY = 0;
	nZ = 1;
	IsPeriodicZ = true;
	dummyCellLayersZ = 0;
	if (nDims > 1) {
		nY = config.nY;
		IsPeriodicY = config.isPeriodicY;
		dummyCellLayersY = dummyCellLayers;
	};
	if (nDims > 2) {
		nZ = config.nZ;
		IsPeriodicZ = config.isPeriodicZ;
		dummyCellLayersZ = dummyCellLayers;
	};
	// degrees of freedom of global grid
	nXAll = nX + 2 * dummyCellLayersX;
	nYAll = nY + 2 * dummyCellLayersY;
	nZAll = nZ + 2 * dummyCellLayersZ;

	// Prepare coordinate arrays
	CoordinateX.resize(nXAll);
	CoordinateY.resize(nYAll);
	CoordinateZ.resize(nZAll);
	hx.resize(nXAll);
	hy.resize(nYAll);
	hz.resize(nZAll);

	//nullify compression if we have uniform grid
	if (config.isUniformAlongX == false) {
		qx = config.qx;
		assert(IsPeriodicX == false);
	}
	else qx = 1;
	if (config.isUniformAlongY == false) {
		qy = config.qy;
		assert(IsPeriodicY == false);
	}
	else qy = 1;
	if (config.isUniformAlongZ == false) {
		qz = config.qz;
		assert(IsPeriodicZ == false);
	}
	else qz = 1;

	return;
};

// Initialize Local Grid
void Grid::InitLocal(KernelConfiguration& config) {
	// Global domain sizes
	double Lx = config.LX;
	double Ly = config.LY;
	double Lz = config.LZ;
	if (nDims < 3) Lz = 0;
	if (nDims < 2) Ly = 0;

	// compute local sizes
	nlocalX = iMax - iMin + 1;
	nlocalY = jMax - jMin + 1;
	nlocalZ = kMax - kMin + 1;
	nlocalXAll = nlocalX + 2 * dummyCellLayersX;
	nlocalYAll = nlocalY + 2 * dummyCellLayersY;
	nlocalZAll = nlocalZ + 2 * dummyCellLayersZ;
	nCellsLocal = nlocalX * nlocalY * nlocalZ;
	nCellsLocalAll = nlocalXAll * nlocalYAll * nlocalZAll;

	// Compressible Grid scheme
	//
	//	_|________|_________|_
	//	 |        |         |
	//	 |        |         |
	//	 |        |         |
	//	 |        |         |					Second layer of Inner cells
	//	 |        |         |
	//	_|________|_________|_
	//	 |        |         |		)
	//	 |        |         |		> H			First layer of Inner cells
	//	_|________|_________|_		)	
	//	______________________		<-- Border is Here
	//	 |        |         |		)
	//	 |        |         |		> H			First lauer of Fictious cells
	//	_|________|_________|_		)
	//	 |        |         |		)
	//	 |        |         |		> H			Second layer of Fictious cells
	//	_|________|_________|_		)
	//	 |		  |			|

	// fill cell centers positions and edges sizes
	double h_x = Lx / nX;			//uniform grid case
	if ((qx != 1) && (nX % 2 == 0)) h_x = 0.5 * Lx * (1.0 - qx) / (1.0 - pow(qx, 0.5 * nX));	// X step around the border for even cells number
	if ((qx != 1) && (nX % 2 == 1)) h_x = Lx * (1.0 - qx) / (2.0 - pow(qx, 0.5 * (nX - 1)) * (1.0 + qx));		// for odd cell numbers
	double xl = -(dummyCellLayersX * h_x) + 0.5 * h_x;				//left cell (global) position
	double xr = Lx + dummyCellLayersX * h_x - 0.5 * h_x;			//right cell (global) position
	for (int i = 0; i < dummyCellLayersX; i++) {
		CoordinateX[i] = xl;
		CoordinateX[nXAll - 1 - i] = xr;
		hx[i] = h_x;
		hx[nXAll - 1 - i] = h_x;
		xl += h_x;
		xr -= h_x;
	};

	for (int i = dummyCellLayersX; i <= 0.5 * (nX - 1) + dummyCellLayersX; i++) {
		CoordinateX[i] = xl;
		CoordinateX[nXAll - 1 - i] = xr;
		hx[i] = h_x;
		hx[nXAll - 1 - i] = h_x;
		xl += 0.5 * h_x * (1.0 + qx);
		xr -= 0.5 * h_x * (1.0 + qx);
		h_x *= qx;
	};

	double h_y = Ly / nY;			//uniform grid case
	if ((qy != 1) && (nY % 2 == 0)) h_y = 0.5 * Ly * (1.0 - qy) / (1.0 - pow(qy, 0.5 * nY));	// Y step around the border for even cells number
	if ((qy != 1) && (nY % 2 == 1)) h_y = Ly * (1.0 - qy) / (2.0 - pow(qy, 0.5 * (nY - 1)) * (1.0 + qy));		// for odd cell numbers in Y direction
	double yl = -(dummyCellLayersY * h_y) + 0.5 * h_y;			//	left cell (global) position
	double yr = Ly + dummyCellLayersY * h_y - 0.5 * h_y;		//	right cell (global) position
	for (int i = 0; i < dummyCellLayersY; i++) {
		CoordinateY[i] = yl;
		CoordinateY[nYAll - 1 - i] = yr;
		hy[i] = h_y;
		hy[nYAll - 1 - i] = h_y;
		yl += h_y;
		yr -= h_y;
	};

	for (int j = dummyCellLayersY; j <= 0.5 * (nY - 1) + dummyCellLayersY; j++) {
		CoordinateY[j] = yl;
		CoordinateY[nYAll - 1 - j] = yr;
		hy[j] = h_y;
		hy[nYAll - 1 - j] = h_y;
		yl += 0.5 * h_y * (1.0 + qy);
		yr -= 0.5 * h_y * (1.0 + qy);
		h_y *= qy;
	};

	double h_z = Lz / nZ;			//uniform grid case
	if ((qz != 1) && (nZ % 2 == 0)) h_z = 0.5 * Lz * (1.0 - qz) / (1.0 - pow(qz, 0.5 * nZ));	// Z step around the border for even cells number
	if ((qz != 1) && (nZ % 2 == 1)) h_z = Lz * (1.0 - qz) / (2.0 - pow(qz, 0.5 * (nZ - 1)) * (1.0 + qz));		// for odd cell numbers in Z direction
	double zl = -(dummyCellLayersZ * h_z) + 0.5 * h_z;				//left cell (global) position
	double zr = Lz + dummyCellLayersZ * h_z - 0.5 * h_z;			//right cell (global) position
	for (int i = 0; i < dummyCellLayersZ; i++) {
		CoordinateZ[i] = zl;
		CoordinateZ[nZAll - 1 - i] = zr;
		hz[i] = h_z;
		hz[nZAll - 1 - i] = h_z;
		zl += h_z;
		zr -= h_z;
	};

	for (int k = dummyCellLayersZ; k <= 0.5 * (nZ - 1) + dummyCellLayersZ; k++) {
		CoordinateZ[k] = zl;
		CoordinateZ[nZAll - 1 - k] = zr;
		hz[k] = h_z;
		hz[nZAll - 1 - k] = h_z;
		zl += 0.5 * h_z * (1.0 + qz);
		zr -= 0.5 * h_z * (1.0 + qz);
		h_z *= qz;
	};

	if (nDims < 2) hy[0] = 1.0;
	if (nDims < 3) hz[0] = 1.0;

	// compute volumes of local cells
	volumes.resize(nCellsLocalAll);
	for (int i = iMin - dummyCellLayersX; i <= iMax + dummyCellLayersX; i++)
	{
		for (int j = jMin - dummyCellLayersY; j <= jMax + dummyCellLayersY; j++)
		{
			for (int k = kMin - dummyCellLayersZ; k <= kMax + dummyCellLayersZ; k++)
			{
				// get serial local index
				int idx = getSerialIndexLocal(i, j, k);
				//Compute cell volume
				volumes[idx] = hx[i] * hy[j] * hz[k];
			};
		};
	};
}

// Create new grid as part of initial grid
Grid Grid::CreateSubGrid(int _iMin, int _iMax, int _jMin, int _jMax, int _kMin, int _kMax) {
	Grid res;

	// define min and max indexes of subgrid in all directions 
	// X first
	if ((_iMin > iMax) || (_iMax < iMin)) {
		// no intersection case
		res.nlocalX = 0;
		res.iMin = iMin;
		res.iMax = iMin - 1;
	} else {
		// default (subgrid includes full grid g )
		res.iMin = iMin;
		res.iMax = iMax;
		if ((_iMin >= iMin) && (_iMin <= iMax)) res.iMin = _iMin;
		if ((_iMax >= iMin) && (_iMax <= iMax)) res.iMax = _iMax;
		res.nlocalX = res.iMax - res.iMin + 1;
	};

	// Y first
	if ((_jMin > jMax) || (_jMax < jMin)) {
		// no intersection case
		res.nlocalY = 0;
		res.jMin = jMin;
		res.jMax = jMin - 1;
	}
	else {
		// default (subgrid includes full grid g )
		res.jMin = jMin;
		res.jMax = jMax;
		if ((_jMin >= jMin) && (_jMin <= jMax)) res.jMin = _jMin;
		if ((_jMax >= jMin) && (_jMax <= jMax)) res.jMax = _jMax;
		res.nlocalY = res.jMax - res.jMin + 1;
	};
	
	// Z first
	if ((_kMin > kMax) || (_kMax < kMin)) {
		// no intersection case
		res.nlocalZ = 0;
		res.kMin = kMin;
		res.kMax = kMin - 1;
	}
	else {
		// default (subgrid includes full grid g )
		res.kMin = kMin;
		res.kMax = kMax;
		if ((_kMin >= kMin) && (_kMin <= kMax)) res.kMin = _kMin;
		if ((_kMax >= kMin) && (_kMax <= kMax)) res.kMax = _kMax;
		res.nlocalZ = res.kMax - res.kMin + 1;
	};
	
	// define coordinates and space steps
	res.CoordinateX = CoordinateX; //Cell center coordinates			// grid par
	res.CoordinateY = CoordinateY; //Cell center coordinates			// grid par
	res.CoordinateZ = CoordinateZ; //Cell center coordinates			// grid par
	res.hx = hx;
	res.hy = hy;
	res.hz = hz;

	return res;
};

// Take information about the cell
CellInfo Grid::CreateCell(int i, int j, int k) {
	CellInfo res;
	res.hx = hx[i];
	res.hy = hy[j];
	res.hz = hz[k];
	res.x = CoordinateX[i];
	res.y = CoordinateY[j];
	res.z = CoordinateZ[k];

	return res;
};

#endif