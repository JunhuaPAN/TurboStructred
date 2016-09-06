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
	//	 |        |         |		> H			First layer of Fictious cells
	//	_|________|_________|_		)
	//	 |        |         |		)
	//	 |        |         |		> H			Second layer of Fictious cells
	//	_|________|_________|_		)
	//	 |		  |			|


	// check if we have uniform on X grid
	config.CompressionX.push_back(BlockNode(Lx));	// add final (empty) node (right border point)
	auto x_stp{ 1.0 };
	if (config.CompressionX[0].N_cells < 1) {
		x_stp = Lx / nX;
		config.CompressionX[0].N_cells = config.nX;		// cells number
	}	else {
		auto first_n = config.CompressionX[0];
		auto bl_len = config.CompressionX[1].pos - first_n.pos;
		auto Nc = first_n.N_cells;		// number of cells in first block

		// compute first space step
		if (first_n.q_com == 1.0) x_stp = bl_len / Nc;
		else x_stp = bl_len * (1.0 - first_n.q_com) / (1.0 - pow(first_n.q_com, Nc));		// first term of geometric progression
	};

	// compute node positions on the left side of the domain
	auto xc = -(dummyCellLayersX * x_stp) + 0.5 * x_stp;			//left cell (global) position
	auto c_idx = 0;
	for (; c_idx < dummyCellLayersX; c_idx++) {
		CoordinateX[c_idx] = xc;
		hx[c_idx] = x_stp;
		xc += x_stp;
	};

	// obtain space between block nodes
	for (auto b = 1; b < config.CompressionX.size(); b++) {
		auto blNode = config.CompressionX[b - 1];	// left node
		for (auto i = 0; i < blNode.N_cells; i++) {
			CoordinateX[c_idx + i] = xc;
			hx[c_idx + i] = x_stp;
			xc += 0.5 * x_stp * (1.0 + blNode.q_com);
			x_stp *= blNode.q_com;
		};
		c_idx += blNode.N_cells;
		// check if we have new not empty node
		if (config.CompressionX[b].N_cells > 0) {
			// compute new x_stp
			auto new_bl_len = config.CompressionX[b + 1].pos - config.CompressionX[b].pos;
			auto new_Nc = config.CompressionX[b].N_cells;
			auto new_q = config.CompressionX[b].q_com;

			// compute first space step for new block and first cell position
			if (new_q == 1.0) x_stp = new_bl_len / new_Nc;
			else x_stp = new_bl_len * (1.0 - new_q) / (1.0 - pow(new_q, new_Nc));
			xc = config.CompressionX[b].pos + 0.5 * x_stp;
		}
	};

	// final layers of dummy cells
	for (auto i = 0; i < dummyCellLayersX; i++) {
		CoordinateX[c_idx + i] = xc;
		hx[c_idx + i] = x_stp;
		xc += x_stp;
	};

	// Repeat that part for Y and Z 
	// Copy and change

	// for Y direction
	if (nDims > 1) {
		config.CompressionY.push_back(BlockNode(Ly));	// add final (empty) node (right border point)
		auto y_stp{ 1.0 };
		if (config.CompressionY[0].N_cells < 1) {
			y_stp = Ly / nY;
			config.CompressionY[0].N_cells = config.nY;		// cells number
		}
		else {
			auto first_n = config.CompressionY[0];
			auto bl_len = config.CompressionY[1].pos - first_n.pos;
			auto Nc = first_n.N_cells;		// number of cells in first block

			// compute first space step
			if (first_n.q_com == 1.0) y_stp = bl_len / Nc;
			else y_stp = bl_len * (1.0 - first_n.q_com) / (1.0 - pow(first_n.q_com, Nc));		// first term of geometric progression
		};

		// compute node positions on the left side of the domain
		auto yc = -(dummyCellLayersY * y_stp) + 0.5 * y_stp;			//left cell (global) position
		auto c_idx = 0;
		for (; c_idx < dummyCellLayersY; c_idx++) {
			CoordinateY[c_idx] = yc;
			hy[c_idx] = y_stp;
			yc += y_stp;
		};

		// obtain space between block nodes
		for (auto b = 1; b < config.CompressionY.size(); b++) {
			auto blNode = config.CompressionY[b - 1];	// left node
			for (auto i = 0; i < blNode.N_cells; i++) {
				CoordinateY[c_idx + i] = yc;
				hy[c_idx + i] = y_stp;
				yc += 0.5 * y_stp * (1.0 + blNode.q_com);
				y_stp *= blNode.q_com;
			};
			c_idx += blNode.N_cells;
			// check if we have new not empty node
			if (config.CompressionY[b].N_cells > 0) {
				// compute new y_stp
				auto new_bl_len = config.CompressionY[b + 1].pos - config.CompressionY[b].pos;
				auto new_Nc = config.CompressionY[b].N_cells;
				auto new_q = config.CompressionY[b].q_com;

				// compute first space step for new block
				if (new_q == 1.0) y_stp = new_bl_len / new_Nc;
				else y_stp = new_bl_len * (1.0 - new_q) / (1.0 - pow(new_q, new_Nc));
				yc = config.CompressionY[b].pos + 0.5 * y_stp;
			}
		};

		// final layers of dummy cells
		for (auto i = 0; i < dummyCellLayersY; i++) {
			CoordinateY[c_idx + i] = yc;
			hy[c_idx + i] = y_stp;
			yc += y_stp;
		};

	}

	// for Z direction
	if (nDims > 2) {
		config.CompressionZ.push_back(BlockNode(Lz));	// add final (empty) node (right border point)
		auto z_stp{ 1.0 };
		if (config.CompressionZ[0].N_cells < 1) {
			z_stp = Lz / nZ;
			config.CompressionZ[0].N_cells = config.nZ;		// cells number
		}
		else {
			auto first_n = config.CompressionZ[0];
			auto bl_len = config.CompressionZ[1].pos - first_n.pos;
			auto Nc = first_n.N_cells;		// number of cells in first block

											// compute first space step
			if (first_n.q_com == 1.0) z_stp = bl_len / Nc;
			else z_stp = bl_len * (1.0 - first_n.q_com) / (1.0 - pow(first_n.q_com, Nc));		// first term of geometric progression
		};

		// compute node positions on the left side of the domain
		auto zc = -(dummyCellLayersZ * z_stp) + 0.5 * z_stp;			//left cell (global) position
		auto c_idx = 0;
		for (; c_idx < dummyCellLayersZ; c_idx++) {
			CoordinateZ[c_idx] = zc;
			hz[c_idx] = z_stp;
			zc += z_stp;
		};

		// obtain space between block nodes
		for (auto b = 1; b < config.CompressionZ.size(); b++) {
			auto blNode = config.CompressionZ[b - 1];	// left node
			for (auto i = 0; i < blNode.N_cells; i++) {
				CoordinateZ[c_idx + i] = zc;
				hz[c_idx + i] = z_stp;
				zc += 0.5 * z_stp * (1.0 + blNode.q_com);
				z_stp *= blNode.q_com;
			};
			c_idx += blNode.N_cells;
			// check if we have new not empty node
			if (config.CompressionZ[b].N_cells > 0) {
				// compute new z_stp
				auto new_bl_len = config.CompressionZ[b + 1].pos - config.CompressionZ[b].pos;
				auto new_Nc = config.CompressionZ[b].N_cells;
				auto new_q = config.CompressionZ[b].q_com;

				// compute first space step for new block
				if (new_q == 1.0) z_stp = new_bl_len / new_Nc;
				else z_stp = new_bl_len * (1.0 - new_q) / (1.0 - pow(new_q, new_Nc));
				zc = config.CompressionZ[b].pos + 0.5 * z_stp;
			}
		};

		// final layers of dummy cells
		for (auto i = 0; i < dummyCellLayersZ; i++) {
			CoordinateZ[c_idx + i] = zc;
			hz[c_idx + i] = z_stp;
			zc += z_stp;
		};

	}

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