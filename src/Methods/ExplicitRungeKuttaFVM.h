#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "KernelConfiguration.h"
#include "utility/Vector.h"
#include "utility/Timer.h"
#include "utility/GradientComputer.h"
#include "utility/Stencil.h"
#include "RiemannSolvers/RiemannSolversList.h"
#include "kernel.h"
#include "Reconstruction/ReconstructorsList.h"


//Base class for all solution methods that desribe iterations process in detail
template <typename ReconstructionType>
class ExplicitRungeKuttaFVM : public Kernel {

	//Current riemann problem solver
	std::unique_ptr<RiemannSolver> _riemannSolver;

	//Method parameters
	int RungeKuttaOrder;
	double CFL;

public:
	//Inherit constructor
	//using Kernel::Kernel;
	ExplicitRungeKuttaFVM(int* argc, char **argv[]) : Kernel(argc, argv) {};

	//Additional data
	std::valarray<double> spectralRadius;		//array for storing spectral radiuses
	std::vector<std::vector<double>> fluxes; // array for storing fluxes
	std::vector<Vector> gradientVelocityX;  //array for storing gradients of u
	std::vector<Vector> gradientVelocityY;  //array for storing gradients of v
	std::vector<Vector> gradientVelocityZ;  //array for storing gradients of w
	std::vector< std::vector< std::vector<ReconstructionType> > > reconstructions; //array for storing 

	// Initialize Method
	virtual void InitializeMethod(KernelConfiguration& kernelConfig) override {
		// Method specific part
		MethodConfiguration config = kernelConfig.methodConfiguration;
		assert(config.RungeKuttaOrder != 0);							// Runge-Kutta order is not initialized
		assert(config.RiemannProblemSolver != RPSolver::NoSolver);		// No solver was choosen

		// Set parameters CFL and so on
		CFL = config.CFL;
		RungeKuttaOrder = config.RungeKuttaOrder;

		// Set RP solver
		if (config.RiemannProblemSolver == RPSolver::RoePikeSolver) {
			_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.Gamma, config.Eps));
		};
		if (config.RiemannProblemSolver == RPSolver::GodunovSolver) {
			_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new GodunovSolverPerfectGasEOS(kernelConfig.Gamma, config.Eps));
		};

		// Allocate memory for all structures
		spectralRadius.resize(grid.nCellsLocalAll);

		// Resize reconstructions 3D array and initialize them
		reconstructions.resize(grid.nlocalX + 2 * min(1, grid.dummyCellLayersX));
		for (auto i = 0; i < reconstructions.size(); i++) {
			size_t sy = grid.nlocalY + 2 * min(1, grid.dummyCellLayersY);
			reconstructions[i].resize(sy);

			for(auto j = 0; j < sy; j++) {
				size_t sz = grid.nlocalZ + 2 * min(1, grid.dummyCellLayersZ);
				reconstructions[i][j].resize(sz);

				for (auto k = 0; k < sz; k++) reconstructions[i][j][k].Init(nVariables, nDims);
			};
		};

		// end of method initialization part
	};


	//=================================   Viscous Part  ======================================

	//Compute all required gradient in each inner cell
	void ComputeVelocityGradients() {
		//Gradients required
		auto getu = [](double *U) { return U[1] / U[0]; };
		auto getv = [](double *U) { return U[2] / U[0]; };
		auto getw = [](double *U) { return U[3] / U[0]; };
		if (gradientVelocityX.size() != grid.nCellsLocalAll) gradientVelocityX.resize(grid.nCellsLocalAll);
		if (gradientVelocityY.size() != grid.nCellsLocalAll) gradientVelocityY.resize(grid.nCellsLocalAll);
		if (gradientVelocityZ.size() != grid.nCellsLocalAll) gradientVelocityZ.resize(grid.nCellsLocalAll);

		//Inside domain
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					int idx = grid.getSerialIndexLocal(i, j, k);

					//Cell values
					double *V = getCellValues(i, j, k);

					//Gradients required
					gradientVelocityX[idx] = ComputeGradientNonUniformStructured(i, j, k, getu);
					gradientVelocityY[idx] = ComputeGradientNonUniformStructured(i, j, k, getv);
					gradientVelocityZ[idx] = ComputeGradientNonUniformStructured(i, j, k, getw);
				};
			};
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Velocity gradients have'n calculated inside the domain\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Interprocessor exchange
		ExchangeGradients(gradientVelocityX);
		if (nDims > 1) ExchangeGradients(gradientVelocityY);
		if (nDims > 2) ExchangeGradients(gradientVelocityZ);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Gradients exchange executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Boundaries
		ComputeDummyCellGradients();

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Gradients calculated in dummy cells\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Exchange gradient values between processors
	void ExchangeGradients(std::vector<Vector>& cellGradients) {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//Get parallel run info
		MPI_Comm comm = pManager->getComm();

		//Buffers
		std::vector<double> bufferToRecv;
		std::vector<double> bufferToSend;

		//Allocate buffers
		int layerSizeX = grid.nlocalYAll * grid.nlocalZAll;
		int layerSizeY = grid.nlocalXAll * grid.nlocalZAll;
		int layerSizeZ = grid.nlocalXAll * grid.nlocalYAll;
		int bufferSize = std::max(std::max(layerSizeX, layerSizeY), layerSizeZ) * nDims;// * nVariables;
		bufferToSend.resize(bufferSize);
		bufferToRecv.resize(bufferSize);

		// TO DO lift from lyabmda to member function
		auto exchangeLayers = [&](Direction direction, int iSend, int rankDest, int iRecv, int rankSource) {
			int nRecv = 0;
			int nSend = 0;

			//Fill buffer with gradient values
			int idxBuffer = 0;
			for (i = grid.iMin - grid.dummyCellLayersX; i <= grid.iMax + grid.dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iSend)) continue; //skip
				for (j = grid.jMin - grid.dummyCellLayersY; j <= grid.jMax + grid.dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iSend)) continue; //skip
					for (k = grid.kMin - grid.dummyCellLayersZ; k <= grid.kMax + grid.dummyCellLayersZ; k++) {
						if ((direction == Direction::ZDirection) && (k != iSend)) continue; //skip
						nRecv += nDims;
						if (rankDest == -1) continue;

						int idxCell = grid.getSerialIndexLocal(i, j, k);
						bufferToSend[idxBuffer * nDims] = cellGradients[idxCell].x;
						if (nDims > 1) bufferToSend[idxBuffer * nDims + 1] = cellGradients[idxCell].y;
						if (nDims > 2) bufferToSend[idxBuffer * nDims + 2] = cellGradients[idxCell].z;
						idxBuffer++;
					};
				};
			};

			if (rankDest != -1) nSend = idxBuffer * nDims; //Determine number of doubles to send

														   //Determine recive number
			if (rankSource == -1) {
				nRecv = 0;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankDest, rankSource, &bufferToSend.front(), nSend, &bufferToRecv.front(), nRecv);

			//Write recieved values back
			idxBuffer = 0;
			for (i = grid.iMin - grid.dummyCellLayersX; i <= grid.iMax + grid.dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iRecv)) continue; //skip
				for (j = grid.jMin - grid.dummyCellLayersY; j <= grid.jMax + grid.dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iRecv)) continue; //skip
					for (k = grid.kMin - grid.dummyCellLayersZ; k <= grid.kMax + grid.dummyCellLayersZ; k++) {
						if ((direction == Direction::ZDirection) && (k != iRecv)) continue; //skip
						int idxCell = grid.getSerialIndexLocal(i, j, k);
						cellGradients[idxCell].x = bufferToRecv[idxBuffer * nDims];
						cellGradients[idxCell].y = 0;
						if (nDims > 1) cellGradients[idxCell].y = bufferToRecv[idxBuffer * nDims + 1];
						cellGradients[idxCell].z = 0;
						if (nDims > 2) cellGradients[idxCell].z = bufferToRecv[idxBuffer * nDims + 2];
						idxBuffer++;
					};
				};
			};

			//Syncronize
			pManager->Barrier();
		};

		//Determine neighbours' ranks
		int rankL = pManager->GetRankByCartesianIndexShift(-1, 0, 0);
		int rankR = pManager->GetRankByCartesianIndexShift(+1, 0, 0);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		for (int layer = 1; layer <= grid.dummyCellLayersX; layer++) {
			// Minus direction exchange
			int iSend = grid.iMin + layer - 1; // layer index to send
			int iRecv = grid.iMax + layer; // layer index to recv
			exchangeLayers(Direction::XDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = grid.iMax - layer + 1; // layer index to send
			iRecv = grid.iMin - layer; // layer index to recv
			exchangeLayers(Direction::XDirection, iSend, rankR, iRecv, rankL);
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange gradients in X-direction was executed\n";
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction exchange

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);

		//Debug info
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		//Y direction exchange
		for (int layer = 1; layer <= grid.dummyCellLayersY; layer++) {
			// Minus direction exchange
			int iSend = grid.jMin + layer - 1; // layer index to send
			int iRecv = grid.jMax + layer; // layer index to recv
			exchangeLayers(Direction::YDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = grid.jMax - layer + 1; // layer index to send
			iRecv = grid.jMin - layer; // layer index to recv
			exchangeLayers(Direction::YDirection, iSend, rankR, iRecv, rankL);
		};

		//Debug info
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange gradients in Y-direction was executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction exchange

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, 0, -1);
		rankR = pManager->GetRankByCartesianIndexShift(0, 0, +1);

		//Debug info
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		//Z direction exchange 
		for (int layer = 1; layer <= grid.dummyCellLayersZ; layer++) {
			// Minus direction exchange
			int iSend = grid.kMin + layer - 1; // layer index to send
			int iRecv = grid.kMax + layer; // layer index to recv
			exchangeLayers(Direction::ZDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = grid.kMax - layer + 1; // layer index to send
			iRecv = grid.kMin - layer; // layer index to recv
			exchangeLayers(Direction::ZDirection, iSend, rankR, iRecv, rankL);
		};

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange gradients in Z-direction was executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

	}; // function

	//Compute gradients inside dummy cells
	void ComputeDummyCellGradients() {
		//Current face and cell information
		Vector faceNormalL;
		Vector faceNormalR;
		Vector faceCenter;
		Vector cellCenter;

		// accumulate all usefull lyambdas
		std::vector<std::function<double(double*)> > val_comp;
		val_comp.push_back(compute.u);
		val_comp.push_back(compute.v);
		val_comp.push_back(compute.w);
		std::vector<std::function< double(std::valarray<double>) > > valc_comp;
		valc_comp.push_back(compute.uc);
		valc_comp.push_back(compute.vc);
		valc_comp.push_back(compute.wc);

		// Create entity for computind gradients by points
		GradientComputer gradcmp(nDims);

		// object for stencil creation
		Stencil t_type_stn(nDims, grid);

		////	X direction		/////
		if (!grid.IsPeriodicX) {
			faceNormalL = Vector(-1.0, 0.0, 0.0);
			faceNormalR = Vector(1.0, 0.0, 0.0);
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					//Inner cell
					cellCenter.y = grid.CoordinateY[j];
					cellCenter.z = grid.CoordinateZ[k];

					//And face
					faceCenter.y = grid.CoordinateY[j];
					faceCenter.z = grid.CoordinateZ[k];

					//Left border
					faceCenter.x = (grid.CoordinateX[grid.iMin] - 0.5 * grid.hx[grid.iMin]);
					cellCenter.x = grid.CoordinateX[grid.iMin];
					if (pManager->rankCart[0] == 0) {
						auto i = grid.iMin;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[xLeftBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalL, faceCenter, cellCenter);
							
						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::XDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i - 1, j, k));	// add dummy one cell

						// write stencil point positions
						std::vector<Vector> points;				
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i - 1, j, k);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};

					//Right border
					faceCenter.x = (grid.CoordinateX[grid.iMax] + 0.5 * grid.hx[grid.iMax]);
					cellCenter.x = grid.CoordinateX[grid.iMax];
					if (pManager->rankCart[0] == pManager->dimsCart[0] - 1) {
						auto i = grid.iMax;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[xRightBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalR, faceCenter, cellCenter);

						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::XDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i + 1, j, k));	// add dummy one cell

						// write stencil point positions
						std::vector<Vector> points;
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i + 1, j, k);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};
				};
			};
		}; //X direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Gradients dummy cells X-direction calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
		if (nDims < 2) return;

		////	Y direction		/////		
		if (!grid.IsPeriodicY) {
			// Normals
			faceNormalL = Vector(0.0, -1.0, 0.0);
			faceNormalR = Vector(0.0, 1.0, 0.0);

			// Obtain top and bottom layers
			for (int i = grid.iMin; i <= grid.iMax; i++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					//Inner cell
					cellCenter.x = grid.CoordinateX[i];
					cellCenter.z = grid.CoordinateZ[k];

					//And face
					faceCenter.x = grid.CoordinateX[i];
					faceCenter.z = grid.CoordinateZ[k];

					//Bottom border
					faceCenter.y = (grid.CoordinateY[grid.jMin] - 0.5 * grid.hy[grid.jMin]);
					cellCenter.y = grid.CoordinateY[grid.jMin];
					if (pManager->rankCart[1] == 0) {
						auto j = grid.jMin;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[yLeftBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalL, faceCenter, cellCenter);

						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::YDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i , j - 1, k));	// add one dummy cell

						// write stencil points positions
						std::vector<Vector> points;
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i, j - 1, k);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};

					//Top border
					faceCenter.y = (grid.CoordinateY[grid.jMax] + 0.5 * grid.hy[grid.jMax]);
					cellCenter.y = grid.CoordinateY[grid.jMax];
					if (pManager->rankCart[1] == pManager->dimsCart[1] - 1) {
						//Right border
						auto j = grid.jMax;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[yRightBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalR, faceCenter, cellCenter);

						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::YDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i, j + 1, k));	// add dummy one cell

						// write stencil point positions
						std::vector<Vector> points;
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i, j + 1, k);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};
				};
			};
		}; // Y direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Gradients dummy cells Y-direction calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
		if (nDims < 3) return;

		////	Z direction		/////		
		if (!grid.IsPeriodicZ) {
			// Normals
			faceNormalL = Vector(0.0, 0.0, -1.0);
			faceNormalR = Vector(0.0, 0.0, 1.0);

			// Obtain top and bottom layers
			for (int i = grid.iMin; i <= grid.iMax; i++) {
				for (int j = grid.jMin; j <= grid.jMax; j++) {
					//Inner cell
					cellCenter.x = grid.CoordinateX[i];
					cellCenter.y = grid.CoordinateY[j];

					//And face
					faceCenter.x = grid.CoordinateX[i];
					faceCenter.y = grid.CoordinateY[j];

					//Bottom border
					faceCenter.z = (grid.CoordinateZ[grid.kMin] - 0.5 * grid.hz[grid.kMin]);
					cellCenter.z = grid.CoordinateZ[grid.kMin];
					if (pManager->rankCart[2] == 0) {
						auto k = grid.kMin;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[zLeftBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalL, faceCenter, cellCenter);

						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::ZDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i, j, k - 1));	// add one dummy cell

						// write stencil points positions
						std::vector<Vector> points;
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i, j, k - 1);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};

					//Top border
					faceCenter.z = grid.CoordinateZ[grid.kMax] + 0.5 * grid.hz[grid.kMax];
					cellCenter.z = grid.CoordinateZ[grid.kMax];
					if (pManager->rankCart[2] == pManager->dimsCart[2] - 1) {
						//Right border
						auto k = grid.kMax;

						// compute conservative values at the face
						auto Vin = getCellValues(i, j, k);
						auto Vf = bConditions[zRightBC.getMarker(faceCenter)]->getFaceValues(Vin, faceNormalR, faceCenter, cellCenter);

						// get all cell indexes for inner area (create stencil in fact)
						auto cells = t_type_stn.InternalBorderStencil(Direction::ZDirection, CellIdx(i, j, k));
						cells.push_back(CellIdx(i, j, k + 1));	// add dummy one cell
						
						// write stencil point positions
						std::vector<Vector> points;
						for (auto r : cells) {
							auto p_stn = Vector(grid.CoordinateX[r.i], grid.CoordinateY[r.j], grid.CoordinateZ[r.k]);
							points.push_back(p_stn);
						};

						// compute gradients for all values given by lyambdas
						std::vector<Vector> grads;
						for (int d = 0; d < 3; d++) {
							std::vector<double> values;
							for (auto r : cells) {
								auto val_stn = getCellValues(r.i, r.j, r.k);
								values.push_back(val_comp[d](val_stn));			// set values in stencil cells
							};
							grads.push_back(gradcmp.ExecuteByLeastSquares(faceCenter, valc_comp[d](Vf), points, values));
						};

						// Extrapolate in dummy area and save results
						int idx = grid.getSerialIndexLocal(i, j , k + 1);
						int idxIn = grid.getSerialIndexLocal(i, j, k);
						gradientVelocityX[idx] = 2.0 * grads[0] - gradientVelocityX[idxIn];
						gradientVelocityY[idx] = 2.0 * grads[1] - gradientVelocityY[idxIn];
						gradientVelocityZ[idx] = 2.0 * grads[2] - gradientVelocityZ[idxIn];
					};
				};
			};
		}; // Z direction

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Gradients dummy cells Z-direction calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Compute viscous Fluxes in all faces
	void ComputeViscousFluxes(std::vector<std::vector<double>> &vfluxesX, std::vector<std::vector<double>> &vfluxesY, std::vector<std::vector<double>> &vfluxesZ) {
		//resize storage for appropriate case
		if (nDims == 1) vfluxesX.resize(grid.nlocalX + 1);
		if (nDims == 2) {
			size_t size = (grid.nlocalX + 1) * grid.nlocalY + (grid.nlocalY + 1) * grid.nlocalX;
			vfluxesX.resize(size);
			vfluxesY.resize(size);
		};
		if (nDims == 3) {
			size_t size = (grid.nlocalX + 1) * grid.nlocalY * grid.nlocalZ + (grid.nlocalY + 1) * grid.nlocalX * grid.nlocalZ + (grid.nlocalZ + 1) * grid.nlocalY * grid.nlocalX;
			vfluxesX.resize(size);
			vfluxesY.resize(size);
			vfluxesZ.resize(size);
		};

		//second viscosity
		double s_viscosity = (-2.0 / 3.0) * gas_prop.viscosity;

		//index of face
		int faceInd = 0;

		//create derivative functions
		auto getu = [](double *U) { return U[1] / U[0]; };
		auto getv = [](double *U) { return U[2] / U[0]; };
		auto getw = [](double *U) { return U[3] / U[0]; };

		//X direction first
		faceInd = 0;			/// Don't change the order in cycle
		for (int k = grid.kMin; k <= grid.kMax; k++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int i = grid.iMin; i <= grid.iMax + 1; i++) {
					//Compute vertical faces fluxes (left from cell ijk)
					//Conservative variables and serial index in left and right cell 
					double* UL = getCellValues(i - 1, j, k);
					double* UR = getCellValues(i, j, k);
					int sL = grid.getSerialIndexLocal(i - 1, j, k);
					int sR = grid.getSerialIndexLocal(i, j, k);
					double dx = grid.CoordinateX[i] - grid.CoordinateX[i - 1];
					double l = grid.hx[i - 1] / dx;	//double weights
					double r = grid.hx[i] / dx;

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = 0.5 * (r * graduL + l * graduR);
					du.x = (getu(UR) - getu(UL)) / dx;			// first order for derrivative

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = 0.5 * (r * gradvL + l * gradvR);
					dv.x = (getv(UR) - getv(UL)) / dx;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = 0.5 * (r * gradwL + l * gradwR);
					dw.x = (getw(UR) - getw(UL)) / dx;

					//compute average velocity vector
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity * (du.x + dv.y + dw.z);
					double tau_xx = tau_diagonal + 2 * gas_prop.viscosity * du.x;
					double tau_xy = gas_prop.viscosity * (du.y + dv.x);
					double tau_xz = gas_prop.viscosity * (du.z + dw.x);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaX = u * tau_xx + v * tau_xy + w * tau_xz;

					//write flux in left vertical face
					std::vector<double> vflux(5, 0);
					vflux[1] = tau_xx;
					vflux[2] = tau_xy;
					vflux[3] = tau_xz;
					vflux[4] = ThettaX;
					vfluxesX[faceInd] = vflux;
					faceInd++;
				};
			};
		};

		if (nDims < 2) return;
		//Y direction then
		faceInd = 0;			/// Don't change the order in cycle
		for (int k = grid.kMin; k <= grid.kMax; k++) {
			for (int i = grid.iMin; i <= grid.iMax; i++) {
				for (int j = grid.jMin; j <= grid.jMax + 1; j++) {

					//Compute horizontal faces fluxes (bottom from cell ijk)
					//Conservative variables and serial index in left (bottom) and right (top) cell 
					double* UL = getCellValues(i, j - 1, k);
					double* UR = getCellValues(i, j, k);
					int sL = grid.getSerialIndexLocal(i, j - 1, k);
					int sR = grid.getSerialIndexLocal(i, j, k);
					double dy = grid.CoordinateY[j] - grid.CoordinateY[j - 1];
					double l = grid.hy[j - 1] / dy;		// weights
					double r = grid.hy[j] / dy;

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = 0.5 * (r * graduL + l * graduR);
					du.y = (getu(UR) - getu(UL)) / dy;
					
					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = 0.5 * (r * gradvL + l * gradvR);
					dv.y = (getv(UR) - getv(UL)) / dy;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = 0.5 * (r * gradwL + l * gradwR);
					dw.y = (getw(UR) - getw(UL)) / dy;

					//compute average velocity vector
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));
					
					//compute stress tensor
					double tau_diagonal = s_viscosity * (du.x + dv.y + dw.z);
					double tau_yy = tau_diagonal + 2 * gas_prop.viscosity * dv.y;
					double tau_xy = gas_prop.viscosity * (du.y + dv.x);
					double tau_yz = gas_prop.viscosity * (dv.z + dw.y);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaY = u * tau_xy + v * tau_yy + w * tau_yz;

					//write flux in left vertical face
					std::vector<double> vflux(5, 0);
					vflux[1] = tau_xy;
					vflux[2] = tau_yy;
					vflux[3] = tau_yz;
					vflux[4] = ThettaY;
					vfluxesY[faceInd] = vflux;
					faceInd++;
				};
			};
		};

		if (nDims < 3) return;
		//Z direction then
		faceInd = 0;			/// Don't change the order in cycle
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax + 1; k++) {

					//Compute frontal faces fluxes (back from cell ijk)
					//Conservative variables and serial index in left (back) and right (front) cell 
					double* UL = getCellValues(i, j, k - 1);
					double* UR = getCellValues(i, j, k);
					int sL = grid.getSerialIndexLocal(i, j, k - 1);
					int sR = grid.getSerialIndexLocal(i, j, k);
					double dz = grid.CoordinateZ[k] - grid.CoordinateZ[k - 1];
					double l = grid.hz[k - 1] / dz;		//weights
					double r = grid.hz[k] / dz;

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = 0.5 * (r * graduL + l * graduR);
					du.z = (getu(UR) - getu(UL)) / dz;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = 0.5 * (r * gradvL + l * gradvR);
					dv.z = (getv(UR) - getv(UL)) / dz;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = 0.5 * (r * gradwL + l * gradwR);
					dw.z = (getw(UR) - getw(UL)) / dz;

					//compute average velocity vector
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity * (du.x + dv.y + dw.z);
					double tau_zz = tau_diagonal + 2 * gas_prop.viscosity * dw.z;
					double tau_xz = gas_prop.viscosity*(du.z + dw.x);
					double tau_yz = gas_prop.viscosity*(dv.z + dw.y);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaZ = u * tau_xz + v * tau_yz + w * tau_zz;

					//write flux in left vertical face
					std::vector<double> vflux(5, 0);
					vflux[1] = tau_xz;
					vflux[2] = tau_yz;
					vflux[3] = tau_zz;
					vflux[4] = ThettaZ;
					vfluxesZ[faceInd] = vflux;
					faceInd++;
				};
			};
		};

		return;
	};

	//================================ End of Viscous Part ====================================

	//Update all reconstruction functions for each cell
	void ComputeSolutionReconstruction() {
		//dummy reconstructions layers number
		int yLayer = 0;
		int zLayer = 0;
		if (nDims > 1) yLayer = 1;
		if (nDims > 2) zLayer = 1;

		// Obtain all inner cells			// TO DO rewrite that part using Stencil structure
		for (int i = grid.iMin; i <= grid.iMax; i++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				for (int k = grid.kMin; k <= grid.kMax; k++) {
					std::vector<std::valarray<double> > st_values;		// vector of stencil values
					std::vector<CellInfo> st_cells;						// vector of cells in stencil 

					// Create structures concerning central cell 
					std::valarray<double> cell_values = std::valarray<double>(getCellValues(i, j, k), nVariables);	// primitive variables in our cell
					CellInfo cell = grid.CreateCell(i, j, k);

					// X direction stencil
					for (int iStencil = -grid.dummyCellLayersX; iStencil <= grid.dummyCellLayersX; iStencil++) {
						if (iStencil == 0) continue;
						st_values.push_back(std::valarray<double>(getCellValues(i + iStencil, j, k), nVariables));
						st_cells.push_back(grid.CreateCell(i + iStencil, j, k));
					};

					//Y direction stencil
					for (int jStencil = -grid.dummyCellLayersY; jStencil <= grid.dummyCellLayersY; jStencil++) {
						if (jStencil == 0) continue;
						st_values.push_back(std::valarray<double>(getCellValues(i, j + jStencil, k), nVariables));
						st_cells.push_back(grid.CreateCell(i, j + jStencil, k));
					};

					//Z direction stencil
					for (int kStencil = -grid.dummyCellLayersZ; kStencil <= grid.dummyCellLayersZ; kStencil++) {
						if (kStencil == 0) continue;
						st_values.push_back(std::valarray<double>(getCellValues(i, j, k + kStencil), nVariables));
						st_cells.push_back(grid.CreateCell(i, j, k + kStencil));
					};

					//we have only one or no external layer of cells reconstructions
					reconstructions[i - grid.iMin + 1][j - grid.jMin + yLayer][k - grid.kMin + zLayer] = ComputeReconstruction<ReconstructionType>(st_values, st_cells, cell_values, cell, nDims);
				};
			};
		};

		return;
	};

	//Exchange reconstructions objects between processors and apply boundary conditions
	void ExchangeReconstructions() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;

		//layers number for dummy reconstructions
		int yLayer = 0;
		int zLayer = 0;
		int xLayer = 1;
		if (nDims > 1) yLayer = 1;
		if (nDims > 2) zLayer = 1;

		// flags wich depend from the direction considered
		int xDummyFlag = 0;
		int yDummyFlag = 0;
		int zDummyFlag = 0;

		//Get parallel run info
		MPI_Comm comm = pManager->getComm();

		//Declare buffers
		std::vector<double> bufferToRecv;
		std::vector<double> bufferToSend;

		//Allocate buffers
		auto msgLen = ReconstructionType::GetBufferLenght(nDims, nVariables);
		auto layerSizeX = grid.nlocalY * grid.nlocalZ;
		auto layerSizeY = grid.nlocalX * grid.nlocalZ;
		auto layerSizeZ = grid.nlocalX * grid.nlocalY;
		size_t bufferSize{ 0 };
		if (nDims < 2) bufferSize = msgLen;
		else if (nDims < 3) bufferSize = std::max(layerSizeX, layerSizeY) * msgLen;
		else bufferSize = std::max(std::max(layerSizeX, layerSizeY), layerSizeZ) * msgLen;
		bufferToSend.resize(bufferSize);
		bufferToRecv.resize(bufferSize);

		//! Main layer exchanging procedure //TO DO lift from lambda to member
		auto exchangeLayers = [&](Direction direction, int iSend, int rankDest, int iRecv, int rankSource) {
			size_t nRecv{ 0 };
			size_t nSend{ 0 };

			//Fill buffer with reconstructions only from inner cells (buffer to SEND reconstruction to another core)
			int idxBuffer = 0;
			for (i = grid.iMin; i <= grid.iMax; i++) {
				if ((direction == Direction::XDirection) && (i != iSend)) continue; //skip
				for (j = grid.jMin; j <= grid.jMax; j++) {
					if ((direction == Direction::YDirection) && (j != iSend)) continue; //skip
					for (k = grid.kMin; k <= grid.kMax; k++) {
						if ((direction == Direction::ZDirection) && (k != iSend)) continue; //skip

																							//Increase aticipating buffer size
						nRecv += msgLen;

						//If destination not set we don't send
						if (rankDest == -1) continue;

						//Get indexes to send reconstruction
						int iRec = i - grid.iMin + xLayer;
						int jRec = j - grid.jMin + yLayer;
						int kRec = k - grid.kMin + zLayer;

						//Serialize to string of doubles
						std::valarray<double> msg = reconstructions[iRec][jRec][kRec].Serialize();

						//Write message to buffer
						for (int i = 0; i < msgLen; i++) bufferToSend[idxBuffer * msgLen + i] = msg[i];

						//Increase object counter
						idxBuffer++;
					};
				};
			};

			//Determine number of doubles to send
			if (rankDest != -1) nSend = idxBuffer * msgLen;

			//Determine recive number											  
			if (rankSource == -1) {
				nRecv = 0;
			};

			//Make exchange
			pManager->SendRecvDouble(comm, rankDest, rankSource, &bufferToSend.front(), int(nSend), &bufferToRecv.front(), int(nRecv));

			//Write recieved values back
			idxBuffer = 0;
			switch (direction) {
			case Direction::XDirection:
				xDummyFlag = 1;
				break;
			case Direction::YDirection:
				yDummyFlag = 1;
				break;
			default:
				zDummyFlag = 1;
				break;
			};

			for (i = grid.iMin - xDummyFlag; i <= grid.iMax + xDummyFlag; i++) {
				if ((direction == Direction::XDirection) && (i != iRecv)) continue; //skip
				for (j = grid.jMin - yDummyFlag; j <= grid.jMax + yDummyFlag; j++) {
					if ((direction == Direction::YDirection) && (j != iRecv)) continue; //skip
					for (k = grid.kMin - zDummyFlag; k <= grid.kMax + zDummyFlag; k++) {
						if ((direction == Direction::ZDirection) && (k != iRecv)) continue; //skip

																							//Get indexes to recv reconstruction
						int iRec = i - grid.iMin + xLayer;
						int jRec = j - grid.jMin + yLayer;
						int kRec = k - grid.kMin + zLayer;

						//Extract message from bufer
						std::valarray<double> msg(msgLen);
						for (int i = 0; i < msgLen; i++) msg[i] = bufferToRecv[idxBuffer * msgLen + i];

						//Deserialize from string of doubles
						reconstructions[iRec][jRec][kRec].Deserialize(msg);

						//Increase object counter
						idxBuffer++;
					};
				};
			};

			// Nullify dummy flags
			xDummyFlag = 0;
			yDummyFlag = 0;
			zDummyFlag = 0;

			//Syncronize
			pManager->Barrier();
		};

		//Determine neighbours' ranks
		int rankL = pManager->GetRankByCartesianIndexShift(-1, 0, 0);
		int rankR = pManager->GetRankByCartesianIndexShift(+1, 0, 0);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		// Minus direction exchange
		int iSend = grid.iMin; // layer index to send
		int iRecv = grid.iMax + 1; // layer index to recv
		exchangeLayers(Direction::XDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = grid.iMax; // layer index to send
		iRecv = grid.iMin - 1; // layer index to recv
		exchangeLayers(Direction::XDirection, iSend, rankR, iRecv, rankL);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange recontructions X-direction executed\n";
		};

		//Sync
		pManager->Barrier();

		//Skip Y direction exchange if not needed
		if (nDims < 2) return;

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		// Minus direction exchange
		iSend = grid.jMin; // layer index to send
		iRecv = grid.jMax + 1; // layer index to recv
		exchangeLayers(Direction::YDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = grid.jMax; // layer index to send
		iRecv = grid.jMin - 1; // layer index to recv
		exchangeLayers(Direction::YDirection, iSend, rankR, iRecv, rankL);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange recontructions Y-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction exchange

		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, 0, -1);
		rankR = pManager->GetRankByCartesianIndexShift(0, 0, +1);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << rank <<
				"; rankL = " << rankL <<
				"; rankR = " << rankR <<
				std::endl << std::flush;
		};

		// Minus direction exchange
		iSend = grid.kMin; // layer index to send
		iRecv = grid.kMax + 1; // layer index to recv
		exchangeLayers(Direction::ZDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = grid.kMax; // layer index to send
		iRecv = grid.kMin - 1; // layer index to recv
		exchangeLayers(Direction::ZDirection, iSend, rankR, iRecv, rankL);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange recontructions Z-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

	}; // function

	   //Compute residual
	void ComputeResidual(const std::valarray<double>& values, std::valarray<double>& residual, std::valarray<double>& spectralRadius) {
		//  init spectral radius storage for each cell
		for (auto& sr : spectralRadius) sr = 0; //Nullify
		for (auto& r : residual) r = 0; //Nullify residual

		// Fluxes temporary storage
		std::vector<double> fl(nVariables, 0); //left flux -1/2
		std::vector<double> fr(nVariables, 0); //right flux +1/2
		std::vector<double> fvisc(nVariables, 0); //right flux +1/2		

		// Viscous part
		int faceInd = 0;
		if (isGradientRequired == true) ComputeVelocityGradients();

		std::vector<std::vector<double> > vfluxesX;
		std::vector<std::vector<double> > vfluxesY;
		std::vector<std::vector<double> > vfluxesZ;
		if (isViscousFlow == true) {
			ComputeViscousFluxes(vfluxesX, vfluxesY, vfluxesZ);
		};

		// arrays of left and right reconstructions
		std::valarray<double> UL;
		std::valarray<double> UR;

		//layers number for dummy reconstructions
		int yLayer = 0;
		int zLayer = 0;
		if (nDims > 1) yLayer = 1;
		if (nDims > 2) zLayer = 1;

		// I step
		faceInd = 0;
		Vector fn = Vector(1.0, 0.0, 0.0); // Don't change the order in cycles
		for (int k = grid.kMin; k <= grid.kMax; k++) {
			for (int j = grid.jMin; j <= grid.jMax; j++) {
				// Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = grid.hy[j];
				if (nDims == 3) fS = grid.hy[j] * grid.hz[k];

				for (int i = grid.iMin; i <= grid.iMax + 1; i++) {
					auto faceCenter = Vector(grid.CoordinateX[i] - 0.5 * grid.hx[i], grid.CoordinateY[j], grid.CoordinateZ[k]);
					auto hl = 0.5 * grid.hx[i - 1];	// distance from the L cell center to the face 
					auto hr = -0.5 * grid.hx[i];		// distance from the R cell center to the face 

					// Apply boundary conditions
					if ((pManager->rankCart[0] == 0) && (i == grid.iMin) && (grid.IsPeriodicX != true))									// Left border
					{
						UR = reconstructions[1][j - grid.jMin + yLayer][k - grid.kMin + zLayer].SampleSolution({ hr,0,0 });
						auto bcMarker = xLeftBC.getMarker(faceCenter);
						UL = bConditions[bcMarker]->getDummyReconstructions(&UR[0], fn);
					}
					else if ((pManager->rankCart[0] == pManager->dimsCart[0] - 1) && (i == grid.iMax + 1) && (grid.IsPeriodicX != true))	// Right border
					{
						UL = reconstructions[grid.iMax - grid.iMin + 1][j - grid.jMin + yLayer][k - grid.kMin + zLayer].SampleSolution({ hl,0,0 });
						auto bcMarker = xRightBC.getMarker(faceCenter);
						UR = bConditions[bcMarker]->getDummyReconstructions(&UL[0], fn);
					}
					else
					{
						UL = reconstructions[i - grid.iMin][j - grid.jMin + yLayer][k - grid.kMin + zLayer].SampleSolution({ hl,0,0 });
						UR = reconstructions[i - grid.iMin + 1][j - grid.jMin + yLayer][k - grid.kMin + zLayer].SampleSolution({ hr,0,0 });
					};

					// Compute convective flux
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
					fr = result.Fluxes;

					// Compute viscous flux
					if (isViscousFlow == true) fvisc = vfluxesX[faceInd];
					for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];

					// Correct spectral radius value from left face
					spectralRadius[ grid.getSerialIndexLocal(i, j, k) ] += fS * result.MaxEigenvalue;
					
					// Update residuals and compute timestep parameters
					if (i > grid.iMin)
					{
						// Fluxes difference equals residual
						int idx = grid.getSerialIndexLocal(i - 1, j, k);
						for (int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] -= (fr[nv] - fl[nv]) * fS;
						
						// Correct spectral radius from right face
						spectralRadius[idx] += fS * result.MaxEigenvalue;

						// TO DO check that part
						// Prepare to compute viscous timestep part
						double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
						double visc_Face = gas_prop.viscosity;
						double PrandtlNumberFace = 1.0;

						// Add up spectral radius estimate
						double vsr = std::max(4.0 / (3.0 * roFace), gas_prop.gamma / roFace);
						vsr *= visc_Face / PrandtlNumberFace;
						vsr *= fS * fS;
						vsr /= grid.volumes[idx];
						spectralRadius[idx] += vsr;
					};

					// Shift stencil
					fl = fr;

					// Increment face index
					faceInd++;
				};
			};
		};

		// J step
		if (nDims > 1)
		{
			faceInd = 0;
			Vector fn = Vector(0.0, 1.0, 0.0); // Don't change the order in cycles
			for (int k = grid.kMin; k <= grid.kMax; k++) {
				for (int i = grid.iMin; i <= grid.iMax; i++) {
					// Compute face square
					double fS = grid.hx[i];
					if (nDims == 3) fS = grid.hx[i] * grid.hz[k];

					for (int j = grid.jMin; j <= grid.jMax + 1; j++) {
						// Set geometric properties
						auto faceCenter = Vector(grid.CoordinateX[i], grid.CoordinateY[j] - 0.5 * grid.hy[j], grid.CoordinateZ[k]);
						auto hl = 0.5 * grid.hy[j - 1];		// distance from the L cell center to the face 
						auto hr = -0.5 * grid.hy[j];		// distance from the R cell center to the face 

						// Apply boundary conditions
						if ((pManager->rankCart[1] == 0) && (j == grid.jMin) && (grid.IsPeriodicY != true))									// Left border
						{
							UR = reconstructions[i - grid.iMin + 1][1][k - grid.kMin + zLayer].SampleSolution({ 0,hr,0 });
							auto bcMarker = yLeftBC.getMarker(faceCenter);
							UL = bConditions[bcMarker]->getDummyReconstructions(&UR[0], fn);
						}
						else if ((pManager->rankCart[1] == pManager->dimsCart[1] - 1) && (j == grid.jMax + 1) && (grid.IsPeriodicY != true))	// Right border
						{
							UL = reconstructions[i - grid.iMin + 1][grid.jMax - grid.jMin + 1][k - grid.kMin + zLayer].SampleSolution({ 0,hl,0 });
							auto bcMarker = yRightBC.getMarker(faceCenter);
							UR = bConditions[bcMarker]->getDummyReconstructions(&UL[0], fn);
						}
						else
						{
							UL = reconstructions[i - grid.iMin + 1][j - grid.jMin][k - grid.kMin + zLayer].SampleSolution({ 0,hl,0 });
							UR = reconstructions[i - grid.iMin + 1][j - grid.jMin + 1][k - grid.kMin + zLayer].SampleSolution({ 0,hr,0 });
						};

						// Compute convective flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
						fr = result.Fluxes;

						// Compute viscous flux
						if (isViscousFlow == true) fvisc = vfluxesY[faceInd];
						for (int nv = 0; nv < nVariables; nv++) fr[nv] -= fvisc[nv];

						// Correct spectral radius value from left face
						spectralRadius[grid.getSerialIndexLocal(i, j, k)] += fS * result.MaxEigenvalue;

						// Update residuals and compute timestep parameters
						if (j > grid.jMin)
						{
							// Fluxes difference equals residual
							int idx = grid.getSerialIndexLocal(i, j - 1, k);
							for (int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] -= (fr[nv] - fl[nv]) * fS;

							// Correct spectral radius from right face
							spectralRadius[idx] += fS * result.MaxEigenvalue;

							// TO DO check that part
							// Prepare to compute viscous timestep part
							double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
							double visc_Face = gas_prop.viscosity;
							double PrandtlNumberFace = 1.0;

							// Add up spectral radius estimate
							double vsr = std::max(4.0 / (3.0 * roFace), gas_prop.gamma / roFace);
							vsr *= visc_Face / PrandtlNumberFace;
							vsr *= fS * fS;
							vsr /= grid.volumes[idx];
							spectralRadius[idx] += vsr;
						};

						// Shift stencil
						fl = fr;

						// Increment face index
						faceInd++;
					};
				};
			};
		};

		// K step		
		if (nDims > 2)
		{
			faceInd = 0;
			Vector fn = Vector(0.0, 0.0, 1.0); // Don't change the order in cycles
			for (int i = grid.iMin; i <= grid.iMax; i++) {
				for (int j = grid.jMin; j <= grid.jMax; j++) {
					// Compute face square
					double fS = grid.hx[i] * grid.hy[j];

					for (int k = grid.kMin; k <= grid.kMax + 1; k++) {
						// Set geometric properties
						auto faceCenter = Vector(grid.CoordinateX[i], grid.CoordinateY[j], grid.CoordinateZ[k] - 0.5 * grid.hz[k]);
						auto hl = 0.5 * grid.hz[k - 1];		// distance from the L cell center to the face 
						auto hr = -0.5 * grid.hz[k];		// distance from the R cell center to the face 

															// Apply boundary conditions
						if ((pManager->rankCart[2] == 0) && (k == grid.kMin) && (grid.IsPeriodicZ != true))									// Left border
						{
							UR = reconstructions[i - grid.iMin + 1][j - grid.jMin + 1][1].SampleSolution({ 0,0,hr });
							auto bcMarker = zLeftBC.getMarker(faceCenter);
							UL = bConditions[bcMarker]->getDummyReconstructions(&UR[0], fn);
						}
						else if ((pManager->rankCart[2] == pManager->dimsCart[2] - 1) && (k == grid.kMax + 1) && (grid.IsPeriodicZ != true))	// Right border
						{
							UL = reconstructions[i - grid.iMin + 1][j - grid.jMin + 1][grid.kMax - grid.kMin + 1].SampleSolution({ 0,0,hl });
							auto bcMarker = zRightBC.getMarker(faceCenter);
							UR = bConditions[bcMarker]->getDummyReconstructions(&UL[0], fn);
						}
						else
						{
							UL = reconstructions[i - grid.iMin + 1][j - grid.jMin + 1][k - grid.kMin].SampleSolution({ 0,0,hl });
							UR = reconstructions[i - grid.iMin + 1][j - grid.jMin + 1][k - grid.kMin + 1].SampleSolution({ 0,0,hr });
						};

						// Compute convective flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
						fr = result.Fluxes;

						// Compute viscous flux
						if (isViscousFlow == true) fvisc = vfluxesZ[faceInd];
						for (int nv = 0; nv < nVariables; nv++) fr[nv] -= fvisc[nv];

						// Correct spectral radius value from left face
						spectralRadius[grid.getSerialIndexLocal(i, j, k)] += fS * result.MaxEigenvalue;

						// Update residuals and compute timestep parameters
						if (k > grid.kMin)
						{
							// Fluxes difference equals residual
							int idx = grid.getSerialIndexLocal(i, j, k - 1);
							for (int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] -= (fr[nv] - fl[nv]) * fS;

							// Correct spectral radius from right face
							spectralRadius[idx] += fS * result.MaxEigenvalue;

							// TO DO check that part
							// Prepare to compute viscous timestep part
							double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
							double visc_Face = gas_prop.viscosity;
							double PrandtlNumberFace = 1.0;

							// Add up spectral radius estimate
							double vsr = std::max(4.0 / (3.0 * roFace), gas_prop.gamma / roFace);
							vsr *= visc_Face / PrandtlNumberFace;
							vsr *= fS * fS;
							vsr /= grid.volumes[idx];
							spectralRadius[idx] += vsr;
						};

						// Shift stencil
						fl = fr;

						// Increment face index
						faceInd++;
					};
				};
			};
		};

		//Source term treatment
		if (isExternalForce == true) ProcessExternalForces();
		if (isExternalAccelaration == true) ProcessExternalAcceleration();
	};

	//Compute global time step
	double ComputeTimeStep(std::valarray<double>& spectralRadius) {

		// Compute local domain timestep
		double dt = (grid.volumes / spectralRadius).min() * CFL;

		//Check if we need to save data
		if ((SaveSolutionTime > 0) && (stepInfo.NextSolutionSnapshotTime < stepInfo.Time + dt)) dt = stepInfo.NextSolutionSnapshotTime - stepInfo.Time;
		if ((SaveSliceTime > 0) && (stepInfo.NextSliceSnapshotTime < stepInfo.Time + dt)) dt = stepInfo.NextSliceSnapshotTime - stepInfo.Time;
		
		return pManager->Min(dt);
	};

	//Explicit time step
	virtual void IterationStep() override {
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Iteration started\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		// Prepare to compute residuals	//

		// Compute all dummy values first //
		// Exchange values between processors
		ExchangeValues();
		// Apply boundary conditions for dummy cells
		ComputeDummyCellValues();

		// Reconstructions computations
		// Compute reconstruction functions for each inner cell and one extern layer
		ComputeSolutionReconstruction();
		// Exchange reconstructions objects between processors
		ExchangeReconstructions();

		// Compute residual
		ComputeResidual(values, residual, spectralRadius);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Residual calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Determine timestep
		stepInfo.TimeStep = ComputeTimeStep(spectralRadius);

		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Timestep calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Compute residuals and update solution
		UpdateSolution(stepInfo.TimeStep);

		// Update time
		stepInfo.Time += stepInfo.TimeStep;
		stepInfo.Iteration++;

		//Sync
		pManager->Barrier();
	};

};

// Global Function that create pointer to the kernel
std::unique_ptr<Kernel> CreateKernel(KernelConfiguration& conf, int argc, char *argv[]) {
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::PiecewiseConstant) {
		return std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<PiecewiseConstant>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::ENO2PointsStencil) {
		return std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM<ENO2PointsStencil>(&argc, &argv));
	};
	if (conf.methodConfiguration.ReconstructionType == Reconstruction::Linear2psLim) {
		// Check limiter
		if (conf.methodConfiguration.GeneralLimitter == LimiterType::Venkatakrishnan)
			return std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limVenkatar> >(&argc, &argv));

		if (conf.methodConfiguration.GeneralLimitter == LimiterType::BarsJespersen)
			return std::unique_ptr<Kernel>(new ExplicitRungeKuttaFVM< Linear2psLim<limBarsJespersen> >(&argc, &argv));
	};

	return std::nullptr_t();
};

#endif
