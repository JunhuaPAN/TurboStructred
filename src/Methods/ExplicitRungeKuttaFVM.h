#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "KernelConfiguration.h"
#include "utility/Vector.h"
#include "utility/Matrix.h"
#include "utility/Timer.h"
#include "utility/GeomFunctions.h"
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
	std::vector<double> spectralRadius;		//array for storing spectral radiuses
	std::vector<std::vector<double>> fluxes; // array for storing fluxes
	std::vector<Vector> gradientVelocityX;  //array for storing gradients of u
	std::vector<Vector> gradientVelocityY;  //array for storing gradients of v
	std::vector<Vector> gradientVelocityZ;  //array for storing gradients of w
	std::vector< std::vector< std::vector<ReconstructionType> > > reconstructions; //array for storing 

	//Initizalization
	void Init(KernelConfiguration& kernelConfig) {
		//Invoke base class init
		Kernel::Init(kernelConfig);

		//Method specific part
		MethodConfiguration config = kernelConfig.methodConfiguration;
		assert(config.RungeKuttaOrder != 0);			// Runge-Kutta order is not initialized
		assert(config.RiemannProblemSolver != RPSolver::NoSolver);		// No solver was choosen
		CFL = config.CFL;
		RungeKuttaOrder = config.RungeKuttaOrder;
		if (config.RiemannProblemSolver == RPSolver::RoePikeSolver) {
			_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.Gamma, config.Eps, config.OperatingPresure));
		};
		if (config.RiemannProblemSolver == RPSolver::GodunovSolver) {
			_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new GodunovSolverPerfectGasEOS(kernelConfig.Gamma, config.Eps, config.OperatingPresure));
		};

		//Allocate memory for values and residual
		spectralRadius.resize(g.nCellsLocalAll);

		// Resize our reconstructions array
		reconstructions.resize(g.nlocalX + 2 * g.dummyCellLayersX);
		for (auto& r : reconstructions) r.resize(g.nlocalY + 2 * g.dummyCellLayersY);
		for (auto& r : reconstructions) {
			for (auto& s : r) s.resize(g.nlocalZ + 2 * g.dummyCellLayersZ);
		};

		// pass N Dimensions and N values to dummy cells reconstructions
		if(g.dummyCellLayersX > 0) {
			int yLayer = 0;
			int zLayer = 0;
			if (g.dummyCellLayersY > 0) yLayer = 1;
			if (g.dummyCellLayersZ > 0) zLayer = 1;
			for (int j = yLayer; j <= g.jMax - g.jMin + yLayer; j++) {
				for (int k = zLayer; k <= g.kMax - g.kMin + zLayer; k++) {
					reconstructions[0][j][k].nValues = nVariables;
					reconstructions[0][j][k].nDims = nDims;
					reconstructions[g.iMax - g.iMin + 2][j][k].nValues = nVariables;
					reconstructions[g.iMax - g.iMin + 2][j][k].nDims = nDims;
				};
			};
		};
		if(g.dummyCellLayersY > 0) {
			int xLayer = 0;
			int zLayer = 0;
			if (g.dummyCellLayersX > 0) xLayer = 1;
			if (g.dummyCellLayersZ > 0) zLayer = 1;
			for (int i = xLayer; i <= g.iMax - g.iMin + xLayer; i++) {
				for (int k = zLayer; k <= g.kMax - g.kMin + zLayer; k++) {
					reconstructions[i][0][k].nValues = nVariables;
					reconstructions[i][0][k].nDims = nDims;
					reconstructions[i][g.jMax - g.jMin + 2][k].nValues = nVariables;
					reconstructions[i][g.jMax - g.jMin + 2][k].nDims = nDims;
				};
			};
		};
		if(g.dummyCellLayersZ > 0) {
			int xLayer = 0;
			int yLayer = 0;
			if (g.dummyCellLayersX > 0) xLayer = 1;
			if (g.dummyCellLayersY > 0) yLayer = 1;
			for (int i = xLayer; i <= g.iMax - g.iMin + xLayer; i++) {
				for (int j = yLayer; j <= g.jMax - g.jMin + yLayer; j++) {
					reconstructions[i][j][0].nValues = nVariables;
					reconstructions[i][j][0].nDims = nDims;
					reconstructions[i][j][g.kMax - g.kMin + 2].nValues = nVariables;
					reconstructions[i][j][g.kMax - g.kMin + 2].nDims = nDims;
				};
			};
		};


		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//=================================   Viscous Part  ======================================

	//! Compute gradient of given function in given cell
	Vector ComputeGradient(int i, int j, int k, std::function<double(double*)> func) {
		std::vector<Vector> points;
		std::vector<double> pointValues;
		Vector center = Vector(CoordinateX[i], CoordinateY[j], CoordinateZ[k]);
		double value = func(getCellValues(i, j, k));

		//Add neighbours
		points.clear();
		pointValues.clear();
		points.push_back( Vector(CoordinateX[i-1], CoordinateY[j], CoordinateZ[k]) );
		pointValues.push_back( func(getCellValues(i-1, j, k)) );
		points.push_back( Vector(CoordinateX[i+1], CoordinateY[j], CoordinateZ[k]) );
		pointValues.push_back( func(getCellValues(i+1, j, k)) );
		if (nDims > 1) {
			points.push_back( Vector(CoordinateX[i], CoordinateY[j-1], CoordinateZ[k]) );
			pointValues.push_back( func(getCellValues(i, j-1, k)) );
			points.push_back( Vector(CoordinateX[i], CoordinateY[j+1], CoordinateZ[k]) );
			pointValues.push_back( func(getCellValues(i, j+1, k)) );
		};
		if (nDims > 2) {
			points.push_back( Vector(CoordinateX[i], CoordinateY[j], CoordinateZ[k-1]) );
			pointValues.push_back( func(getCellValues(i, j, k-1)) );
			points.push_back( Vector(CoordinateX[i], CoordinateY[j], CoordinateZ[k+1]) );
			pointValues.push_back( func(getCellValues(i, j, k+1)) );
		};
		Vector grad = ComputeGradientByPoints(nDims, center, value, points, pointValues);
		return grad;
	};

	Vector ComputeGradientUniformStructured(int i, int j, int k, std::function<double(double*)> func) {\
		Vector grad(0,0,0);
		double dux = func(getCellValues(i+1, j, k)) - func(getCellValues(i-1, j, k));
		double dudx = dux / (2.0 * hx[0]);					
		grad.x = dudx;

		if (nDims > 1) {
			double duy = func(getCellValues(i, j+1, k)) - func(getCellValues(i, j-1, k));
			double dudy = duy / (2.0 * hy[0]);
			grad.y = dudy;
		};
		if (nDims > 2) {
			double duz = func(getCellValues(i, j, k+1)) - func(getCellValues(i, j, k-1));
			double dudz = duz / (2.0 * hz[0]);
			grad.z = dudz;
		};
		return grad;
	};
	Vector ComputeGradientNonUniformStructured(int i, int j, int k, std::function<double(double*)> func) {\
		Vector grad(0,0,0);
		double dux = func(getCellValues(i + 1, j, k)) - func(getCellValues(i - 1, j, k));
		double dudx = dux / (g.CoordinateX[i + 1] - g.CoordinateX[i - 1]);
		grad.x = dudx;

		if (nDims > 1) {
			double duy = func(getCellValues(i, j + 1, k)) - func(getCellValues(i, j - 1, k));
			double dudy = duy / (g.CoordinateY[j + 1] - g.CoordinateY[j - 1]);
			grad.y = dudy;
		};
		if (nDims > 2) {
			double duz = func(getCellValues(i, j, k + 1)) - func(getCellValues(i, j, k - 1));
			double dudz = duz / (g.CoordinateZ[k + 1] - g.CoordinateZ[k - 1]);
			grad.z = dudz;
		};
		return grad;
	};

	//Compute all required gradient in each inner cell
	void ComputeVelocityGradients() {
		//Gradients required
		auto getu = [](double *U) { return U[1]/U[0]; };
		auto getv = [](double *U) { return U[2]/U[0]; };
		auto getw = [](double *U) { return U[3]/U[0]; };
		if (gradientVelocityX.size() != g.nCellsLocalAll) gradientVelocityX.resize(g.nCellsLocalAll);
		if (gradientVelocityY.size() != g.nCellsLocalAll) gradientVelocityY.resize(g.nCellsLocalAll);
		if (gradientVelocityZ.size() != g.nCellsLocalAll) gradientVelocityZ.resize(g.nCellsLocalAll);

		//Inside domain
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					int idx = getSerialIndexLocal(i, j, k);

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
			std::cout<<"rank = "<<pManager->getRank()<<", Gradients calculated inside domain\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
		
		//Interprocessor exchange
		ExchangeGradients(gradientVelocityX);
		if (nDims > 1) ExchangeGradients(gradientVelocityY);
		if (nDims > 2) ExchangeGradients(gradientVelocityZ);

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Gradients exchange executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Boundaries
		ComputeDummyCellGradients();

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Gradients calculated in dummy cells\n";
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
		int layerSizeX = g.nlocalYAll * g.nlocalZAll;
		int layerSizeY = g.nlocalXAll * g.nlocalZAll;
		int bufferSize = std::max(layerSizeX, layerSizeY) * nDims;// * nVariables;
		bufferToSend.resize(bufferSize);
		bufferToRecv.resize(bufferSize);

		SubDirection subDirection = SubDirection::Left;

		auto exchangeLayers = [&](Direction direction, int iSend, int rankDest, int iRecv, int rankSource) {
			int nRecv = 0;
			int nSend = 0;

			//Fill buffer with gradient values
			int idxBuffer = 0;
			for (i = g.iMin - g.dummyCellLayersX; i <= g.iMax + g.dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iSend)) continue; //skip
				for (j = g.jMin - g.dummyCellLayersY; j <= g.jMax + g.dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iSend)) continue; //skip
					for (k = g.kMin - g.dummyCellLayersZ; k <= g.kMax + g.dummyCellLayersZ; k++) {
						if ((direction == Direction::ZDirection) && (k != iSend)) continue; //skip
						nRecv+=nDims;
						if (rankDest == -1) continue;

						int idxCell = getSerialIndexLocal(i, j, k);
						bufferToSend[idxBuffer * nDims] = cellGradients[idxCell].x;
						if (nDims > 1) bufferToSend[idxBuffer * nDims + 1] =  cellGradients[idxCell].y;
						if (nDims > 2) bufferToSend[idxBuffer * nDims + 2] =  cellGradients[idxCell].z;
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
			for (i = g.iMin - g.dummyCellLayersX; i <= g.iMax + g.dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iRecv)) continue; //skip
				for (j = g.jMin - g.dummyCellLayersY; j <= g.jMax + g.dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iRecv)) continue; //skip
					for (k = g.kMin - g.dummyCellLayersZ; k <= g.kMax + g.dummyCellLayersZ; k++) {
						if ((direction == Direction::ZDirection) && (k != iRecv)) continue; //skip
						int idxCell = getSerialIndexLocal(i, j, k);
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
			std::cout<<"rank = "<<rank<<
				"; rankL = "<<rankL<<
				"; rankR = "<<rankR<<
				std::endl<<std::flush;
		};

		for (int layer = 1; layer <= g.dummyCellLayersX; layer++) {
			// Minus direction exchange
			int iSend = g.iMin + layer - 1; // layer index to send
			int iRecv = g.iMax + layer; // layer index to recv
			exchangeLayers(Direction::XDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = g.iMax - layer + 1; // layer index to send
			iRecv = g.iMin - layer; // layer index to recv
			exchangeLayers(Direction::XDirection, iSend, rankR, iRecv, rankL);
		};

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values X-direction executed\n";
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction exchange		
		//Determine neighbours' ranks
		rankL = pManager->GetRankByCartesianIndexShift(0, -1, 0);
		rankR = pManager->GetRankByCartesianIndexShift(0, +1, 0);
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<rank<<
				"; rankL = "<<rankL<<
				"; rankR = "<<rankR<<
				std::endl<<std::flush;
		};

		
		for (int layer = 1; layer <= g.dummyCellLayersY; layer++) {
			// Minus direction exchange
			int iSend = g.jMin + layer - 1; // layer index to send
			int iRecv = g.jMax + layer; // layer index to recv
			exchangeLayers(Direction::YDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = g.jMax - layer + 1; // layer index to send
			iRecv = g.jMin - layer; // layer index to recv
			exchangeLayers(Direction::YDirection, iSend, rankR, iRecv, rankL);
		};

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values Y-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 3) return;
		//Z direction exchange		

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Exchange values Z-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

	}; // function

	//Compute gradients inside dummy cells
	void ComputeDummyCellGradients() {
		//Index variables
		int i = 0;
		int j = 0;
		int k = 0;		

		//Current face and cell information
		Vector faceNormalL;
		Vector faceNormalR;
		Vector faceCenter;
		Vector cellCenter;

		//Gradients required
		auto getu = [](double *U) { return U[1]/U[0]; };
		auto getv = [](double *U) { return U[2]/U[0]; };
		auto getw = [](double *U) { return U[3]/U[0]; };
		
		//X direction
		faceNormalL = Vector(-1.0, 0.0, 0.0);
		faceNormalR = Vector(1.0, 0.0, 0.0);
		if (!g.IsPeriodicX) {
			for (j = g.jMin; j <= g.jMax; j++) {
				for (k = g.kMin; k <= g.kMax; k++) {
					//Inner cell
					cellCenter.y = g.CoordinateY[j];
					cellCenter.z = g.CoordinateZ[k];

					//And face
					faceCenter.y = g.CoordinateY[j];
					faceCenter.z = g.CoordinateZ[k];

					for (int layer = 1; layer <= g.dummyCellLayersX; layer++) {
						if (pManager->rankCart[0] == 0) {
							//Left border
							i = g.iMin - layer; // layer index
							int iIn = g.iMin + layer - 1; // opposite index
							cellCenter.x = g.CoordinateX[iIn];
							faceCenter.x = (g.CoordinateX[iIn] + g.CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);
					
							//Apply left boundary conditions						
							double *V = getCellValues(iIn,j,k);
							Vector dGrad;
							dGrad = xLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityX].GetDummyGradient(getu(V), gradientVelocityX[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityX[idx] = dGrad;
							dGrad = xLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityY].GetDummyGradient(getv(V), gradientVelocityY[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityY[idx] = dGrad;
							dGrad = xLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityZ].GetDummyGradient(getw(V), gradientVelocityZ[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityZ[idx] = dGrad;
						};

						if (pManager->rankCart[0] == pManager->dimsCart[0]-1) {
							//Right border
							i = g.iMax + layer; // layer index
							int iIn = g.iMax - layer + 1; // opposite index
							cellCenter.x = g.CoordinateX[iIn];
							faceCenter.x = (g.CoordinateX[iIn] + g.CoordinateX[i]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(iIn, j, k);
					
							//Apply right boundary conditions						
							double *V = getCellValues(iIn,j,k);
							Vector dGrad;
							dGrad = xRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityX].GetDummyGradient(getu(V), gradientVelocityX[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityX[idx] = dGrad;
							dGrad = xRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityY].GetDummyGradient(getv(V), gradientVelocityY[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityY[idx] = dGrad;
							dGrad = xRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityZ].GetDummyGradient(getw(V), gradientVelocityZ[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityZ[idx] = dGrad;	
						};
					};
				};
			};
		}; //X direction

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Gradients dummy cells X-direction calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		if (nDims < 2) return;
		//Y direction		
		faceNormalL = Vector(0.0, -1.0, 0.0);
		faceNormalR = Vector(0.0, 1.0, 0.0);
		for (i = g.iMin; i <= g.iMax; i++) {
			for (k = g.kMin; k <= g.kMax; k++) {
				//Inner cell
				cellCenter.x = g.CoordinateX[i];
				cellCenter.z = g.CoordinateZ[k];

				//And face
				faceCenter.x = g.CoordinateX[i];
				faceCenter.z = g.CoordinateZ[k];

				for (int layer = 1; layer <= g.dummyCellLayersY; layer++) {
					if (!g.IsPeriodicY) {
						if (pManager->rankCart[1] == 0) {
							//Left border
							j = g.jMin - layer; // layer index
							int jIn = g.jMin + layer - 1; // opposite index
							cellCenter.y = g.CoordinateY[jIn];
							faceCenter.y = (g.CoordinateY[jIn] + g.CoordinateY[j]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);
											
							//Apply left boundary conditions						
							double *V = getCellValues(i,jIn,k);
							Vector dGrad;
							dGrad = yLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityX].GetDummyGradient(getu(V), gradientVelocityX[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityX[idx] = dGrad;
							dGrad = yLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityY].GetDummyGradient(getv(V), gradientVelocityY[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityY[idx] = dGrad;
							dGrad = yLeftBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityZ].GetDummyGradient(getw(V), gradientVelocityZ[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityZ[idx] = dGrad;
						};

						if (pManager->rankCart[1] == pManager->dimsCart[1] - 1.0) {
							//Right border
							j = g.jMax + layer; // layer index
							int jIn = g.jMax - layer + 1; // opposite index
							cellCenter.y = g.CoordinateY[jIn];
							faceCenter.y = (g.CoordinateY[jIn] + g.CoordinateY[j]) / 2.0;
							int idx = getSerialIndexLocal(i, j, k);
							int idxIn = getSerialIndexLocal(i, jIn, k);
					
							//Apply right boundary conditions						
							double *V = getCellValues(i,jIn,k);
							Vector dGrad;
							dGrad = yRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityX].GetDummyGradient(getu(V), gradientVelocityX[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityX[idx] = dGrad;
							dGrad = yRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityY].GetDummyGradient(getv(V), gradientVelocityY[idxIn], faceNormalL, faceCenter, cellCenter);
							gradientVelocityY[idx] = dGrad;
							dGrad = yRightBC->boundaryConditions[BoundaryConditions::BoundaryVariableType::VelocityZ].GetDummyGradient(getw(V), gradientVelocityZ[idxIn], faceNormalL, faceCenter, cellCenter);
						};

					};
				};
			};
		}; //Y direction

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Gradients dummy cells Y-direction calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Compute viscous Fluxes in all faces
	void ComputeViscousFluxes(std::vector<std::vector<double>> &vfluxesX, std::vector<std::vector<double>> &vfluxesY, std::vector<std::vector<double>> &vfluxesZ) {
		//resize storage for appropriate case
		if(nDims == 1) vfluxesX.resize(g.nlocalX + 1 );
		if(nDims == 2) {
			vfluxesX.resize( (g.nlocalX + 1) * g.nlocalY + (g.nlocalY + 1) * g.nlocalX );
			vfluxesY.resize( (g.nlocalX + 1) * g.nlocalY + (g.nlocalY + 1) * g.nlocalX );
		};
		if(nDims == 3) {
			vfluxesX.resize( (g.nlocalX + 1) * g.nlocalY * g.nlocalZ + (g.nlocalY + 1) * g.nlocalX * g.nlocalZ + (g.nlocalZ + 1) * g.nlocalY * g.nlocalX );
			vfluxesY.resize( (g.nlocalX + 1) * g.nlocalY * g.nlocalZ + (g.nlocalY + 1) * g.nlocalX * g.nlocalZ + (g.nlocalZ + 1) * g.nlocalY * g.nlocalX );
			vfluxesZ.resize( (g.nlocalX + 1) * g.nlocalY * g.nlocalZ + (g.nlocalY + 1) * g.nlocalX * g.nlocalZ + (g.nlocalZ + 1) * g.nlocalY * g.nlocalX );
		};

		//second viscosity
		double s_viscosity = (-2.0/3.0) * viscosity;
		
		//index of face
		int faceInd = 0;
		
		//create derivative functions
		auto getu = [](double *U) { return U[1]/U[0]; };
		auto getv = [](double *U) { return U[2]/U[0]; };
		auto getw = [](double *U) { return U[3]/U[0]; };

		//X direction first
		faceInd = 0;
		for(int k = g.kMin; k <= g.kMax; k++) {
			for(int j = g.jMin; j <= g.jMax; j++) {
				for(int i = g.iMin; i <= g.iMax + 1; i++) {

					//Compute vertical faces fluxes (left from cell ijk)
					//Conservative variables and serial index in left and right cell 
					double* UL = getCellValues(i - 1, j, k);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i - 1, j, k);
					int sR = getSerialIndexLocal(i, j, k);
					double dx = g.CoordinateX[i - 1] - g.CoordinateX[i];

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.x = (getu(UR) - getu(UL)) / dx;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.x = (getv(UR) - getv(UL)) / dx;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.x = (getw(UR) - getw(UL)) / dx;

					//compute average velocity vector
					double l = 0.5 * g.hx[i - 1] / dx;	//weights
					double r = 0.5 * g.hx[i] / dx;
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));
					
					//compute stress tensor
					double tau_diagonal = s_viscosity * (du.x + dv.y + dw.z);
					double tau_xx = tau_diagonal + 2 * viscosity*du.x;
					double tau_xy = viscosity * (du.y + dv.x);
					double tau_xz = viscosity * (du.z + dw.x);

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

		if(nDims < 2) return;
		//Y direction then
		faceInd = 0;
		for(int k = g.kMin; k <= g.kMax; k++) {
			for(int i = g.iMin; i <= g.iMax; i++) {
				for(int j = g.jMin; j <= g.jMax + 1; j++) {

					//Compute horizontal faces fluxes (bottom from cell ijk)
					//Conservative variables and serial index in left (bottom) and right (top) cell 
					double* UL = getCellValues(i, j - 1, k);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i , j - 1, k);
					int sR = getSerialIndexLocal(i, j, k);
					double dy = g.CoordinateY[j] - g.CoordinateY[j - 1];

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.y = (getu(UR) - getu(UL)) / dy;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.y = (getv(UR) - getv(UL)) / dy;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.y = (getw(UR) - getw(UL)) / dy;

					//compute average velocity vector
					double l = 0.5 * g.hy[j - 1] / dy;
					double r = 0.5 * g.hy[j] / dy;
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity * (du.x + dv.y + dw.z);
					double tau_yy = tau_diagonal + 2 * viscosity*dv.y;
					double tau_xy = viscosity * (du.y + dv.x);
					double tau_yz = viscosity * (dv.z + dw.y);

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

		if(nDims < 3) return;
		//Z direction then
		faceInd = 0;
		for(int i = g.iMin; i <= g.iMax; i++) {
			for(int j = g.jMin; j <= g.jMax; j++) {
				for(int k = g.kMin; k <= g.kMax + 1; k++) {

					//Compute frontal faces fluxes (back from cell ijk)
					//Conservative variables and serial index in left (back) and right (front) cell 
					double* UL = getCellValues(i, j , k - 1);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i , j, k - 1);
					int sR = getSerialIndexLocal(i, j, k);
					double dz = g.CoordinateZ[k] - g.CoordinateZ[k - 1];

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.z = (getu(UR) - getu(UL)) / dz;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.z = (getv(UR) - getv(UL)) / dz;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.z = (getw(UR) - getw(UL)) / dz;

					//compute average velocity vector
					double l = 0.5 * g.hz[k - 1] / dz;		//weights
					double r = 0.5 * g.hz[k] / dz;
					double u = 0.5 * (l * getu(UR) + r * getu(UL));
					double v = 0.5 * (l * getv(UR) + r * getv(UL));
					double w = 0.5 * (l * getw(UR) + r * getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity*(du.x + dv.y + dw.z);
					double tau_zz = tau_diagonal + 2*viscosity*dw.z;
					double tau_xz = viscosity*(du.z + dw.x);
					double tau_yz = viscosity*(dv.z + dw.y);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaZ = u*tau_xz + v*tau_yz + w*tau_zz;

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

		// Obtain all inner cells
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					std::vector<std::valarray<double> > stencil_values;		// vector of stencil values
					std::vector<Vector> points;								// vector of stencil points
					std::valarray<double> cell_values(std::move(ConservativeToPrimitive(getCellValues(i, j, k))));	// primitive variables in our cell
					Vector cell_center = Vector(g.CoordinateX[i], g.CoordinateY[j], g.CoordinateZ[k]);

					// X direction stencil
					for (int iStencil = -g.dummyCellLayersX; iStencil <= g.dummyCellLayersX; iStencil++) {
						if (iStencil == 0) continue;
						stencil_values.push_back(std::move(ConservativeToPrimitive(getCellValues(i + iStencil, j, k))));
						points.push_back(std::move(Vector(g.CoordinateX[i + iStencil], g.CoordinateY[j], g.CoordinateZ[k])));
					};

					//Y direction stencil
					for (int jStencil = -g.dummyCellLayersY; jStencil <= g.dummyCellLayersY; jStencil++) {
						if (jStencil == 0) continue;
						stencil_values.push_back(std::move(ConservativeToPrimitive(getCellValues(i, j + jStencil, k))));
						points.push_back(std::move(Vector(g.CoordinateX[i], g.CoordinateY[j + jStencil], g.CoordinateZ[k])));
					};

					//Z direction stencil
					for (int kStencil = -g.dummyCellLayersZ; kStencil <= g.dummyCellLayersZ; kStencil++) {
						if (kStencil == 0) continue;
						stencil_values.push_back(std::move(ConservativeToPrimitive(getCellValues(i, j, k + kStencil))));
						points.push_back(std::move(Vector(g.CoordinateX[i], g.CoordinateY[j], g.CoordinateZ[k + kStencil])));
					};

					//we have only one or no external layer of cells reconstructions
					reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][k - g.kMin + zLayer] = ComputeReconstruction<ReconstructionType>(stencil_values, points, cell_values, cell_center, nDims);
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
		auto layerSizeX = g.nlocalY * g.nlocalZ;
		auto layerSizeY = g.nlocalX * g.nlocalZ;
		auto layerSizeZ = g.nlocalX * g.nlocalY;
		size_t bufferSize{ 0 };
		if (nDims < 2) bufferSize = msgLen;
		else if (nDims < 3) bufferSize = std::max(layerSizeX, layerSizeY) * msgLen;
		else bufferSize = std::max(std::max(layerSizeX, layerSizeY), layerSizeZ) * msgLen;
		bufferToSend.resize(bufferSize);
		bufferToRecv.resize(bufferSize);

		//SubDirection subDirection = SubDirection::Left;

		//! Main layer exchanging procedure //TO DO lift from lambda to member
		auto exchangeLayers = [&](Direction direction, int iSend, int rankDest, int iRecv, int rankSource) {
			size_t nRecv{ 0 };
			size_t nSend{ 0 };

			//Fill buffer with reconstructions only from inner cells (buffer to SEND reconstruction to another core)
			int idxBuffer = 0;
			for (i = g.iMin; i <= g.iMax; i++) {
				if ((direction == Direction::XDirection) && (i != iSend)) continue; //skip
				for (j = g.jMin; j <= g.jMax; j++) {
					if ((direction == Direction::YDirection) && (j != iSend)) continue; //skip
					for (k = g.kMin; k <= g.kMax; k++) {
						if ((direction == Direction::ZDirection) && (k != iSend)) continue; //skip

						//Increase aticipating buffer size
						nRecv += msgLen;

						//If destination not set we don't send
						if (rankDest == -1) continue;

						//Get indexes to send reconstruction
						int iRec = i - g.iMin + xLayer;
						int jRec = j - g.jMin + yLayer;
						int kRec = k - g.kMin + zLayer;
						
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

			for (i = g.iMin - xDummyFlag; i <= g.iMax + xDummyFlag; i++) {
				if ((direction == Direction::XDirection) && (i != iRecv)) continue; //skip
				for (j = g.jMin - yDummyFlag; j <= g.jMax + yDummyFlag; j++) {
					if ((direction == Direction::YDirection) && (j != iRecv)) continue; //skip
					for (k = g.kMin - zDummyFlag; k <= g.kMax + zDummyFlag; k++) {
						if ((direction == Direction::ZDirection) && (k != iRecv)) continue; //skip
						
						//Get indexes to recv reconstruction
						int iRec = i - g.iMin + xLayer;
						int jRec = j - g.jMin + yLayer;
						int kRec = k - g.kMin + zLayer;

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
		int iSend = g.iMin; // layer index to send
		int iRecv = g.iMax + 1; // layer index to recv
		exchangeLayers(Direction::XDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = g.iMax; // layer index to send
		iRecv = g.iMin - 1; // layer index to recv
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
		iSend = g.jMin; // layer index to send
		iRecv = g.jMax + 1; // layer index to recv
		exchangeLayers(Direction::YDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = g.jMax; // layer index to send
		iRecv = g.jMin - 1; // layer index to recv
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
		iSend = g.kMin; // layer index to send
		iRecv = g.kMax + 1; // layer index to recv
		exchangeLayers(Direction::ZDirection, iSend, rankL, iRecv, rankR);

		// Plus direction exchange
		iSend = g.kMax; // layer index to send
		iRecv = g.kMin - 1; // layer index to recv
		exchangeLayers(Direction::ZDirection, iSend, rankR, iRecv, rankL);
		if (DebugOutputEnabled) {
			std::cout << "rank = " << pManager->getRank() << ", Exchange recontructions Z-direction executed\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

	}; // function

	//Compute residual
	void ComputeResidual(const std::valarray<double>& values, std::valarray<double>& residual, std::vector<double>& spectralRadius) {
		//  init spectral radius storage for each cell
		for (double& sr : spectralRadius) sr = 0; //Nullify
		for (double& r : residual) r = 0; //Nullify residual

		// Fluxes temporary storage
		std::vector<double> fl(nVariables,0); //left flux -1/2
		std::vector<double> fr(nVariables,0); //right flux +1/2
		std::vector<double> fvisc(nVariables,0); //right flux +1/2		

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
		Vector fn = Vector(1.0, 0.0, 0.0); // x direction
		for (int k = g.kMin; k <= g.kMax; k++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				// Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = g.hy[j];
				if (nDims == 3) fS = g.hy[j] * g.hz[k];

				for (int i = g.iMin; i <= g.iMax + 1; i++) {

					// Set Riemann Problem arguments
					Vector faceCenter = Vector(g.CoordinateX[i] - 0.5 * g.hx[i], g.CoordinateY[j], g.CoordinateZ[k]);

					// Apply boundary conditions
					if ((pManager->rankCart[0] == 0) && (i == g.iMin) && (g.IsPeriodicX != true))									// Left border
					{
						UR = PrimitiveToConservative(reconstructions[1][j - g.jMin + yLayer][k - g.kMin + zLayer].SampleSolution(CubeFaces::xL), getCellValues(i, j, k));
						UL = xLeftBC->getDummyReconstructions(&UR[0]);
					}
					else if ((pManager->rankCart[0] == pManager->dimsCart[0] - 1) && (i == g.iMax + 1) && (g.IsPeriodicX != true))	// Right border
					{
						UL = PrimitiveToConservative(reconstructions[g.iMax - g.iMin + 1][j - g.jMin + yLayer][k - g.kMin + zLayer].SampleSolution(CubeFaces::xR), getCellValues(i - 1, j, k));
						UR = xRightBC->getDummyReconstructions(&UL[0]);
					} else
					{
						UL = PrimitiveToConservative(reconstructions[i - g.iMin][j - g.jMin + yLayer][k - g.kMin + zLayer].SampleSolution(CubeFaces::xR), getCellValues(i - 1, j, k));
						UR = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][k - g.kMin + zLayer].SampleSolution(CubeFaces::xL), getCellValues(i, j, k));
					};

					// TO DO DELETE
					std::valarray<double> X = PrimitiveToConservative(ConservativeToPrimitive(getCellValues(i, j, k)), getCellValues(i, j, k));
					double* x = getCellValues(i, j, k);
					X -= std::valarray<double> { x[0], x[1], x[2], x[3], x[4] };
					double sum = abs(X.sum());
					if (sum > 2.0e-10) {
						std::cout << "ERROR TRANSITION TO CHARACTERISTIC VARIABLES" << std::endl;
					};

					// Compute convective flux
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
					fr = result.Fluxes;

					// Compute viscous flux
					int sL = getSerialIndexLocal(i - 1, j, k);
					int sR = getSerialIndexLocal(i, j, k);
					if (isViscousFlow == true) fvisc = vfluxesX[faceInd];
					for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];

					// Update residuals
					if(i > g.iMin)
					{
						// Fluxes difference equals residual
						int idx = getSerialIndexLocal(i - 1, j, k);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

						// Compute volume of cell
						double volume = fS * g.hx[i - 1];

						// Add up spectral radius estimate
						spectralRadius[idx] += fS * result.MaxEigenvalue; //Convective part
						double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
						double gammaFace = gamma;
						double vsr = std::max(4.0 / (3.0 * roFace), gammaFace/roFace);
						double viscosityFace = viscosity;
						double PrandtlNumberFace = 1.0;
						vsr *= viscosityFace / PrandtlNumberFace; 
						vsr *= fS*fS;
						vsr /= volume;
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
			Vector fn = Vector(0.0, 1.0, 0.0); // y direction
			for (int k = g.kMin; k <= g.kMax; k++) {
				for (int i = g.iMin; i <= g.iMax; i++) {
					// Compute face square
					double fS = g.hx[i];
					if (nDims == 3) fS = g.hx[i] * g.hz[k];

					for (int j = g.jMin; j <= g.jMax + 1; j++) {

						// Set Riemann Problem arguments
						Vector faceCenter = Vector(g.CoordinateX[i], g.CoordinateY[j] - 0.5 * g.hy[j], g.CoordinateZ[k]);

						// Apply boundary conditions
						if ((pManager->rankCart[1] == 0) && (j == g.jMin) && (g.IsPeriodicY != true))									// Left border
						{
							UR = PrimitiveToConservative(reconstructions[i - g.iMin + 1][1][k - g.kMin + zLayer].SampleSolution(CubeFaces::yL), getCellValues(i, j, k));
							UL = yLeftBC->getDummyReconstructions(&UR[0]);
						}
						else if ((pManager->rankCart[1] == pManager->dimsCart[1] - 1) && (j == g.jMax + 1) && (g.IsPeriodicY != true))	// Right border
						{
							UL = PrimitiveToConservative(reconstructions[i - g.iMin + 1][g.jMax - g.jMin + 1][k - g.kMin + zLayer].SampleSolution(CubeFaces::yR), getCellValues(i, j - 1, k));
							UR = yRightBC->getDummyReconstructions(&UL[0]);
						}
						else
						{
							UL = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin][k - g.kMin + zLayer].SampleSolution(CubeFaces::yR), getCellValues(i, j - 1, k));
							UR = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + 1][k - g.kMin + zLayer].SampleSolution(CubeFaces::yL), getCellValues(i, j, k));
						};

						// Compute convective flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
						fr = result.Fluxes;

						// Compute viscous flux
						int sL = getSerialIndexLocal(i, j - 1, k);
						int sR = getSerialIndexLocal(i, j, k);
						if (isViscousFlow == true) fvisc = vfluxesY[faceInd];
						for (int nv = 0; nv < nVariables; nv++) fr[nv] -= fvisc[nv];

						// Update residuals
						if (j > g.jMin)
						{
							// Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j - 1, k);
							for (int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += -(fr[nv] - fl[nv]) * fS;

							// Compute volume of cell
							double volume = fS * g.hy[j - 1];

							// Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue; //Convective part
							double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
							double gammaFace = gamma;
							double vsr = std::max(4.0 / (3.0 * roFace), gammaFace / roFace);
							double viscosityFace = viscosity;
							double PrandtlNumberFace = 1.0;
							vsr *= viscosityFace / PrandtlNumberFace;
							vsr *= fS*fS;
							vsr /= volume;
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
			Vector fn = Vector(0.0, 0.0, 1.0); // y direction
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int i = g.iMin; i <= g.iMax; i++) {
					// Compute face square
					double fS = g.hx[i] * g.hy[j];

					for (int k = g.kMin; k <= g.kMax + 1; k++) {

						// Set Riemann Problem arguments
						Vector faceCenter = Vector(g.CoordinateX[i], g.CoordinateY[j], g.CoordinateZ[k] - 0.5 * g.hz[k]);

						// Apply boundary conditions
						if ((pManager->rankCart[1] == 0) && (k == g.kMin) && (g.IsPeriodicZ != true))									// Left border
						{
							UR = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][1].SampleSolution(CubeFaces::zL), getCellValues(i, j, k));
							UL = zLeftBC->getDummyReconstructions(&UR[0]);
						}
						else if ((pManager->rankCart[1] == pManager->dimsCart[1] - 1) && (k == g.kMax + 1) && (g.IsPeriodicZ != true))	// Right border
						{
							UL = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][g.kMax - g.kMin + 1].SampleSolution(CubeFaces::zR), getCellValues(i, j, k - 1));
							UR = zRightBC->getDummyReconstructions(&UL[0]);
						}
						else
						{
							UL = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][k - g.kMin].SampleSolution(CubeFaces::zR), getCellValues(i, j, k - 1));
							UR = PrimitiveToConservative(reconstructions[i - g.iMin + 1][j - g.jMin + yLayer][k - g.kMin + 1].SampleSolution(CubeFaces::zL), getCellValues(i, j, k));
						};

						// Compute convective flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(UL, UR, fn);
						fr = result.Fluxes;

						// Compute viscous flux
						int sL = getSerialIndexLocal(i, j, k - 1);
						int sR = getSerialIndexLocal(i, j, k);
						if (isViscousFlow == true) fvisc = vfluxesZ[faceInd];
						for (int nv = 0; nv < nVariables; nv++) fr[nv] -= fvisc[nv];

						// Update residuals
						if (k > g.kMin)
						{
							// Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j, k - 1);
							for (int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += -(fr[nv] - fl[nv]) * fS;

							// Compute volume of cell
							double volume = fS * g.hz[k - 1];

							// Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue; //Convective part
							double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
							double gammaFace = gamma;
							double vsr = std::max(4.0 / (3.0 * roFace), gammaFace / roFace);
							double viscosityFace = viscosity;
							double PrandtlNumberFace = 1.0;
							vsr *= viscosityFace / PrandtlNumberFace;
							vsr *= fS*fS;
							vsr /= volume;
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
		if(isExternalForce == true) ProcessExternalForces();
		if(isExternalAccelaration == true) ProcessExternalAcceleration();
	};

	//Compute global time step
	double ComputeTimeStep(std::vector<double>& spectralRadius) {
		double dt = std::numeric_limits<double>::max();
		for (int i = g.iMin; i <= g.iMax; i++) {
			for (int j = g.jMin; j <= g.jMax; j++) {
				for (int k = g.kMin; k <= g.kMax; k++) {
					//Compute cell volume
					double volume = g.hx[i] * g.hy[j] * g.hz[k];
					int idx = getSerialIndexLocal(i, j, k);
					double sR = spectralRadius[idx];
					double localdt = CFL * volume / sR; //Blazek f. 6.20

					//Find minimum
					if (dt > localdt) dt = localdt;
				};
			};
		};

		if (SaveSolutionSnapshotTime != 0) dt = std::min(stepInfo.NextSnapshotTime - stepInfo.Time, dt);
		dt = pManager->Min(dt);

		return dt;
	};

	//Explicit time step
	virtual void IterationStep() override {	
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Iteration started\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		// Prepare to compute residuals

		//Exchange values between processors
		ExchangeValues();

		// Apply boundary conditions for dummy cells
		ComputeDummyCellValues();

		// Compute reconstruction functions for each inner cell and one extern layer
		ComputeSolutionReconstruction();

		//Exchange reconstructions objects between processors
		ExchangeReconstructions();

		//Compute residual TO DO second argument remove
		ComputeResidual(values, residual, spectralRadius);

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Residual calculated\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();

		//Determine timestep
		stepInfo.TimeStep = ComputeTimeStep(spectralRadius);
		
		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Timestep calculated\n";
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

#endif
