#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "KernelConfiguration.h"
#include "utility\Vector.h"
#include "utility\Matrix.h"
#include "utility\Timer.h"
#include "utility\GeomFunctions.h"
#include "RiemannSolvers\RoeSolverPerfectGasEOS.h"
#include "kernel.h"

//Base class for all solution methods that desribe iterations process in detail
class ExplicitRungeKuttaFVM : public Kernel {
	//Current riemann problem solver
	std::unique_ptr<RiemannSolver> _riemannSolver;
public:
	//Inherit constructor
	//using Kernel::Kernel;
	ExplicitRungeKuttaFVM(int* argc, char **argv[]) : Kernel(argc, argv) {};

	//Method parameters
	int RungeKuttaOrder;
	double CFL;	

	//Additional data
	std::vector<double> spectralRadius;		//array for storing spectral radiuses
	std::vector<std::vector<double>> fluxes; // array for storing fluxes
	std::vector<Vector> gradientVelocityX;  //array for storing gradients of u
	std::vector<Vector> gradientVelocityY;  //array for storing gradients of v
	std::vector<Vector> gradientVelocityZ;  //array for storing gradients of w


	//Number of required dummy layers
	virtual int GetDummyCellLayerSize() override {
		return 1;
	};

	//Initizalization
	void Init(KernelConfiguration& kernelConfig) {
		//Invoke base class init
		Kernel::Init(kernelConfig);

		//Method specific part
		MethodConfiguration config = kernelConfig.methodConfiguration;
		_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.Gamma, config.Eps, 0.0));
		CFL = config.CFL;
		RungeKuttaOrder = config.RungeKuttaOrder;
		
		//Allocate memory for values and residual
		spectralRadius.resize(nCellsLocalAll);	

		if (DebugOutputEnabled) {
			std::cout<<"rank = "<<pManager->getRank()<<", Kernel initialized\n";
			std::cout.flush();
		};
		//Sync
		pManager->Barrier();
	};

	//Viscous Part
	//Compute gradient of given function in given cell
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
		double dudx = dux / (CoordinateX[i + 1] - CoordinateX[i - 1]);					
		grad.x = dudx;

		if (nDims > 1) {
			double duy = func(getCellValues(i, j + 1, k)) - func(getCellValues(i, j - 1, k));
			double dudy = duy / (CoordinateY[j + 1] - CoordinateY[j - 1]);
			grad.y = dudy;
		};
		if (nDims > 2) {
			double duz = func(getCellValues(i, j, k + 1)) - func(getCellValues(i, j, k - 1));
			double dudz = duz / (CoordinateZ[k + 1] - CoordinateZ[k - 1]);
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
		if (gradientVelocityX.size() != nCellsLocalAll) gradientVelocityX.resize(nCellsLocalAll);
		if (gradientVelocityY.size() != nCellsLocalAll) gradientVelocityY.resize(nCellsLocalAll);
		if (gradientVelocityZ.size() != nCellsLocalAll) gradientVelocityZ.resize(nCellsLocalAll);

		//Inside domain
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
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
		int rankLeft; // if equals -1 no left neighbour
		int rankRight; // if equals -1 no right neighbour
		
		//X direction exchange		

		//Allocate buffers
		int layerSizeX = nlocalYAll * nlocalZAll;
		int layerSizeY = nlocalXAll * nlocalZAll;
		int bufferSize = std::max(layerSizeX, layerSizeY) * nDims;// * nVariables;
		bufferToSend.resize(bufferSize);
		bufferToRecv.resize(bufferSize);

		SubDirection subDirection = SubDirection::Left;

		auto exchangeLayers = [&](Direction direction, int iSend, int rankDest, int iRecv, int rankSource) {
			int nRecv = 0;
			int nSend = 0;

			//Fill buffer with gradient values
			int idxBuffer = 0;
			for (i = iMin - dummyCellLayersX; i <= iMax + dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iSend)) continue; //skip
				for (j = jMin - dummyCellLayersY; j <= jMax + dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iSend)) continue; //skip
					for (k = kMin - dummyCellLayersZ; k <= kMax + dummyCellLayersZ; k++) {
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
			for (i = iMin - dummyCellLayersX; i <= iMax + dummyCellLayersX; i++) {
				if ((direction == Direction::XDirection) && (i != iRecv)) continue; //skip
				for (j = jMin - dummyCellLayersY; j <= jMax + dummyCellLayersY; j++) {
					if ((direction == Direction::YDirection) && (j != iRecv)) continue; //skip
					for (k = kMin - dummyCellLayersZ; k <= kMax + dummyCellLayersZ; k++) {
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

		for (int layer = 1; layer <= dummyCellLayersX; layer++) {
			// Minus direction exchange
			int iSend = iMin + layer - 1; // layer index to send
			int iRecv = iMax + layer; // layer index to recv
			exchangeLayers(Direction::XDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = iMax - layer + 1; // layer index to send
			iRecv = iMin - layer; // layer index to recv
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

		
		for (int layer = 1; layer <= dummyCellLayersY; layer++) {
			// Minus direction exchange
			int iSend = jMin + layer - 1; // layer index to send
			int iRecv = jMax + layer; // layer index to recv
			exchangeLayers(Direction::YDirection, iSend, rankL, iRecv, rankR);

			// Plus direction exchange
			iSend = jMax - layer + 1; // layer index to send
			iRecv = jMin - layer; // layer index to recv
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
		if (!IsPeriodicX) {
			for (j = jMin; j <= jMax; j++) {
				for (k = kMin; k <= kMax; k++) {				
					//Inner cell
					cellCenter.y = CoordinateY[j];
					cellCenter.z = CoordinateZ[k];

					//And face
					faceCenter.y = CoordinateY[j];
					faceCenter.z = CoordinateZ[k];

					for (int layer = 1; layer <= dummyCellLayersX; layer++) {		
						if (pManager->rankCart[0] == 0) {
							//Left border
							i = iMin - layer; // layer index
							int iIn = iMin + layer - 1; // opposite index
							cellCenter.x = CoordinateX[iIn];
							faceCenter.x = (CoordinateX[iIn] + CoordinateX[i]) / 2.0;
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
							i = iMax + layer; // layer index
							int iIn = iMax - layer + 1; // opposite index
							cellCenter.x = CoordinateX[iIn];
							faceCenter.x = (CoordinateX[iIn] + CoordinateX[i]) / 2.0;
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
		for (i = iMin; i <= iMax; i++) {
			for (k = kMin; k <= kMax; k++) {				
				//Inner cell
				cellCenter.x = CoordinateX[i];
				cellCenter.z = CoordinateZ[k];

				//And face
				faceCenter.x = CoordinateX[i];
				faceCenter.z = CoordinateZ[k];

				for (int layer = 1; layer <= dummyCellLayersY; layer++) {
					if (!IsPeriodicY) {
						if (pManager->rankCart[1] == 0) {
							//Left border
							j = jMin - layer; // layer index
							int jIn = jMin + layer - 1; // opposite index
							cellCenter.y = CoordinateY[jIn];
							faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;
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
							j = jMax + layer; // layer index
							int jIn = jMax - layer + 1; // opposite index
							cellCenter.y = CoordinateY[jIn];
							faceCenter.y = (CoordinateY[jIn] + CoordinateY[j]) / 2.0;
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
		if(nDims == 1) vfluxesX.resize( nlocalX + 1 );
		if(nDims == 2) {
			vfluxesX.resize( (nlocalX + 1) * nlocalY + (nlocalY + 1) * nlocalX );
			vfluxesY.resize( (nlocalX + 1) * nlocalY + (nlocalY + 1) * nlocalX );
		};
		if(nDims == 3) {
			vfluxesX.resize( (nlocalX + 1) * nlocalY * nlocalZ + (nlocalY + 1) * nlocalX * nlocalZ + (nlocalZ + 1) * nlocalY * nlocalX );
			vfluxesY.resize( (nlocalX + 1) * nlocalY * nlocalZ + (nlocalY + 1) * nlocalX * nlocalZ + (nlocalZ + 1) * nlocalY * nlocalX );
			vfluxesZ.resize( (nlocalX + 1) * nlocalY * nlocalZ + (nlocalY + 1) * nlocalX * nlocalZ + (nlocalZ + 1) * nlocalY * nlocalX );
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
		for(int k = kMin; k <= kMax; k++) {
			for(int j = jMin; j <= jMax; j++) {
				for(int i = iMin; i <= iMax + 1; i++) {

					//Compute vertical faces fluxes (left from cell ijk)
					//Conservative variables and serial index in left and right cell 
					double* UL = getCellValues(i - 1, j, k);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i - 1, j, k);
					int sR = getSerialIndexLocal(i, j, k);
					double dx = CoordinateX[i - 1] - CoordinateX[i];

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
					double l = 0.5 * hx[i - 1] / dx;	//weights
					double r = 0.5 * hx[i] / dx;
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
		for(int k = kMin; k <= kMax; k++) {
			for(int i = iMin; i <= iMax; i++) {
				for(int j = jMin; j <= jMax + 1; j++) {

					//Compute horizontal faces fluxes (bottom from cell ijk)
					//Conservative variables and serial index in left (bottom) and right (top) cell 
					double* UL = getCellValues(i, j - 1, k);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i , j - 1, k);
					int sR = getSerialIndexLocal(i, j, k);
					double dy = CoordinateY[j] - CoordinateY[j - 1];

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
					double l = 0.5 * hy[j - 1] / dy;
					double r = 0.5 * hy[j] / dy;
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
		for(int i = iMin; i <= iMax; i++) {
			for(int j = jMin; j <= jMax; j++) {
				for(int k = kMin; k <= kMax + 1; k++) {

					//Compute frontal faces fluxes (back from cell ijk)
					//Conservative variables and serial index in left (back) and right (front) cell 
					double* UL = getCellValues(i, j , k - 1);
					double* UR = getCellValues(i, j, k);
					int sL = getSerialIndexLocal(i , j, k - 1);
					int sR = getSerialIndexLocal(i, j, k);
					double dz = CoordinateZ[k] - CoordinateZ[k - 1];

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
					double l = 0.5 * hz[k - 1] / dz;		//weights
					double r = 0.5 * hz[k] / dz;
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

	//Compute residual
	void ComputeResidual(const std::valarray<double>& values, std::valarray<double>& residual, std::vector<double>& spectralRadius) {
		//  init spectral radius storage for each cell
		for (double& sr : spectralRadius) sr = 0; //Nullify
		for (double& r : residual) r = 0; //Nullify residual

		// array of pointers to cell values
		std::vector<double*> U(2);
		
		// Fluxes temporary storage
		std::vector<double> fl(5,0); //left flux -1/2
		std::vector<double> fr(5,0); //right flux +1/2
		std::vector<double> fvisc(5,0); //right flux +1/2

		//viscous part
		int faceInd = 0;
		if (isGradientRequired == true) ComputeVelocityGradients();

		std::vector<std::vector<double> > vfluxesX;
		std::vector<std::vector<double> > vfluxesY;
		std::vector<std::vector<double> > vfluxesZ;
		if (isViscousFlow == true) {
			ComputeViscousFluxes(vfluxesX, vfluxesY, vfluxesZ);
		};

		// I step
		faceInd = 0;
		Vector fn = Vector(1.0, 0.0, 0.0); // x direction
		for (int k = kMin; k <= kMax; k++) {
			for (int j = jMin; j <= jMax; j++) {
				//Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = hy[j];	//todo
				if (nDims == 3) fS = hy[j] * hz[k];	//todo

				//Initial load of values
				U[0] = getCellValues(iMin - 1, j, k);

				for (int i = iMin; i <= iMax + 1; i++) {
					//Update stencil values
					U[1] = getCellValues(i, j, k);

					//Compute convective flux		
					std::vector<double> UL(U[0], U[0] + nVariables);
					std::vector<double> UR(U[1], U[1] + nVariables);
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
					fr = result.Fluxes;

					//Compute viscous flux
					int sL = getSerialIndexLocal(i - 1, j, k);
					int sR = getSerialIndexLocal(i, j, k);
					if (isViscousFlow == true) fvisc = vfluxesX[faceInd];
					for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
					//Update residuals
					if(i > iMin)
					{
						//Fluxes difference equals residual
						int idx = getSerialIndexLocal(i - 1, j, k);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

						//Compute volume of cell
						double volume = fS * hx[i - 1];

						//Add up spectral radius estimate
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
          
					//Shift stencil
					fl = fr;
					U[0] = U[1];

					//Increment face index
					faceInd++;
				 };
			};
		};
		
		// J step
		if (nDims > 1)
		{
			faceInd = 0;
			fn = Vector(0.0, 1.0, 0.0); // y direction
			for (int k = kMin; k <= kMax; k++) {
				for (int i = iMin; i <= iMax; i++) {
					//Compute face square
					double fS = 0;
					if (nDims == 2) fS = hx[i];
					if (nDims == 3) fS = hx[i] * hz[k];

					//Initial load of values
					U[0] = getCellValues(i, jMin - 1, k);

					for (int j = jMin; j <= jMax + 1; j++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						std::vector<double> UL(U[0], U[0] + nVariables);
						std::vector<double> UR(U[1], U[1] + nVariables);
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;

						//Compute viscous flux
						int sL = getSerialIndexLocal(i, j - 1, k);
						int sR = getSerialIndexLocal(i, j, k);
						if (isViscousFlow == true) fvisc = vfluxesY[faceInd];
						for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
						//Update residuals
						if(j > jMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j - 1, k);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Compute volume of cell
							double volume = fS * hy[j - 1]; // todo

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
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
          
						//Shift stencil
						fl = fr;
						U[0] = U[1];

						//Increment face index
						faceInd++;
					};
				};
			};
		};
		
		// K step
		if (nDims > 2)
		{
			faceInd = 0;
			fn = Vector(0.0, 0.0, 1.0); // y direction
			for (int i = iMin; i <= iMax; i++) {
				for (int j = jMin; j <= jMax; j++) {
					//Compute face square
					double fS = hx[i] * hy[j];

					//Initial load of values
					U[0] = getCellValues(i, j, kMin - 1);

					for (int k = kMin; k <= kMax + 1; k++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						std::vector<double> UL(U[0], U[0] + nVariables);
						std::vector<double> UR(U[1], U[1] + nVariables);
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;

						//Compute viscous flux
						int sL = getSerialIndexLocal(i, j, k - 1);
						int sR = getSerialIndexLocal(i, j, k);
						if (isViscousFlow == true) fvisc = vfluxesZ[faceInd];
						for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
						//Update residuals
						if(k > kMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j, k - 1);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Compute volume of cell
							double volume = fS * hz[k - 1];

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
							double roFace = sqrt(UL[0] * UR[0]); // (Roe averaged) density
							double gammaFace = gamma;
							double vsr = std::max(4.0 / (3.0 * roFace), gammaFace/roFace);
							double viscosityFace = viscosity;
							double PrandtlNumberFace = 1.0;
							vsr *= viscosityFace / PrandtlNumberFace; 
							vsr *= fS * fS;
							vsr /= volume;
							spectralRadius[idx] += vsr;
						};
          
						//Shift stencil
						fl = fr;
						U[0] = U[1];

						//Increment face index
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
		for (int i = iMin; i <= iMax; i++) {
			for (int j = jMin; j <= jMax; j++) {
				for (int k = kMin; k <= kMax; k++) {
					//Compute cell volume
					double volume = hx[i] * hy[j] * hz[k];
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

		//Compute residual   TO DO second argument remove
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