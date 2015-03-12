#ifndef TurboStructured_Methods_ExplicitRungeKuttaFVM
#define TurboStructured_Methods_ExplicitRungeKuttaFVM

#include "Methods\Method.h"
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
		_riemannSolver = (std::unique_ptr<RiemannSolver>)std::move(new RoeSolverPerfectGasEOS(kernelConfig.Gamma, 0.05, 0.0));
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

	//Compute all required gradient in each inner cell
	void ComputeGradients() {
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
					gradientVelocityX[idx] = ComputeGradient(i, j, k, getu);
					gradientVelocityY[idx] = ComputeGradient(i, j, k, getv);
					gradientVelocityZ[idx] = ComputeGradient(i, j, k, getw);
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
		ExchangeGradients();

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
	void ExchangeGradients() {

	};

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
		double s_viscosity = (-2.0/3.0)*viscosity;
		
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

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.x = (getu(UR) - getu(UL)) / hx;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.x = (getv(UR) - getv(UL)) / hx;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.x = (getw(UR) - getw(UL)) / hx;

					//compute average velocity vector
					double u = 0.5 * (getu(UR) + getu(UL));
					double v = 0.5 * (getv(UR) + getv(UL));
					double w = 0.5 * (getw(UR) + getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity*(du.x + dv.y + dw.z);
					double tau_xx = tau_diagonal + 2*viscosity*du.x;
					double tau_xy = viscosity*(du.y + dv.x);
					double tau_xz = viscosity*(du.z + dw.x);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaX = u*tau_xx + v*tau_xy + w*tau_xz;

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

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.y = (getu(UR) - getu(UL)) / hy;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.y = (getv(UR) - getv(UL)) / hy;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.y = (getw(UR) - getw(UL)) / hy;

					//compute average velocity vector
					double u = 0.5 * (getu(UR) + getu(UL));
					double v = 0.5 * (getv(UR) + getv(UL));
					double w = 0.5 * (getw(UR) + getw(UL));

					//compute stress tensor
					double tau_diagonal = s_viscosity*(du.x + dv.y + dw.z);
					double tau_yy = tau_diagonal + 2*viscosity*dv.y;
					double tau_xy = viscosity*(du.y + dv.x);
					double tau_yz = viscosity*(dv.z + dw.y);

					//work of viscous stresses and heat conduction (not implemented yet)
					double ThettaY = u*tau_xy + v*tau_yy + w*tau_yz;

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

					//U component derivatives
					Vector& graduL = gradientVelocityX[sL];
					Vector& graduR = gradientVelocityX[sR];
					Vector du = (graduL + graduR) / 2.0;
					du.z = (getu(UR) - getu(UL)) / hz;

					//V component derivatives
					Vector& gradvL = gradientVelocityY[sL];
					Vector& gradvR = gradientVelocityY[sR];
					Vector dv = (gradvL + gradvR) / 2.0;
					dv.z = (getv(UR) - getv(UL)) / hz;

					//W component derivatives
					Vector& gradwL = gradientVelocityZ[sL];
					Vector& gradwR = gradientVelocityZ[sR];
					Vector dw = (gradwL + gradwR) / 2.0;
					dw.z = (getw(UR) - getw(UL)) / hz;

					//compute average velocity vector
					double u = 0.5 * (getu(UR) + getu(UL));
					double v = 0.5 * (getv(UR) + getv(UL));
					double w = 0.5 * (getw(UR) + getw(UL));

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

	std::vector<double> ComputeViscousFlux(int sL, int sR, double *UL, double *UR, Vector fn) {
		std::vector<double> flux(5,0);

		//Compute gradient average on face
		double l = fn * Vector(hx, hy, hz);
		Vector t = fn;

		auto getu = [](double *U) { return U[1]/U[0]; };
		double dudl = (getu(UR) - getu(UL)) / l;
		Vector& graduL = gradientVelocityX[sL];
		Vector& graduR = gradientVelocityX[sR];
		Vector graduAvg = (graduL + graduR) / 2.0;
		Vector gradu = graduAvg - (graduAvg * t - dudl) * t;

		auto getv = [](double *U) { return U[2]/U[0]; };
		double dvdl = (getv(UR) - getv(UL)) / l;
		Vector& gradvL = gradientVelocityY[sL];
		Vector& gradvR = gradientVelocityY[sR];
		Vector gradvAvg = (gradvL + gradvR) / 2.0;
		Vector gradv = gradvAvg - (gradvAvg * t - dvdl) * t;

		Vector gradw = Vector(0,0,0);

		//Compute face velocity
		double u =  (getu(UR) + getu(UL)) / 2.0;
		double v =  (getv(UR) + getv(UL)) / 2.0;
		double w =  0; //(getw(UR) + getw(UL)) / 2.0;
		Vector velocity = Vector(u, v, w); //Face velocity
		
		//Compute viscous stress tensor
		double div = gradu.x + gradv.y;
		Matrix Tau(3, 3);
		Tau[0][0] = 2*viscosity*(gradu.x - div / 3.0); //xx
		Tau[1][1] = 2*viscosity*(gradv.y - div / 3.0); //yy
		Tau[2][2] = 2*viscosity*(gradw.z - div / 3.0); //yy
		Tau[0][1] = Tau[1][0] = viscosity*(gradu.y - gradv.x); //xy yx
		Tau[0][2] = Tau[2][0] = 0; //mu*(dudz - dwdx); //xz zx
		Tau[1][2] = Tau[2][1] = 0; //mu*(dvdz - dwdy); //yz zy

		//Compute flux
		flux[0] = 0;
		flux[1] = Tau[0] * fn;
		flux[2] = Tau[1] * fn;
		flux[3] = Tau[2] * fn;
		flux[4] = (Tau * velocity) * fn;

		return flux;
	};

	//Compute residual
	void ComputeResidual(const std::vector<double> values, std::vector<double>& residual, std::vector<double>& spectralRadius) {

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
		isGradientRequired = true;
		if (isGradientRequired == true) ComputeGradients();
		std::vector<std::vector<double> > vfluxesX;
		std::vector<std::vector<double> > vfluxesY;
		std::vector<std::vector<double> > vfluxesZ;
		ComputeViscousFluxes(vfluxesX, vfluxesY, vfluxesZ);

		// I step
		faceInd = 0;
		Vector fn = Vector(1.0, 0.0, 0.0); // x direction
		for (int k = kMin; k <= kMax; k++) {
			for (int j = jMin; j <= jMax; j++) {
				//Compute face square
				double fS = 0;
				if (nDims == 1) fS = 1.0;
				if (nDims == 2) fS = hy;
				if (nDims == 3) fS = hy * hz;

				//Initial load of values
				U[0] = getCellValues(iMin - 1, j, k);

				for (int i = iMin; i <= iMax + 1; i++) {
					//Update stencil values
					U[1] = getCellValues(i, j, k);

					//Compute convective flux
					RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
					fr = result.Fluxes;

					//Compute viscous flux
					int sL = getSerialIndexLocal(i - 1, j, k);
					int sR = getSerialIndexLocal(i, j, k);
					//fvisc = ComputeViscousFlux(sL, sR, U[0], U[1], fn);
					fvisc = vfluxesX[faceInd];
					for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
					//Update residuals
					if(i > iMin)
					{
						//Fluxes difference equals residual
						int idx = getSerialIndexLocal(i - 1, j, k);
						for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

						//Add up spectral radius estimate
						spectralRadius[idx] += fS * result.MaxEigenvalue;
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
					if (nDims == 2) fS = hx;
					if (nDims == 3) fS = hx * hz;

					//Initial load of values
					U[0] = getCellValues(i, jMin - 1, k);

					for (int j = jMin; j <= jMax + 1; j++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;

						//Compute viscous flux
						int sL = getSerialIndexLocal(i, j - 1, k);
						int sR = getSerialIndexLocal(i, j, k);
						//fvisc = ComputeViscousFlux(sL, sR, U[0], U[1], fn);
						fvisc = vfluxesY[faceInd];
						for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
						//Update residuals
						if(j > jMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j - 1, k);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
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
			fn = Vector(0.0, 0.0, 1.0); // z direction
			for (int i = iMin; i <= iMax; i++) {
				for (int j = jMin; j <= jMax; j++) {
					//Face square
					double fS = hx * hy;

					//Initial load of values
					U[0] = getCellValues(i, j, kMin - 1);

					for (int k = kMin; k <= kMax + 1; k++) {
						//Update stencil values
						U[1] = getCellValues(i, j, k);

						//Compute flux
						RiemannProblemSolutionResult result = _riemannSolver->ComputeFlux(U, fn);
						fr = result.Fluxes;

						//Compute viscous flux
						int sL = getSerialIndexLocal(i, j, k - 1);
						int sR = getSerialIndexLocal(i, j, k);
						//fvisc = ComputeViscousFlux(sL, sR, U[0], U[1], fn);
						fvisc = vfluxesZ[faceInd];
						for (int nv = 0; nv<nVariables; nv++) fr[nv] -= fvisc[nv];
	
						//Update residuals
						if(k > kMin)
						{
							//Fluxes difference equals residual
							int idx = getSerialIndexLocal(i, j, k - 1);
							for(int nv = 0; nv < nVariables; nv++) residual[idx * nVariables + nv] += - (fr[nv] - fl[nv]) * fS;

							//Add up spectral radius estimate
							spectralRadius[idx] += fS * result.MaxEigenvalue;
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

		
		double volume = hx * hy * hz;


		//right handside part
		for (int i = iMin; i <= iMax; i++)
		{
			for (int j = jMin; j <= jMax; j++)
			{
				for (int k = kMin; k <= kMax; k++)
				{
					int idx = getSerialIndexLocal(i, j, k);

					double *V = getCellValues(i, j, k);
					double u = V[1]/V[0];
					double v = V[2]/V[0];
					double w = V[3]/V[0];
					Vector velocity = Vector(u, v, w); //Cell velocity

					//Mass forces represents body accelerations (per unit mass)
					Vector mForce = Vector(0, 0, 0);
					double ro = V[0];
					residual[idx * nVariables + 1] += ro * volume * mForce.x;	//rou
					residual[idx * nVariables + 2] += ro * volume * mForce.y;	//rov
					residual[idx * nVariables + 3] += ro * volume * mForce.z;		//row
					residual[idx * nVariables + 4] += ro * mForce * velocity * volume ;	//roE

					//Potential forces
					Vector pForce = Sigma;

					//Compute total residual
					residual[idx * nVariables];										//ro
					residual[idx * nVariables + 1] += volume * pForce.x;			//rou
					residual[idx * nVariables + 2] += volume * pForce.y;			//rov
					residual[idx * nVariables + 3] += volume * pForce.z;			//row
					residual[idx * nVariables + 4] += pForce * velocity * volume;	//roE
				};
			};
		};
		//end right hand side part
	};

	//Compute global time step
	double ComputeTimeStep(std::vector<double>& spectralRadius) {
		//Compute cell volume
		double volume = hx * hy * hz;
		double dt = std::numeric_limits<double>::max();
		for (int cellIndex = 0; cellIndex< nCellsLocalAll; cellIndex++)
		{		
			double sR = spectralRadius[cellIndex];
			double localdt = CFL * volume / sR; //Blazek f. 6.20

			//Find minimum
			if (dt > localdt) dt = localdt;
		}

		dt = std::min(stepInfo.NextSnapshotTime - stepInfo.Time, dt);
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

		//Compute residual
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