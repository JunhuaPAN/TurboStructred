#ifndef TurboStructured_Methods_GeneralEosMethods_HybridBarotropicEOSTwoPhase
#define TurboStructured_Methods_GeneralEosMethods_HybridBarotropicEOSTwoPhase

#include "Methods\HybridFVM.h"
#include "EquationsOfState\EquationsOfState.h"
#include "BCGeneralBarotropic.h"

//Base class for Hybrid methods (works with one phase ideal gas model)
class HybridBarotropicEOSTwoPhase : public HybridFVM {
public:
	//pointer to equation of state
	BarotropicTwoPhaseEOS* eos;

	//use exact expression of pressure derivative
	bool UseExactPressureDerivative;

	void SetEOS(BarotropicTwoPhaseEOS* _eos) {
		eos = _eos;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) override  {
		//accuracy of computer
		double eps = 1.0e-13;
		
		//Left cell
		double rol = UL[0];
		double ul  = UL[1]/rol;
		double vl = UL[2]/rol;
		double wl  = UL[3]/rol;
		double alphal = UL[4]/rol;

		//Right cell
		double ror = UR[0];
		double ur  = UR[1]/ror;
		double vr = UR[2]/ror;
		double wr  = UR[3]/ror;
		double alphar = UL[4]/rol;

		//Values averaged on faces
		double roa = sqrt(rol*ror);
		double sl  = sqrt(rol)/(sqrt(rol) + sqrt(ror));			//Density based weight coefficients
		double sr  = sqrt(ror)/(sqrt(rol) + sqrt(ror));
		double ua = ul*sl + ur*sr;
		double va = vl*sl + vr*sr;
		double wa = wl*sl + wr*sr;
		double alpha = alphal*sl + alphar*sr;

		//Pressure derivatives
		double dpdalpha = 0;
		double dpdro = 0;
		if(UseExactPressureDerivative == true) {
			dpdalpha = eos->GetPressureVolumeFractionDerivative(roa, alpha);
			dpdro = eos->GetPressureDensityDerivative(roa, alpha);
		} else {
			//Use Glaister exrpessions of approximate derivations
			//Derivative by energy
			double delta_alpha = alphar - alphal;
			if(delta_alpha < eps) {
				dpdalpha = 0.5*(eos->GetPressureVolumeFractionDerivative(rol, alphal) + eos->GetPressureVolumeFractionDerivative(ror, alphar));
			} else {
				double prr = eos->GetPressure(ror, alphar);
				double prl = eos->GetPressure(ror, alphal);
				double plr = eos->GetPressure(rol, alphar);
				double pll = eos->GetPressure(rol, alphal);
				dpdalpha = 0.5*(prr + plr - prl - pll);
				dpdalpha /= delta_alpha;
			};
			//Derivative by density
			double delta_ro = ror - rol;
			if(delta_ro < eps) {
				dpdro = 0.5*(eos->GetPressureDensityDerivative(rol, alphal) + eos->GetPressureDensityDerivative(ror, alphar));
			} else {
				double prr = eos->GetPressure(ror, alphar);
				double prl = eos->GetPressure(ror, alphal);
				double plr = eos->GetPressure(rol, alphar);
				double pll = eos->GetPressure(rol, alphal);
				dpdro = 0.5*(prr + prl - plr - pll);
				dpdalpha /= delta_ro;
			};
		};
		double c2 = dpdro;
		double c  = sqrt(c2);
		double gr = dpdalpha / roa;	//some variable similar to grunisen coefficient
		double gc2 = gr / c2;	//just usefull constant

		//write eigenvalues
		eigenvals[0] = ua - c;
		eigenvals[1] = ua;
		eigenvals[2] = ua;
		eigenvals[3] = ua;
		eigenvals[4] = ua + c;

		//write right eigenvectors
		R.element[0][0] = 1;
		R.element[1][0] = ua - c;
		R.element[2][0] = va;
		R.element[3][0] = wa;
		R.element[4][0] = alpha;

		R.element[0][1] = gr;
		R.element[1][1] = gr * ua;
		R.element[2][1] = 0;
		R.element[3][1] = 0;
		R.element[4][1] = gr * alpha - c2;

	    R.element[0][2] = 0;
		R.element[1][2] = 0;
		R.element[2][2] = 1;
		R.element[3][2] = 0;
		R.element[4][2] = 0;

		R.element[0][3] = 0;
		R.element[1][3] = 0;
		R.element[2][3] = 0;
		R.element[3][3] = 1;
		R.element[4][3] = 0;
		
		R.element[0][4] = 1;
		R.element[1][4] = ua + c;
		R.element[2][4] = va;
		R.element[3][4] = wa;
		R.element[4][4] = alpha;
		
		//write inverse matrix (left eigenvectors are rows)
		Rinv.element[0][0] = 0.5 * (c2 + ua * c - alpha * gr) / c2;
		Rinv.element[0][1] = -0.5 / c;
		Rinv.element[0][2] = 0;
		Rinv.element[0][3] = 0;
		Rinv.element[0][4] = 0.5 * gc2;

		Rinv.element[1][0] = alpha / c2;
		Rinv.element[1][1] = 0;
		Rinv.element[1][2] = 0;
		Rinv.element[1][3] = 0;
		Rinv.element[1][4] = -1.0 / c2;

		Rinv.element[2][0] = va * ( alpha * gc2 - 1.0 );
		Rinv.element[2][1] = 0;
		Rinv.element[2][2] = 1;
		Rinv.element[2][3] = 0;
		Rinv.element[2][4] = -va * gc2;

		Rinv.element[3][0] = wa * ( alpha * gc2 - 1.0 );
		Rinv.element[3][1] = 0;
		Rinv.element[3][2] = 0;
		Rinv.element[3][3] = 1;
		Rinv.element[3][4] = -wa * gc2;

		Rinv.element[4][0] = 0.5 * (c2 - ua * c - alpha * gr) / c2;
		Rinv.element[4][1] = 0.5 / c;
		Rinv.element[4][2] = 0;
		Rinv.element[4][3] = 0;
		Rinv.element[4][4] = 0.5 * gc2;

		Matrix test = R*Rinv;
	};
	 
	//compute full tine step
	virtual void ComputeTimeStep() override {
		//variables to collect maximum values
		double dmax  = 0;
		double ccmax = 0;	//maximum of speed of sound (not used yet)

		//conservative variable temporal container
		std::vector<double> U;
	
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				for (int iz = kMin; iz <= kMax; iz++)
				{
					int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
					U =  std::vector<double>(values.begin() + sBegin, values.begin() + sBegin + nVariables);
					double ro = U[0];
					double rou = U[1];
					double rov = U[2];
					double row = U[3];
					double uu = rou/ro;
					double vv = rov/ro;
					double ww = row/ro;
					double alpha = U[4]/ro;
					
					//Compute sound speed by eos
					double c = eos->GetSoundSpeed(ro, alpha);

					double um = fabs(uu) + c;
					ccmax = std::max(ccmax, c);
					um /= hx;
					if(um >= dmax) dmax = um;
          
					if (nDims > 1) {
						um = fabs(vv) + c;
						um /= hy;
						if(um >= dmax) dmax = um;
					};

					if (nDims > 2) {
						um = fabs(ww) + c;
						um /= hz;
						if(um >= dmax) dmax = um; 
					};
				};
			};
		};

		double dt = CFL/dmax;
		if(stepInfo.NextSnapshotTime>0) dt = std::min(stepInfo.NextSnapshotTime - stepInfo.Time, dt);
		stepInfo.TimeStep = dt;
	};
	 
	//Constuctor inherited
	HybridBarotropicEOSTwoPhase(int* argc, char **argv[]) : HybridFVM(argc, argv) {};

	//Initizalization
	virtual void Init(KernelConfiguration& config) override {
		Kernel::Init(config);
		CFL = config.methodConfiguration.CFL;
		UseExactPressureDerivative = config.methodConfiguration.UseExactPressureDerivative;
		SetEOS(dynamic_cast<BarotropicTwoPhaseEOS*>(config.methodConfiguration.eos));
	};

	//Init boundary conditions
	virtual void InitBoundaryConditions(KernelConfiguration& config) override {
		if (!IsPeriodicX) {
			xLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			xRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			yRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};
		if ((!IsPeriodicZ) && (nDims > 2)) {
			zLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			zRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropicTwoPhase>(new BoundaryConditions::BCGeneralBarotropicTwoPhase());
			zLeftBC->loadConfiguration(config.zLeftBoundary);
			zRightBC->loadConfiguration(config.zRightBoundary);
		};
	};

	//Save solution to TecPlot
	virtual void SaveSolutionSega(std::string fname) override {
		//Tecplot version
		std::ofstream ofs(fname);

		//1D tecplot style
		if(nDims == 1)
		{
			//Header				
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"X"<<"\" ";
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			ofs<<"\""<<"v"<<"\" ";
			ofs<<"\""<<"w"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"alpha1"<<"\" ";
			ofs<<"\""<<"alpha2"<<"\" ";
			ofs<<std::endl;
			
			//Solution
			for (int i = iMin; i <= iMax; i++) {
				//Obtain cell data
				double x = CoordinateX[i];
				double* U = getCellValues(i, jMin, kMin);
				double ro = U[0];
				double u = U[1] / ro;
				double v = U[2] / ro;
				double w = U[3] / ro;
				double alpha = U[4] / ro;
				double P = eos->GetPressure(ro, alpha);

				//Write to file
				ofs << x << " ";
				ofs << ro << " ";
				ofs << u << " ";
				ofs << v << " ";
				ofs << w << " ";
				ofs << P << " ";
				ofs << alpha << " ";
				ofs << 1.0 - alpha;
				ofs << std::endl;
			};	//end cycle

			ofs.close();
			return;
		};	//end if

		//2D/3D tecplot style
		if(nDims > 1)
		{
			//Header				
			ofs<<"VARIABLES = ";
			ofs<<"\""<<"X"<<"\" ";
			ofs<<"\""<<"Y"<<"\" ";
			ofs<<"\""<<"Z"<<"\" ";
			ofs<<"\""<<"ro"<<"\" ";
			ofs<<"\""<<"u"<<"\" ";
			ofs<<"\""<<"v"<<"\" ";
			ofs<<"\""<<"w"<<"\" ";
			ofs<<"\""<<"P"<<"\" ";
			ofs<<"\""<<"alpha1"<<"\" ";
			ofs<<"\""<<"alpha2"<<"\" ";
			ofs<<std::endl;
			
			ofs << "ZONE T=\"1\"\nI=" << nX <<"J=" << nY << "K=" << nZ << "F=POINT\n";
			ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
		
			//Solution
			for (int k = kMin; k <= kMax; k++) {
				for (int j = jMin; j <= jMax; j++) {
					for (int i = iMin; i <= iMax; i++) {
						//Obtain cell data
						double x = CoordinateX[i];
						double y = CoordinateY[j];
						double z = CoordinateZ[k];
						double* U = getCellValues(i,j,k);
						double ro = U[0];
						double u = U[1] / ro;
						double v = U[2] / ro;
						double w = U[3] / ro;
						double alpha = U[4] / ro;
						double P = eos->GetPressure(ro, alpha);

						//Write to file
						ofs << x << " ";
						ofs << y << " ";
						ofs << z << " ";
						ofs << ro << " ";
						ofs << u << " ";
						ofs << v << " ";
						ofs << w << " ";
						ofs << P << " ";
						ofs << alpha << " ";
						ofs << 1.0 - alpha;
						ofs << std::endl;
					};
				};
			};	//end for ijk
		};	//end if

		ofs.close();
		return;
	};

};

#endif