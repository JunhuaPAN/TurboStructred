#ifndef TurboStructured_Methods_GeneralEosMethods_HybridBarotropicEOSOnePhase
#define TurboStructured_Methods_GeneralEosMethods_HybridBarotropicEOSOnePhase

#include "Methods\HybridFVM.h"
#include "EquationsOfState\EquationsOfState.h"
#include "BCGeneralBarotropic.h"

//Base class for Hybrid methods (works with one phase ideal gas model)
class HybridBarotropicEOSOnePhase : public HybridFVM {
public:
	//pointer to equation of state
	BarotropicEOS* eos;

	void SetEOS(BarotropicEOS* _eos) {
		eos = _eos;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) override  {
		//Left cell
		double rol = UL[0];
		double ul  = UL[1]/rol;
		double vl = UL[2]/rol;
		double wl  = UL[3]/rol;

		//Right cell
		double ror = UR[0];
		double ur  = UR[1]/ror;
		double vr = UR[2]/ror;
		double wr  = UR[3]/ror;

		//values averaged on faces
		double roa = sqrt(rol*ror);
		double sl  = sqrt(rol)/(sqrt(rol) + sqrt(ror));			//density based weight coefficients
		double sr  = sqrt(ror)/(sqrt(rol) + sqrt(ror));
		double ua = ul*sl + ur*sr;
		double va = vl*sl + vr*sr;
		double wa = wl*sl + wr*sr;

		double c = eos->GetSoundSpeed(roa,  0);
		double c2  = c * c;	

		//write eigenvalues
		eigenvals[0] = ua - c;
		eigenvals[1] = ua;
		eigenvals[2] = ua;
		eigenvals[3] = ua + c;

		//write right eigenvectors
		R.element[0][0] = 1;
		R.element[1][0] = ua - c;
		R.element[2][0] = va;
		R.element[3][0] = wa;

		R.element[0][1] = 0;
		R.element[1][1] = 0;
		R.element[2][1] = 1;
		R.element[3][1] = 0;

	    R.element[0][2] = 0;
		R.element[1][2] = 0;
		R.element[2][2] = 0;
		R.element[3][2] = 1;
		
		R.element[0][3] = 1;
		R.element[1][3] = ua + c;
		R.element[2][3] = va;
		R.element[3][3] = wa;
		
		//write inverse matrix (left eigenvectors are rows)
		Rinv.element[0][0] = 0.5 * (c + ua) / c;
		Rinv.element[0][1] = -0.5 / c;
		Rinv.element[0][2] = 0;
		Rinv.element[0][3] = 0;

		Rinv.element[1][0] = -va;
		Rinv.element[1][1] = 0;
		Rinv.element[1][2] = 1;
		Rinv.element[1][3] = 0;

		Rinv.element[2][0] = -wa;
		Rinv.element[2][1] = 0;
		Rinv.element[2][2] = 0;
		Rinv.element[2][3] = 1;

		Rinv.element[3][0] = 0.5 * (c - ua) / c;
		Rinv.element[3][1] = 0.5 / c;
		Rinv.element[3][2] = 0;
		Rinv.element[3][3] = 0;
	};
	 
	//compute full tine step
	virtual void ComputeTimeStep() override {
		//variables to collect maximum values
		double dmax  = 0;
		double ccmax = 0;	//maximum of speed of sound (not used yet)

		//conservative variable temporal container		
		for (int ix = iMin; ix <= iMax; ix++)
		{
			for (int iy = jMin; iy <= jMax; iy++)
			{
				for (int iz = kMin; iz <= kMax; iz++)
				{
					int sBegin = getSerialIndexLocal(ix, iy, iz) * nVariables;
					std::valarray<double>&& U = values[std::slice(sBegin, nVariables, 1)];
					double ro = U[0];
					double rou = U[1];
					double rov = U[2];
					double row = U[3];
					double uu = rou/ro;
					double vv = rov/ro;
					double ww = row/ro;
					
					//Compute sound speed by eos
					double c = eos->GetSoundSpeed(ro, 0);

					double um = fabs(uu) + c;
					ccmax = std::max(ccmax, c);
					um /= hx[1];	//todo
					if(um >= dmax) dmax = um;
          
					if (nDims > 1) {
						um = fabs(vv) + c;
						um /= hy[1];	//todo
						if(um >= dmax) dmax = um;
					};

					if (nDims > 2) {
						um = fabs(ww) + c;
						um /= hz[1];	//todo
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
	HybridBarotropicEOSOnePhase(int* argc, char **argv[]) : HybridFVM(argc, argv) {
		nVariables = 4;
	};

	//Initizalization
	virtual void Init(KernelConfiguration& config) override {
		Kernel::Init(config);
		CFL = config.methodConfiguration.CFL;
		SetEOS(dynamic_cast<BarotropicEOS*>(config.methodConfiguration.eos));
	};

	//Init boundary conditions
	virtual void InitBoundaryConditions(KernelConfiguration& config) override {
		if (!IsPeriodicX) {
			xLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
			xRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
			xLeftBC->loadConfiguration(config.xLeftBoundary);
			xRightBC->loadConfiguration(config.xRightBoundary);
		};
		if ((!IsPeriodicY) && (nDims > 1)) {
			yLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
			yRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
			yLeftBC->loadConfiguration(config.yLeftBoundary);
			yRightBC->loadConfiguration(config.yRightBoundary);
		};
		if ((!IsPeriodicZ) && (nDims > 2)) {
			zLeftBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
			zRightBC = std::unique_ptr<BoundaryConditions::BCGeneralBarotropic>(new BoundaryConditions::BCGeneralBarotropic());
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
			ofs<<"\""<<"e"<<"\" ";
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
				double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
				double P = eos->GetPressure(ro, e);

				//Write to file
				ofs<<x<<" ";
				ofs<<ro<<" ";
				ofs<<u<<" ";
				ofs<<v<<" ";
				ofs<<w<<" ";
				ofs<<P<<" ";
				ofs<<e<<" ";
				ofs<<std::endl;
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
			ofs<<"\""<<"e"<<"\" ";
			ofs<<std::endl;
			
			ofs << "ZONE T=\"1\"\nI=" << nX <<"J=" << nY << "K=" << nZ << "F=POINT\n";
			ofs << "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)\n";
		
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
						double e = U[4] / ro - 0.5*(u*u + v*v + w*w);
						double P = eos->GetPressure(ro, e);

						//Write to file
						ofs<<x<<" ";
						ofs<<y<<" ";
						ofs<<z<<" ";
						ofs<<ro<<" ";
						ofs<<u<<" ";
						ofs<<v<<" ";
						ofs<<w<<" ";
						ofs<<P<<" ";
						ofs<<e<<" ";
						ofs<<std::endl;
					};
				};
			};	//end for ijk
		};	//end if

		ofs.close();
		return;
	};

};

#endif