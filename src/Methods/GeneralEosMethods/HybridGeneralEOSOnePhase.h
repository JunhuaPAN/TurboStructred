#ifndef TurboStructured_Methods_GeneralEosMethods_HybridGeneralEOSOnePhase
#define TurboStructured_Methods_GeneralEosMethods_HybridGeneralEOSOnePhase

#include "Methods\HybridFVM.h"
#include "EquationsOfState\EquationsOfState.h"

//Base class for Hybrid methods (works with one phase ideal gas model)
class HybridGeneralEOSOnePhase : public HybridFVM {
public:
	//pointer to equation of state
	GeneralEOS* eos;

	//use exact expression of pressure derivative
	bool UseExactPressureDerivative;

	void SetEOS(GeneralEOS* _eos) {
		eos = _eos;
	};

	//Prepare conservative variables form left and right (relatively edge) cells
	//Prepare right eigenvectors matrix R, inverse to it - Rinv (has left eigenvectors rows) and eigenvalues
	void PrepareEigenMatrix(std::vector<double> &UL, std::vector<double> &UR, Matrix &R, Matrix &Rinv, std::vector<double> &eigenvals) override  {
		//accuracy of computer
		double eps = 1.0e-11;
		
		//Left cell
		double rol = UL[0];
		double ul  = UL[1]/rol;
		double vl = UL[2]/rol;
		double wl  = UL[3]/rol;
		double El  = UL[4]/rol;
		double el = El - 0.5*(ul*ul + vl*vl + wl*wl);	//specific internal energy
		
		//compute pressure and entalphy by eos
		double pl  = eos->GetPressure(rol, el);
		double hl = El + pl/rol;	//total entalpy

		//Right cell
		double ror = UR[0];
		double ur  = UR[1]/ror;
		double vr = UR[2]/ror;
		double wr  = UR[3]/ror;
		double Er  = UR[4]/ror;
		double er = Er - 0.5*(ur*ur + vr*vr + wr*wr);	//specific internal energy

		//compute pressure and entalphy by eos
		double pr  = eos->GetPressure(ror, er);
		double hr = Er + pr/ror;	//total entalpy

		//compute averaged values

		//values avaraged on faces
		double roa = sqrt(rol*ror);
		double sl  = sqrt(rol)/(sqrt(rol) + sqrt(ror));			//density based weight coefficients
		double sr  = sqrt(ror)/(sqrt(rol) + sqrt(ror));
		double ua = ul*sl + ur*sr;
		double va = vl*sl + vr*sr;
		double wa = wl*sl + wr*sr;
		double ha  = hl*sl + hr*sr;
		double ea = el*sl + er*sr;
		double q2  = ua*ua + va*va + wa*wa;
		double Ka = 0.5 * q2;						//kinetic energy
		double pa = roa*(ha - ea - Ka);
		double dpdea = 0;
		double dpdroa = 0;
		if(UseExactPressureDerivative == true) {
			dpdea = eos->GetPressureEnergyDerivative(roa, ea);
			dpdroa = eos->GetPressureDensityDerivative(roa, ea);
		} else {
			//Use Glaister exrpessions of approximate derivations
			//Derivative by energy
			double delta_e = er - el;
			if(delta_e < eps) {
				dpdea = 0.5*(eos->GetPressureEnergyDerivative(rol, el) + eos->GetPressureEnergyDerivative(ror, er));
			} else {
				double prr = eos->GetPressure(ror, er);
				double prl = eos->GetPressure(ror, el);
				double plr = eos->GetPressure(rol, er);
				double pll = eos->GetPressure(rol, el);
				dpdea = 0.5*(prr + plr - prl - pll);
				dpdea /= delta_e;
			};
			//Derivative by density
			double delta_ro = ror - rol;
			if(delta_ro < eps) {
				dpdroa = 0.5*(eos->GetPressureDensityDerivative(rol, el) + eos->GetPressureDensityDerivative(ror, er));
			} else {
				double prr = eos->GetPressure(ror, er);
				double prl = eos->GetPressure(ror, el);
				double plr = eos->GetPressure(rol, er);
				double pll = eos->GetPressure(rol, el);
				dpdroa = 0.5*(prr + prl - plr - pll);
				dpdroa /= delta_ro;
			};
		};
		double c2 = pa * dpdea / (roa * roa) + dpdroa;
		double c  = sqrt(c2);	

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
		R.element[4][0] = ha - c*ua;

		R.element[0][1] = 1;
		R.element[1][1] = ua;
		R.element[2][1] = va;
		R.element[3][1] = wa;
		R.element[4][1] = ea + Ka - roa * dpdroa / dpdea;		// H - c^2 / gr
		
		R.element[0][2] = 0;
		R.element[1][2] = 0;
		R.element[2][2] = 1;
		R.element[3][2] = 0;
		R.element[4][2] = va;

	    R.element[0][3] = 0;
		R.element[1][3] = 0;
		R.element[2][3] = 0;
		R.element[3][3] = 1;
		R.element[4][3] = wa;
		
		R.element[0][4] = 1;
		R.element[1][4] = ua + c;
		R.element[2][4] = va;
		R.element[3][4] = wa;
		R.element[4][4] = ha + c*ua;
		
		//grunaisen to c^2 ratio
		double grc2 = ha - ea - Ka + roa * dpdroa / dpdea;
		grc2 = 1.0 / grc2;

		//grc2 and tetta (tetta = q2 - H + 1/grc2) product
		double tetta = grc2 * q2 - grc2 * ha + 1;

		//write inverse matrix (left eigenvectors are rows)
		Rinv.element[0][0] = 0.5 * (tetta + ua / c);
		Rinv.element[0][1] = (-0.5) * (ua * grc2 + 1.0 / c);
		Rinv.element[0][2] = -0.5 * va * grc2;
		Rinv.element[0][3] = -0.5 * wa * grc2;
		Rinv.element[0][4] = 0.5 * grc2;

		Rinv.element[1][0] = grc2 * (ha - q2);
		Rinv.element[1][1] = ua*grc2;
		Rinv.element[1][2] = va*grc2;
		Rinv.element[1][3] = wa*grc2;
		Rinv.element[1][4] = -grc2;

		Rinv.element[2][0] = -va;
		Rinv.element[2][1] = 0;
		Rinv.element[2][2] = 1;
		Rinv.element[2][3] = 0;
		Rinv.element[2][4] = 0;

		Rinv.element[3][0] = -wa;
		Rinv.element[3][1] = 0;
		Rinv.element[3][2] = 0;
		Rinv.element[3][3] = 1;
		Rinv.element[3][4] = 0;

		Rinv.element[4][0] = 0.5 * (tetta - ua / c);
		Rinv.element[4][1] = (-0.5) * (ua * grc2 - 1.0 / c);
		Rinv.element[4][2] = -0.5 * va * grc2;
		Rinv.element[4][3] = -0.5 * wa * grc2;
		Rinv.element[4][4] = 0.5 * grc2;
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
					std::valarray<double>&& U = values[std::slice(sBegin, nVariables, 1)]; //std::vector<double>(values.begin() + sBegin, values.begin() + sBegin + nVariables);
					double ro = U[0];
					double rou = U[1];
					double rov = U[2];
					double row = U[3];
					double roE = U[4];
					double uu = rou/ro;
					double vv = rov/ro;
					double ww = row/ro;
					double E = roE/ro;
					double ek = (uu*uu + vv*vv + ww*ww)/2;
					double e = E - ek;
					
					//Compute sound speed by eos
					double c = eos->GetSoundSpeed(ro, e);

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
	HybridGeneralEOSOnePhase(int* argc, char **argv[]) : HybridFVM(argc, argv) {};

	//Initizalization
	virtual void Init(KernelConfiguration& config) override {
		Kernel::Init(config);
		CFL = config.methodConfiguration.CFL;
		UseExactPressureDerivative = config.methodConfiguration.UseExactPressureDerivative;
		SetEOS(config.methodConfiguration.eos);
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