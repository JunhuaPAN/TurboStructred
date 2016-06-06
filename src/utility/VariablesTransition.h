#ifndef TurboStructured_Utility_VariablesTransition
#define TurboStructured_Utility_VariablesTransition

#include <cmath>
#include <iostream>
#include <vector>
#include "utility/Direction.h"
#include "GasProperties.h"

// Collect usefull functions that computes primitive values
class ValuesTransition {
private:
	//Gas model parameters
	std::unique_ptr<GasProperties> gas_prop;

public:
	// lambdas for velocity computing
	std::function<double(double*)> u{ [](double* U) {
		return U[1]/U[0];
	} };
	std::function<double(double*)> v{ [](double* U) {
		return U[2] / U[0];
	} };
	std::function<double(double*)> w{ [](double* U) {
		return U[3] / U[0];
	} };

	// from conservative variables
	std::function<double(std::valarray<double>)> uc{ [](const std::valarray<double>& U) {
		return U[1] / U[0];
	} };
	std::function<double(std::valarray<double>)> vc{ [](const std::valarray<double>& U) {
		return U[2] / U[0];
	} };
	std::function<double(std::valarray<double>)> wc{ [](const std::valarray<double>& U) {
		return U[3] / U[0];
	} };

	// initialise lymbdas that depends on gas properties
	void BindGasProperties(GasProperties& _gas_prop) {
		gas_prop = std::make_unique<GasProperties>(_gas_prop);
	};
};

// Transition to characteristic variables averaged by Roe procedure
class RoeLineariser {
private:
	Direction myDir;

	// variables in cell number
	const size_t size = 5;

	// double specific gas constant
	double gamma;

	// Roe average state variables
	double ro, H, c;
	Vector vel;

	// Left eigenvetors matrix (rows are left eigenvetcors)
	std::vector<std::valarray<double>> Rinv;

	// Right one (columns are rignt eigenvectors)
	std::vector<std::valarray<double>> R;

	// Rotate frame of reference for velocity vector
	inline std::valarray<double> ExchangeVelocityComponents(std::valarray<double>& value) {
		// default case (don't need to exchange anything)
		if (myDir == Direction::XDirection) {
			return value;
		};

		// space case
		std::valarray<double> res = value;
		if (myDir == Direction::YDirection) {
			res[1] = value[2];
			res[2] = value[3];
			res[3] = value[1];
		};
		if (myDir == Direction::ZDirection) {
			res[1] = value[3];
			res[2] = value[1];
			res[3] = value[2];
		};
		return res;
	};
	inline std::valarray<double> InverseExchangeVelocityComponents(std::valarray<double>& value) {
		// default case (don't need to exchange anything)
		if (myDir == Direction::XDirection) {
			return value;
		};

		// space case
		std::valarray<double> res = value;
		if (myDir == Direction::YDirection) {
			res[1] = value[3];
			res[2] = value[1];
			res[3] = value[2];
		};
		if (myDir == Direction::ZDirection) {
			res[1] = value[2];
			res[2] = value[3];
			res[3] = value[1];
		};
		return res;
	};

	// compute averaged variables
	void ComputeRoeAverageState(std::valarray<double>& _UL, std::valarray<double>& _UR) {
		std::valarray<double> UL = ExchangeVelocityComponents(_UL);
		std::valarray<double> UR = ExchangeVelocityComponents(_UR);

		// left and right density
		double ro_l = UL[0];
		double ro_r = UR[0];

		// velocities
		Vector velocity_l, velocity_r;
		velocity_l.x = UL[1] / ro_l;
		velocity_l.y = UL[2] / ro_l;
		velocity_l.z = UL[3] / ro_l;
		velocity_r.x = UR[1] / ro_r;
		velocity_r.y = UR[2] / ro_r;
		velocity_r.z = UR[3] / ro_r;
		
		// enthalpy
		double e_l = UL[4] / ro_l;
		double e_r = UR[4] / ro_r;
		double k;
		k = 0.5 * (velocity_l * velocity_l);	// kinetik energy
		double h_l = (e_l - k) * gamma + k;		// left enthalpy
		k = 0.5 * (velocity_r * velocity_r);	// kinetik energy
		double h_r = (e_r - k) * gamma + k;		// right enthalpy

		// averaging weights
		double ql = sqrt(ro_l) / (sqrt(ro_l) + sqrt(ro_r));
		double qr = sqrt(ro_r) / (sqrt(ro_l) + sqrt(ro_r));

		// Roe averaged state 
		ro = sqrt(ro_l*ro_r);								// (Roe averaged) density		
		vel = ql*velocity_l + qr*velocity_r;				// (Roe averaged) velocity	
		H = ql*h_l + qr*h_r;								// (Roe averaged) total enthalpy
		c = sqrt((gamma - 1) * (H - 0.5 * (vel * vel)));	// acoustic velocity
		assert(c>0);

		return;
	};

	// compute the matrix that inverse to R matrix 
	void ComputeLeftEigenVectors() {
		double K = 0.5 * vel * vel;
		double b1 = (gamma - 1.0) / (c * c);
		double b2 = K * b1;

		// compute left eigenvectors of Jacobian
		Rinv.resize(size);
		Rinv[0] = { 0.5 * (b2 + vel.x / c), -0.5 * (b1 * vel.x + 1.0 / c), -0.5 * b1 * vel.y, -0.5 * b1 * vel.z, 0.5 * b1 };
		Rinv[1] = { -vel.y, 0.0, 1.0, 0.0, 0.0 };
		Rinv[2] = { -vel.z, 0.0, 0.0, 1.0, 0.0 };
		Rinv[3] = { (1.0 - b2), b1 * vel.x, b1 * vel.y, b1 * vel.z, -b1 };
		Rinv[4] = { 0.5 * (b2 - vel.x / c), -0.5 * (b1 * vel.x - 1.0 / c), -0.5 * b1 * vel.y, -0.5 * b1 * vel.z, 0.5 * b1 };
	};

	// compute R matrix
	void ComputeRightEigenVectors() {
		double K = 0.5 * (vel * vel);
		// compute rows of R matrix
		R.resize(size);

		R[0] = { 1.0, 0.0, 0.0, 1.0, 1.0 };
		R[1] = { vel.x - c, 0.0, 0.0, vel.x, vel.x + c };
		R[2] = { vel.y, 1.0, 0.0, vel.y, vel.y };
		R[3] = { vel.z, 0.0, 1.0, vel.z, vel.z };
		R[4] = { H - vel.x * c, vel.y, vel.z, K, H + vel.x * c };
	};

public:

	//constructor
	RoeLineariser(size_t _size, double _gamma) : size(_size), gamma(_gamma) {};

	// Set the direction of transition 
	void SetDirection(Direction dir) {
		myDir = dir;
	};

	// Compute average values and R matrix and the inverse one
	void PrepareTransitionMatrix(std::valarray<double>& UL, std::valarray<double>& UR) {
		ComputeRoeAverageState(UL, UR);
		ComputeLeftEigenVectors();
		ComputeRightEigenVectors();
	};

	// complete transition to characteristic variables for several cells
	std::vector<std::valarray<double> > GroupTransition(std::vector<std::valarray<double> >& values) {
		// characteristic variables in stencil
		std::vector<std::valarray<double> > res;
		for (auto _r : values) {
			std::valarray<double> CharVal(size);
			auto r = ExchangeVelocityComponents(_r);
			for (int i = 0; i < size; i++) CharVal[i] = (Rinv[i] * r).sum();
			res.push_back(CharVal);
		};
		return std::move(res);
	};

	// do the same for one cell
	std::valarray<double> SingleTransition(std::valarray<double>& _value) {
		// convert values in our cell
		std::valarray<double> CharValue(size);
		auto value = ExchangeVelocityComponents(_value);
		for (int i = 0; i < size; i++) CharValue[i] = (Rinv[i] * value).sum();
		return std::move(CharValue);
	};

	// complete inverse transition
	std::valarray<double> InverseTransition(std::valarray<double>& value) {
		// convert values in our cell
		std::valarray<double> CharValue(size);
		for (int i = 0; i < size; i++) CharValue[i] = (R[i] * value).sum();
		CharValue = InverseExchangeVelocityComponents(CharValue);

		return std::move(CharValue);
	};
};

#endif
