#ifndef TurboStructured_Utility_Limiters
#define TurboStructured_Utility_Limiters

#include "utility\Constants.h"

/// Enumerate all available types of limiters
/*enum class LimiterType {
	BarsJespersen,		///< from Blazek CFD principles and apps 
	Venkatakrishnan		///< from Blazek
};*/

//! Base class for limiters
class Limiter {
public:
	std::valarray<double> lim_values;
	double umax, umin;			// use these vars computing lim_values

	// get limiters value vor i-th variable
	double inline operator[] (int i) {
		return lim_values[i];
	};
	double& operator() (int i) {
		return lim_values[i];
	};

	// resize function
	void inline resize(int nV) {
		lim_values.resize(nV);
	};

	// compute limiter values
	virtual void ComputeLimiterValues(std::vector<std::valarray<double> > const& values,
		std::valarray<double> const& value,
		std::vector<std::valarray<double> > const& projects,
		int nDim) {

		size_t size = value.size();
		lim_values.resize(size, 1.0);	// use 1 as default value

		// for all variables compute limiter coefficient
		double d2;
		for (auto i = 0; i < size; i++) {
			umax = value[i];
			umin = value[i];

			// find min and max values of i variable
			for (auto neig = 0; neig < 2 * nDim; neig++) {
				auto val_n = values[neig][i];
				if (umax < val_n) umax = val_n;
				if (umin > val_n) umin = val_n;
			};

			// check if here is no variation of the variable
			if (umax - umin < 3.0 * std::numeric_limits<double>::min()) continue;

			// find limiter coefficient for i variable
			auto lim_min = 1.0;
			for (auto j = 0; j < projects.size(); j++) {
				auto lim_v = 1.0;
				d2 = projects[j][i];
				if (d2 > 0) lim_v = l_function(umax, value[i], d2);
				if (d2 < 0) lim_v = l_function(umin, value[i], d2);

				// find the minimal one
				if (lim_min > lim_v) lim_min = lim_v;
			};

			// save limiter coefficient
			lim_values[i] = lim_min;
		};
	};

	// limiter function
	virtual double l_function(double um, double ui, double d2) = 0;
};

// For details see Blazek book (for example)
class limBarsJespersen : public Limiter {
public:
	// main function of that limiter:
	// um - Umin or Umax;
	// ui - U[i] in the cell;
	// d2 is difference between the reconstructed value and ui
	virtual double l_function(double um, double ui, double d2) override {
		return min(1.0, (um - ui) / d2);
	};
};

// One of the mosst popular limiter Venkatakrishnan
class limVenkatar : public Limiter {
private:
	double K;							// limiter constant
	double eps2;						// special parameter eps^2 (amount of limiting)
public:

	// Compute eps^2 value
	void SendCellInfo(CellInfo& cell, int nDim) {
		auto dx = std::pow(cell.hx * cell.hy * cell.hz, 1.0 / nDim);
		eps2 = std::pow(K * dx, 3.0);
	};

	// Main function of that limiter:
	// um - Umin or Umax;
	// ui - U[i] in the cell;
	// d2 is difference between the reconstructed value and ui
	virtual double l_function(double um, double ui, double d2) override {

		// First perform normalization procedure
		auto un = 0.5 * (umax - umin);		// normalization coefficient
		um /= un;
		ui /= un;
		d2 /= un;

		// Apply Venkatakrishnan function
		auto dm = um - ui;
		auto res =  ( (dm * dm + eps2) + 2 * d2 * dm) / (dm * dm + 2 * d2 * d2 + dm * d2 + eps2);	// (see Blazek p. 167)

		return res;
	};

	// Constructor
	limVenkatar() : K(Constants::VENKATAKR_KOEFF) {};
};


	

#endif
