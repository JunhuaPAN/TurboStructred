#ifndef TurboStructured_Utility_Limiters
#define TurboStructured_Utility_Limiters

//#include <cmath>
//#include <iostream>
//#include <vector>

/// Enumerate all available types of limiters
/*enum class LimiterType {
	BarsJespersen,		///< from Blazek CFD principles and apps 
	Venkatakrishnan		///< from Blazek
};*/

//! Base class for limiters
class Limiter {
public:
	std::valarray<double> lim_values;

	// get limiters value vor i-th variable
	int inline operator[] (int i) {
		return lim_values[i];
	};

	// compute limiter values
	virtual void ComputeLimiterValues(std::vector<std::valarray<double> > const& values,
		std::valarray<double> const& value,
		std::vector<std::valarray<double> > const& projects,
		int nDim) = 0;
};

// For details see for example Blazek book 
class limBarsJespersen : public Limiter {
public:
	// compute limiter values
	virtual void ComputeLimiterValues(std::vector<std::valarray<double> > const& values,
		std::valarray<double> const& value,
		std::vector<std::valarray<double> > const& projects,
		int nDim) {

		size_t size = value.size();
		lim_values.resize(size);	// use 1 as default value
		
		// for all variables compute limiter coefficient
		double d2, umax, umin;
		for (auto i = 0; i < size; i++) {
			auto umax = value[i], umin = value[i];

			// find min and max values of i variable
			for (auto neig = 0; neig < 2 * nDim; neig++) {
				auto val_n = values[neig][i];
				if (umax < val_n) umax = val_n;
				if (umin > val_n) umin = val_n;
			};

			// find limiter coefficient for i variable
			auto lim_min = 1.0;
			for (auto j = 0; j < projects.size(); j++) {
				auto lim_v = 1.0;
				d2 = projects[j][i];
				if (d2 > 0) lim_v = min(1.0, (umax - value[i]) / d2);
				if (d2 < 0) lim_v = min(1.0, (umin - value[i]) / d2);

				// find the minimal one
				if (lim_min > lim_v) lim_min = lim_v;
			};

			// save limiter coefficient
			lim_values[i] = lim_min;
		};
	};
};


#endif
