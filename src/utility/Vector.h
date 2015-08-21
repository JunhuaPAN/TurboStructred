#ifndef TurboStructured_utility_Vector
#define TurboStructured_utility_Vector

#include <cmath>
#include <valarray>

//Implementation of simple 3D vector and related operations
//Vector type
class Vector {
public:
	double x;
	double y;
	double z;
	Vector(const Vector& v) : x{ v.x }, y{ v.y }, z{ v.z } {};
	Vector(): x(0), y(0), z(0)		
	{
	};
	Vector(double _x, double _y, double _z) :
		x(_x),
		y(_y),
		z(_z)
		{
		};
	double mod() const
	{
		return std::sqrt(x*x+y*y+z*z);
	};

	double& operator[](int i) {
		if (i == 0) return x;
		if (i == 1) return y;
		if (i == 2) return z;
	};

	inline const Vector& operator+=(const Vector& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	};
	inline const Vector& operator-=(const Vector& a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	};
	inline const Vector& operator*=(const double& a)
	{
		x *= a;
		y *= a;
		z *= a;
		return *this;
	};	
	inline const Vector& operator/=(const double& a)
	{
		x /= a;
		y /= a;
		z /= a;
		return *this;
	};

	operator std::valarray<double>() const { return{ x, y, z }; };
};

inline double operator*(const Vector& a, const Vector& b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
};

inline Vector operator+(const Vector& a, const Vector& b)
{
	return Vector(a.x+b.x, a.y+b.y, a.z+b.z);
};

inline Vector operator-(const Vector& a, const Vector& b)
{
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
};

inline Vector operator*(const double& a, const Vector& b)
{
	return Vector(a*b.x, a*b.y, a*b.z);
};

inline Vector operator*(const Vector& b, const double& a)
{
	return Vector(a*b.x, a*b.y, a*b.z);
};

inline Vector operator/(const Vector& b, const double& a)
{
	return Vector(b.x/a, b.y/a, b.z/a);
};

inline Vector operator-(const Vector& b)
{
	return Vector(-b.x, -b.y, -b.z);
};

inline Vector operator&(const Vector& a, const Vector &b)
{
	return Vector(a.y*b.z-b.y*a.z, b.x*a.z-a.x*b.z, a.x*b.y-a.y*b.x);
};


#endif