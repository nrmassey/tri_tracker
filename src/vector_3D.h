/******************************************************************************
** Program : vector_3D.h
** Author  : Neil Massey
** Date    : 23/04/09
** Purpose : 3D vector class
******************************************************************************/

#ifndef VECTOR_3D
#define VECTOR_3D

#include <iostream>

class vector_3D
{
	friend bool point_in_tri(const vector_3D& p, const class tri_3D& T);
	friend std::ostream& operator<<(std::ostream& os, const vector_3D& V);
	friend std::istream& operator>>(std::istream& is, vector_3D& V);

	public:
		vector_3D(void);
		vector_3D(FP_TYPE x, FP_TYPE y, FP_TYPE z);
		vector_3D(const vector_3D& rhs);

		FP_TYPE& operator[](const int i);
		const FP_TYPE& operator[](const int i) const;

		vector_3D operator+(const vector_3D& rhs) const;
		vector_3D operator+(const FP_TYPE sc) const;
		void operator+=(const vector_3D& rhs);
		void operator+=(const FP_TYPE sc);

		vector_3D operator-(const vector_3D& rhs) const;
		vector_3D operator-(const FP_TYPE sc) const;
		void operator-=(const vector_3D& rhs);
		void operator-=(const FP_TYPE sc);

		vector_3D operator*(const FP_TYPE& sc) const;
		void operator*=(const FP_TYPE& sc);

		vector_3D operator/(const FP_TYPE& sc) const;
		void operator/=(const FP_TYPE& sc);

		bool operator==(const vector_3D& rhs) const;
		bool operator!=(const vector_3D& rhs) const;

		FP_TYPE mag(void) const;					// magnitude
		FP_TYPE sqmag(void) const;					// square of the magnitude
		FP_TYPE dp(const vector_3D& rhs) const;		// dot product
		vector_3D xp(const vector_3D& rhs) const;	// cross product
		void zero(void);							// reset to zero

	protected:
		FP_TYPE x;
		FP_TYPE y;
		FP_TYPE z;
};

#endif
