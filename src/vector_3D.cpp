/******************************************************************************
** Program : vector_3D.cpp
** Author  : Neil Massey
** Date    : 23/04/09
** Purpose : 3D vector class
******************************************************************************/

#include <assert.h>
#include <math.h>
#include "vector_3D.h"

vector_3D::vector_3D(void)
{
	x = y = z = 0.0;
}

/*****************************************************************************/

vector_3D::vector_3D(FP_TYPE px, FP_TYPE py, FP_TYPE pz)
{
	x = px;
	y = py;
	z = pz;
}

/*****************************************************************************/

vector_3D::vector_3D(const vector_3D& rhs)
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
}

/*****************************************************************************/

FP_TYPE& vector_3D::operator[](const int i)	
{
	assert(i >=0 && i < 3);
	switch(i)
	{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
	}
	return z;
}

/*****************************************************************************/

const FP_TYPE& vector_3D::operator[](const int i) const
{
	assert(i >=0 && i < 3);
	switch(i)
	{
		case 0:
			return x;
		case 1:
			return y;
		case 2:
			return z;
	}
	return z;
}

/*****************************************************************************/

vector_3D vector_3D::operator+(const vector_3D& rhs) const
{
	return vector_3D(x+rhs.x, y+rhs.y, z+rhs.z);
}

/*****************************************************************************/

vector_3D vector_3D::operator+(const FP_TYPE sc) const
{
	return vector_3D(sc+x, sc+y, sc+z);
}

/*****************************************************************************/

void vector_3D::operator+=(const vector_3D& rhs)
{
	x += rhs.x;
	y += rhs.y;
	z += rhs.z;
}

/*****************************************************************************/

void vector_3D::operator+=(const FP_TYPE sc)
{
	x += sc;
	y += sc;
	z += sc;
}

/*****************************************************************************/

vector_3D vector_3D::operator-(const vector_3D& rhs) const
{
	return vector_3D(x-rhs.x, y-rhs.y, z-rhs.z);
}

/*****************************************************************************/

vector_3D vector_3D::operator-(const FP_TYPE sc) const
{
	return vector_3D(x-sc, y-sc, z-sc);
}

/*****************************************************************************/

void vector_3D::operator-=(const vector_3D& rhs)
{
	x -= rhs.x;
	y -= rhs.y;
	z -= rhs.z;
}

/*****************************************************************************/

void vector_3D::operator-=(const FP_TYPE sc)
{
	x -= sc;
	y -= sc;
	z -= sc;
}

/*****************************************************************************/

vector_3D vector_3D::operator*(const FP_TYPE& sc) const
{
	return vector_3D(x*sc, y*sc, z*sc);
}

/*****************************************************************************/

void vector_3D::operator*=(const FP_TYPE& sc)
{
	x *= sc;
	y *= sc;
	z *= sc;
}

/*****************************************************************************/

vector_3D vector_3D::operator/(const FP_TYPE& sc) const
{
	return vector_3D(x/sc, y/sc, z/sc);
}

/*****************************************************************************/

void vector_3D::operator/=(const FP_TYPE& sc)
{
	x /= sc;
	y /= sc;
	z /= sc;
}

/*****************************************************************************/

bool vector_3D::operator==(const vector_3D& rhs) const
{
	return x == rhs.x && y == rhs.y && z == rhs.z;
}

/*****************************************************************************/

bool vector_3D::operator!=(const vector_3D& rhs) const
{
	return x != rhs.x && y != rhs.y && z != rhs.z;
}

/*****************************************************************************/

FP_TYPE vector_3D::mag(void) const
{
	FP_TYPE l = x*x + y*y + z*z;
	return sqrt(l);
}

/*****************************************************************************/

FP_TYPE vector_3D::sqmag(void) const
{
	FP_TYPE l = x*x + y*y + z*z;
	return l;
}

/*****************************************************************************/

FP_TYPE vector_3D::dp(const vector_3D& rhs) const
{
	return x*rhs.x + y*rhs.y + z*rhs.z;
}

/*****************************************************************************/

vector_3D vector_3D::xp(const vector_3D& rhs) const
{
	return vector_3D(y*rhs.z - z*rhs.y, 
					 z*rhs.x - x*rhs.z,
					 x*rhs.y - y*rhs.x);
}

/*****************************************************************************/

void vector_3D::zero(void)
{
	x = y = z = 0.0f;
}

/*****************************************************************************/

std::ostream& operator<<(std::ostream& os, const vector_3D& V)
{
	os << V.x << " " << V.y << " " << V.z;
	return os;
}

/*****************************************************************************/

std::istream& operator>>(std::istream& is, vector_3D& V)
{
	is >> V.x >> V.y >> V.z;
	return is;
}
