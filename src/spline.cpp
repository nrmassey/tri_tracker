/******************************************************************************
** Program : spline.cpp
** Author  : Neil Massey
** Date    : 24/09/09
** Purpose : Spline class with control points that vary in t
******************************************************************************/

#include "spline.h"
#include <math.h>
#include <iostream>

/*****************************************************************************/

spline::spline(std::vector<FP_TYPE> v, std::vector<FP_TYPE> ix, FP_TYPE imv)
{
	// resize the coefficient vectors
	a.resize(v.size());
	b.resize(v.size());
	c.resize(v.size());
	d.resize(v.size());
	x = ix;
	calc_coefficients(v);
	mv = imv;
}

/*****************************************************************************/

FP_TYPE spline::evaluate(FP_TYPE t)
{
	int i = find_piece(t);
	FP_TYPE v;
	// check for missing value
	if (a[i] == mv)
		v = mv;
	else
	{
		FP_TYPE t1 = t - x[i];
		v = pow(t1, 3.0)*d[i] + pow(t1, 2.0)*c[i] + t1*b[i] + a[i];
	}
	return v;
}

/*****************************************************************************/

FP_TYPE spline::evaluate_dx(FP_TYPE t)
{
	int i = find_piece(t);
	FP_TYPE v;
	// check for mv
	if (a[i] == mv)
		v = mv;
	else
	{
		FP_TYPE t1 = t - x[i];
		v = 3.0*pow(t1, 2.0)*d[i] + 2.0*t1*c[i] + b[i];
	}
	return v;
}

/*****************************************************************************/

FP_TYPE spline::evaluate_d2x(FP_TYPE t)
{
	int i = find_piece(t);
	FP_TYPE v;
	// check for mv
	if (a[i] == mv)
		v = mv;
	else
	{
		FP_TYPE t1 = t - x[i];
		v = 6.0*t1*d[i]+2.0*t1*c[i];
	}
	return v;
}

/*****************************************************************************/

int spline::find_piece(FP_TYPE t)
{
	int ti=0;
	for (int i=0; i<x.size()-1; i++)
	{
		if (x[i] <= t && t <= x[i+1])
			break;
		ti++;
	}
	return ti;
}

/*****************************************************************************/

void spline::calc_coefficients(std::vector<FP_TYPE> v)
{
	// check that the size of the point and value vectors match
	if (v.size() != x.size())
		throw("Intervals and value vector lengths do not match");

	int n = x.size();
	// calculate the interval values
	std::vector<FP_TYPE> h(n, 0.0);
	for (int i=0; i<n-1; i++)
		h[i] = x[i+1] - x[i];

	// assign the starting offset
	for (int i=0; i<n; i++)
		a[i] = v[i];

	std::vector<FP_TYPE> alpha(n, 0.0);
	for (int i=1; i<n-1; i++)
		alpha[i] = 3.0/h[i]*(a[i+1]-a[i]) - 3.0/h[i-1]*(a[i]-a[i-1]);

	// solve the tridiagonal system
	std::vector<FP_TYPE>  l(n, 0.0);
	std::vector<FP_TYPE> mu(n, 0.0);
	std::vector<FP_TYPE>  z(n, 0.0);
	l[0] = 1.0;

	// elimination phase
	for (int i=1; i<n-1; i++)
	{
		l[i]  = 2.0*(x[i+1] - x[i-1]) - h[i-1]*mu[i-1];
		mu[i] = h[i]/l[i];
		z[i]  = (alpha[i] - h[i-1] * z[i-1])/l[i];
	}

	// substitution phase
	l[n-1] = 1.0;
	z[n-1] = 0.0;
	c[n-1] = 0.0;

	for (int j=n-2; j>=0; j--)
	{
		c[j] = z[j] - mu[j] * c[j+1];
		b[j] = (a[j+1] - a[j])/h[j] - h[j]*(c[j+1] + 2.0*c[j])/3.0;
		d[j] = (c[j+1] - c[j]) / (3*h[j]);
	}
}

/*****************************************************************************/

