#ifndef _CUBIC_B_SPLINE_H
#define _CUBIC_B_SPLINE_H

#include "SplineStroke.h"

class CubicBSpline
{
public:
	CubicBSpline(void){}
	~CubicBSpline(void){}

	static void cubic_b_spline(const SplineStroke &stroke, double t, double *x, double *y);
};

#endif
