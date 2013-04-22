#ifndef _CUBIC_B_SPLINE_H_
#define _CUBIC_B_SPLINE_H_
#pragma once

#include "SplineStroke.h"

class CubicBSpline
{
public:
	CubicBSpline(void){}
	~CubicBSpline(void){}

	static void cubic_b_spline(const SplineStroke * stroke, double t, double * x, double * y);

	static void cubic_b_spline_new(SplineStroke * stroke, double t, double * x, double * y);
};

#endif
