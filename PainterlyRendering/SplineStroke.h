#ifndef _SPLINE_STROKE_H_
#define _SPLINE_STROKE_H_
#pragma once

#include "Queue.h"

typedef struct Point
{
	int x, y;

} Point;

typedef struct SplineStroke
{
	int num_points;
	int s;
	unsigned int r ,g,b;
	int start_point_x,start_point_y;
	Queue * points;
	float sumGradient;

} SplineStroke;

class SplineStrokeService
{
public:
	SplineStrokeService(void){}
	~SplineStrokeService(void){}

	/******************************************************************************/
	/*** Creates an empty spline stroke                                         ***/
	/***                                                                        ***/
	/*** Input: s: radius of control points                                     ***/
	/***        r: red                                                          ***/
	/***        g: green                                                        ***/
	/***        b: blue                                                         ***/
	/*** Returns: a new empty spline stroke                                     ***/
	/******************************************************************************/
	static SplineStroke * spline_stroke_create(int s, int r, int g, int b,int x0,int y0);//加到起点，为了计算 brush中每bristle颜色；


	/******************************************************************************/
	/*** Destroy a spline stroke                                                ***/
	/***                                                                        ***/
	/*** Input: spline stroke:                                                  ***/
	/******************************************************************************/
	static void spline_stroke_destroy(SplineStroke * spline_stroke);


	/******************************************************************************/
	/*** Add Control Point to Cubic Spline Stroke at (x, y)                     ***/
	/******************************************************************************/
	static void spline_stroke_add(SplineStroke * spline_stroke, int x, int y);


	/******************************************************************************/
	/*** Get Control Point from Cubic Spline Stroke                             ***/
	/******************************************************************************/
	static Point * spline_stroke_get(SplineStroke * s, int index);
};
#endif
