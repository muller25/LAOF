#include "CubicBSpline.h"
#include<cv.h>
#include"PainterlyService.h"
#define POINT(x) ((Point *) QueueService::queue_get(stroke->points, x))

void CubicBSpline::cubic_b_spline(const SplineStroke * stroke, double t, double * x, double * y)
{
	double x0, y0;
	double x1, y1;
	double x2, y2;
	double x3, y3;
    double tt, tt2, tt3, omtt, omtt3;
	int index;
//	static double ts;
	if (stroke->num_points == 0) return;
	if(t>1) t=1;

	index = (int) (t * (stroke->num_points - 3.0));
	tt = (double) (t - (float) index / ((float) stroke->num_points - 3.0)) * (stroke->num_points - 3.0); 
	//tt=t*t;
	//index -= 3;

	if (index < 0) index = 0;
	//if (index == stroke->num_points - 4||index == stroke->num_points)
	//	return;
	//	ts=t;
	if (index > stroke->num_points - 4)
	{
		index = stroke->num_points - 4;
		//t=(t-ts)/(1.0-ts);
	}

	if (stroke->num_points == 1)
	{
		*x = ((Point *) QueueService::queue_get(stroke->points, 0))->x+(0.5-t)*stroke->s;
		*y = ((Point *) QueueService::queue_get(stroke->points, 0))->y;
		return;
	}

	if (stroke->num_points == 2)
	{
		*x = (((Point *) QueueService::queue_get(stroke->points, 0))->x * (1.0 - t) + ((Point *) QueueService::queue_get(stroke->points, 1))->x * t);
		*y = (((Point *) QueueService::queue_get(stroke->points, 0))->y * (1.0 - t) + ((Point *) QueueService::queue_get(stroke->points, 1))->y * t);
		return;
	}

	if (stroke->num_points == 3)
	{
		x0 = ((Point *) QueueService::queue_get(stroke->points, 0))->x;
		y0 = ((Point *) QueueService::queue_get(stroke->points, 0))->y;
		x1 = ((Point *) QueueService::queue_get(stroke->points, 1))->x;
		y1 = ((Point *) QueueService::queue_get(stroke->points, 1))->y;
		x2 = ((Point *) QueueService::queue_get(stroke->points, 1))->x;
		y2 = ((Point *) QueueService::queue_get(stroke->points, 1))->y;
		x3 = ((Point *) QueueService::queue_get(stroke->points, 2))->x;
		y3 = ((Point *) QueueService::queue_get(stroke->points, 2))->y;

		int ax=POINT(0)->x-2 * POINT(1)->x+POINT(2)->x;
		int bx=2*(POINT(1)->x-POINT(0)->x);
		int ay=POINT(0)->y-2 * POINT(1)->y+POINT(2)->y;
		int by=2*(POINT(1)->y-POINT(0)->y);
		*x=ax*t*t+bx*t+POINT(0)->x;
		*y=ay*t*t+by*t+POINT(0)->y;
		return;
	}
	else
	{
		x0 = ((Point *) QueueService::queue_get(stroke->points, index))->x;
		y0 = ((Point *) QueueService::queue_get(stroke->points, index))->y;
		x1 = ((Point *) QueueService::queue_get(stroke->points, index + 1))->x;
		y1 = ((Point *) QueueService::queue_get(stroke->points, index + 1))->y;
		x2 = ((Point *) QueueService::queue_get(stroke->points, index + 2))->x;
		y2 = ((Point *) QueueService::queue_get(stroke->points, index + 2))->y;
		x3 = ((Point *) QueueService::queue_get(stroke->points, index + 3))->x;
		y3 = ((Point *) QueueService::queue_get(stroke->points, index + 3))->y;
	}

	tt2 = tt * tt;
	tt3 = tt * tt2;
	omtt = 1.0 - tt;
	omtt3 = omtt * omtt * omtt;

	//cvCircle(PainterlyService::map_image, cvPoint(x0,y0), 2, CV_RGB(255,0,0), -1);
	//cvCircle(PainterlyService::map_image, cvPoint(x1,y1), 2, CV_RGB(255,0,0), -1);
	//cvCircle(PainterlyService::map_image, cvPoint(x2,y2), 2, CV_RGB(255,0,0), -1);
	//cvCircle(PainterlyService::map_image, cvPoint(x3,y3), 2, CV_RGB(255,0,0), -1);

	*x = (omtt3 * x0 + (3 * tt3 - 6 * tt2 + 4) * x1 + (-3 * tt3 + 3 * tt2 + 3 * tt + 1) * x2 + tt3 * x3) / 6.0; 
	*y = (omtt3 * y0 + (3 * tt3 - 6 * tt2 + 4) * y1 + (-3 * tt3 + 3 * tt2 + 3 * tt + 1) * y2 + tt3 * y3) / 6.0; 
}

void CubicBSpline::cubic_b_spline_new(SplineStroke * stroke, double t, double * x, double * y)
{

#define POINT(x) ((Point *) QueueService::queue_get(stroke->points, x))

	double x0, y0;
	double x1, y1;
	double x2, y2;
	double x3, y3;
//	double tt, tt2, tt3, omtt, omtt3;
	int index;

	

	if (stroke->num_points == 0)    
		return;

	index=t*(stroke->num_points-1);
	if(index>stroke->num_points-4) 
		index=stroke->num_points-4;

	if (stroke->num_points == 1)
	{
		*x = POINT(0)->x;
		*y = POINT(0)->y;
		return;
	}

	if (stroke->num_points == 2)
	{
		*x = (POINT(0)->x * t  + POINT(1)->x * (1.0 - t));
		*y = (POINT(0)->y * t  + POINT(1)->y * (1.0 - t));
		return;
	}
	if (stroke->num_points == 3)
	{
		int ax=POINT(0)->x-2 * POINT(1)->x+POINT(2)->x;
		int bx=2*(POINT(1)->x-POINT(0)->x);
		int ay=POINT(0)->y-2 * POINT(1)->y+POINT(2)->y;
		int by=2*(POINT(1)->y-POINT(0)->y);
		*x=ax*t*t+bx*t+POINT(0)->x;
		*y=ay*t*t+by*t+POINT(0)->y;
		return;
	}
	else
	{
		x0 = ((Point *) QueueService::queue_get(stroke->points, index))->x;
		y0 = ((Point *) QueueService::queue_get(stroke->points, index))->y;
		x1 = ((Point *) QueueService::queue_get(stroke->points, index + 1))->x;
		y1 = ((Point *) QueueService::queue_get(stroke->points, index + 1))->y;
		x2 = ((Point *) QueueService::queue_get(stroke->points, index + 2))->x;
		y2 = ((Point *) QueueService::queue_get(stroke->points, index + 2))->y;
		x3 = ((Point *) QueueService::queue_get(stroke->points, index + 3))->x;
		y3 = ((Point *) QueueService::queue_get(stroke->points, index + 3))->y;

		float   ax, bx, cx;
		float   ay, by, cy;
		float   tSquared, tCubed;
 
		/*計算多項式係數*/
 
		cx = 3.0 * (x1 - x0);
		bx = 3.0 * (x2 - x1) - cx;
		ax = x3 - x0 - cx - bx;
 
		cy = 3.0 * (y1 - y0);
		by = 3.0 * (y2 - y1) - cy;
		ay = y3 - y0 - cy - by;
 
		/*計算位於參數值t的曲線點*/
 
		tSquared = t * t;
		tCubed = tSquared * t;
 
		*x = (ax * tCubed) + (bx * tSquared) + (cx * t) + x0;
		*y = (ay * tCubed) + (by * tSquared) + (cy * t) + y0;
	}
}
