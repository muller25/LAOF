#include "SplineStroke.h"
/******************************************************************************/
/*** C Headers                                                              ***/
/******************************************************************************/
#include <stdlib.h>

/******************************************************************************/
/*** Create Empty Cubic Spline Stroke                                       ***/
/******************************************************************************/
SplineStroke * SplineStrokeService::spline_stroke_create(int s, int r, int g, int b,int x0,int y0)
{
	SplineStroke * spline = (SplineStroke *) malloc(sizeof(SplineStroke));
	Queue * points = QueueService::queue_create();
	spline->r = r;
	spline->g = g;
	spline->b = b;
	spline->s = s;
	spline->start_point_x = x0;
	spline->start_point_y = y0;
	spline->num_points = 0;
	spline->points = points;
	spline->sumGradient = 0.0;
	return spline;
}


/******************************************************************************/
/*** Destroy Cubic Spline Stroke                                            ***/
/******************************************************************************/
void SplineStrokeService::spline_stroke_destroy(SplineStroke * s)
{
	while (s->num_points > 0)
	{
		Point * p = (Point*)QueueService::queue_remove(s->points);
		free(p);
		s->num_points--;
	}

	free(s->points);
	free(s);
}


/******************************************************************************/
/*** Add Control Point to Cubic Spline Stroke                               ***/
/******************************************************************************/
void SplineStrokeService::spline_stroke_add(SplineStroke * spline_stroke, int x, int y)
{
	Point * point = (Point *) malloc(sizeof(Point));
	point->x = x;
	point->y = y;

//	#pragma omp critical(queue_add)//这里并行修改
	{
		QueueService::queue_add(spline_stroke->points, point);
	}
	spline_stroke->num_points++;
}

/******************************************************************************/
/*** Get Control Point from Cubic Spline Stroke                             ***/
/******************************************************************************/
Point * SplineStrokeService::spline_stroke_get(SplineStroke * s, int index)
{
	if (index < 0 || index >= s->num_points) 
		return NULL;

	return (Point *) QueueService::queue_get(s->points, index);
}
