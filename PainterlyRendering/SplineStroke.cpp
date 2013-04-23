#include "SplineStroke.h"

void SplineStroke::cubic_b_spline(Point2d &p, double t) const
{
	double x0, y0, x1, y1, x2, y2, x3, y3;
    double tt, tt2, tt3, omtt, omtt3;
	int index, npoints = m_points.size();
    
	if (npoints == 0) return;
	if (t > 1) t = 1;

	index = (int)(t * (npoints - 3.0));
	tt = (double) (t - (double)index / ((double)npoints - 3.0)) * (npoints - 3.0);

	if (index < 0) index = 0;
	if (index > npoints - 4) index = npoints - 4;

    Point p0 = m_points[0];
    if (npoints == 1){
		p.x = (double)p0.x + (0.5-t) * m_radius;
		p.y = (double)p0.y;
        return;
    }

    Point p1 = m_points[1];
    if (npoints == 2){
        p.x = (double)p0.x * (1.0 - t) + (p1.x * t);
        p.y = (double)p0.y * (1.0 - t) + (p1.y * t);
        return;
    }
    
    Point p2 = m_points[2];
    if (npoints == 3){
		int ax = p0.x - 2 * p1.x + p2.x;
		int bx = 2 * (p1.x - p0.x);
		int ay = p0.y - 2 * p1.y + p2.y;
		int by = 2 * (p1.y - p0.y);
		p.x = ax * t * t + bx * t + p0.x;
		p.y = ay * t * t + by * t + p0.y;
        return;
    }

    Point p3;
    p0 = m_points[index];
    x0 = p0.x, y0 = p0.y;
    p1 = m_points[index + 1];
    x1 = p1.x, y1 = p1.y;
    p2 = m_points[index + 2];
    x2 = p2.x, y2 = p2.y;
    p3 = m_points[index + 3];
    x3 = p3.x, y3 = p3.y;
	tt2 = tt * tt;
	tt3 = tt * tt2;
	omtt = 1.0 - tt;
	omtt3 = omtt * omtt * omtt;
	p.x = (omtt3 * x0 + (3 * tt3 - 6 * tt2 + 4) * x1 + (-3 * tt3 + 3 * tt2 + 3 * tt + 1) * x2 + tt3 * x3) / 6.0; 
	p.y = (omtt3 * y0 + (3 * tt3 - 6 * tt2 + 4) * y1 + (-3 * tt3 + 3 * tt2 + 3 * tt + 1) * y2 + tt3 * y3) / 6.0; 
}
