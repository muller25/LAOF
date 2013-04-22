#ifndef _SplineStroke_H
#define _SplineStroke_H

#include <cv.h>
using namespace cv;

#include <vector>
using std::vector;

class SplineStroke
{
public:
    SplineStroke(){}
    inline void set(int s, int r, int g, int b, int x0, int y0);
    inline void add(int x, int y);
    void cubic_b_spline(Point2d &p, double t) const;

    inline Point& get(int index);
    inline const Point& get(int index) const;
    inline const Point SplineStroke::getStartPoint() const;
    inline int nPoints() const;
    inline int SplineStroke::nRadius() const;
    inline int SplineStroke::getColorR() const;
    inline int SplineStroke::getColorG() const;
    inline int SplineStroke::getColorB() const;

private:
	int m_radius;
	unsigned int m_r, m_g, m_b;
    Point m_start_point;
	vector<Point> m_points;
};

#endif
