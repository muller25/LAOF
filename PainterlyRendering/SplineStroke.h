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
    inline void set(int radius, int r, int g, int b, int x0, int y0)
    {
        m_r = r;
        m_g = g;
        m_b = b;
        m_radius = radius;

        //加到起点，为了计算 brush中每bristle颜色；
        m_start_point = Point(x0, y0);
    }

    inline int nPoints() const{return m_points.size();}
    inline Point getStartPoint(){return m_start_point;}
    inline const Point getStartPoint() const{return m_start_point;}

    inline int nRadius() const{return m_radius;}
    inline int getColorR() const{return m_r;}
    inline int getColorG() const{return m_g;}
    inline int getColorB() const{return m_b;}
    inline void add(int x, int y){
        m_points.push_back(Point(x, y));
    }

    inline Point& get(int index)
    {
        assert(index >= 0 && index < m_points.size());
        return m_points[index];
    }

    inline const Point& get(int index) const
    {
        assert(index >= 0 && index < m_points.size());
        return m_points[index];
    }

    void cubic_b_spline(Point2d &p, double t) const;

private:
	int m_radius;
	unsigned int m_r, m_g, m_b;
    Point m_start_point;
	vector<Point> m_points;
};

#endif
