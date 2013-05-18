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
    inline void set(int radius, Scalar color, int x, int y)
    {
        set(radius, color.val[2], color.val[1], color.val[0], x, y);
    }
    
    inline void set(int radius, int r, int g, int b, int x0, int y0)
    {
        m_r = r;
        m_g = g;
        m_b = b;
        m_radius = radius;
        m_alpha = 1.0;

        //加到起点，为了计算 brush中每bristle颜色；
        m_points.push_back(Point(x0, y0));
    }

    inline void add(int x, int y){m_points.push_back(Point(x, y));}
    
    inline int nPoints() const{return m_points.size();}
    inline Point getStartPoint(){return m_points[0];}
    inline const Point getStartPoint() const{return m_points[0];}
    inline int nRadius() const{return m_radius;}
    inline double nAlpha() const{return m_alpha;}
    inline Scalar color() const{
        Scalar color(m_b, m_g, m_r);
        return color;
    }
    
    inline int ColorR() const{return m_r;}
    inline int ColorG() const{return m_g;}
    inline int ColorB() const{return m_b;}
    inline void setAlpha(double alpha){m_alpha = alpha;}
    inline bool isTransparent() const{return fabs(m_alpha) < 1e-6;}

    inline void fadeOut(double step){
        m_alpha -= step;
        if (m_alpha < 0) m_alpha = 0;
        if (m_alpha > 1) m_alpha = 1;
    }
    
    inline Point& get(int index)
    {
        assert(index >= 0 && index < (int)m_points.size());
        return m_points[index];
    }

    inline const Point& get(int index) const
    {
        assert(index >= 0 && index < (int)m_points.size());
        return m_points[index];
    }

    inline void setControlPoint(int index, Point &p)
    {
        assert(index >= 0 && index < (int)m_points.size());
        m_points[index] = p;
    }

    inline void resize(int size){m_points.resize(size);}
    
    void cubic_b_spline(Point2d &p, double t) const;

private:
	int m_radius;
    double m_alpha;
	unsigned int m_r, m_g, m_b;
	vector<Point> m_points;
};

#endif
