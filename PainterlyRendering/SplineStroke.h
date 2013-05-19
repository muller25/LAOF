#ifndef _SplineStroke_H
#define _SplineStroke_H

#include <cv.h>
using cv::Vec3b;
using cv::Point;

#include <vector>
using std::vector;

class SplineStroke
{
public:
    SplineStroke(){}

    inline void set(int radius, double angle, Vec3b color)
    {
        m_color = color;
        m_radius = radius;
        m_angle = angle;
        m_alpha = 1;
        m_points.clear();
    }
    
    inline void set(int radius, double angle, int r, int g, int b)
    {
        set(radius, angle, Vec3b(b, g, r));
    }

    inline void add(int x, int y){m_points.push_back(Point(x, y));}
    
    inline int nPoints() const{return m_points.size();}
    inline Point getStartPoint() const{return m_points[0];}
    inline int getRadius() const{return m_radius;}
    inline double getAngle() const{return m_angle;}
    inline double setAngle(double angle){m_angle = angle;}
    inline Vec3b getColor() const{return color;}
    inline double getAlpha() const{return m_alpha;}
    inline void setAlpha(double alpha){m_alpha = alpha;}
    inline bool isTransparent() const{return fabs(m_alpha) < 1e-6;}
    inline void changeOpacity(double step){
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
    double m_alpha, m_angle;
    Vec3b m_color; // BGR
	vector<Point> m_points;
};

#endif
