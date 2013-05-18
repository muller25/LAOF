#ifndef _Brush_H
#define _Brush_H

#include <cv.h>
using cv::Size;
using cv::Point;
using cv::Vec3b;

class Brush
{
public:
    Brush(){}
    
    Brush(Point center, Vec3b color, Size size, double angle, double opacity=1)
    {
        set(center, color, size, angle, opacity);
    }
        
    void set(Point center, Vec3b color, Size size, double angle, double opacity)
    {
        m_center = center;
        m_color = color;
        m_size = size;
        m_angle = angle;
        m_opacity = opacity;
    }
    
    void setColor(Vec3b color){m_color = color;}
    Vec3b getColor() const{return m_color;}
    void setCenter(Point p){m_center = p;}
    Point getCenter() const{return m_center;}
    void setAngle(double angle){m_angle = angle;}
    double getAngle() const{return m_angle;}
    void setOpacity(double opacity){m_opacity = opacity;}
    double getOpacity() const{return m_opacity;}
    void setSize(Size size){m_size = size;}
    Size getSize() const{return m_size;}
    
private:
    Vec3b m_color;
    Point m_center;
    Size m_size;
    double m_angle, m_opacity;
};

#endif
