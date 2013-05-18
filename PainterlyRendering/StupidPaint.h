#ifndef _StupidPaint_H
#define _StupidPaint_H

#include "Brush.h"
#include <cv.h>
using cv::Mat;
using cv::Vec3f;

#include <vector>
using std::vector;

void rotate(Mat &dst, const Mat &src, double radian);

class StupidPaint
{
public:
    static void loadTexture();
    static void generate_strokes(vector<Brush> &strokes, const Mat &src, int radius);
    static void strokes_placement(Mat &dst, const Mat &src, const vector<Brush> &strokes);

    static Vec3b extract_color(const Mat &src, const Mat &mask, int x, int y);
    static void stroke_orient(Mat &orient, const Mat &src);

    static void paint(Mat &dst, const Mat &src);
private:
    static Mat m_height_map, m_alpha_map;
};

#endif
