#ifndef _RBF_H
#define _RBF_H

#include <cv.h>
#include <vector>

using cv::Mat;
using cv::Point;
using std::vector;

typedef double PERCISION;
#define CV_PERCISION CV_64F

class RBF
{
public:
    static void rbf_center(vector<Point> &centers, const Mat &gm);
    static void rbf_solver(Mat &rbfres, const vector<Point> &centers, const Mat &gm);
    static void rbf_interpolate(Mat &rbfx, Mat &rbfy, const Mat &src);
    static void rbf_interpolate(Mat &rbfx, Mat &rbfy, vector<Point> &centers, const Mat &src);
    static void rbf_interpolate(Mat &orient, vector<Point> &centers, const Mat &im);
    static void rbf_interpolate(Mat &orient, const Mat &im);
    
    static double phi(double x){return fabs(x);}
    static void plot(Mat &show, const vector<Point> &centers,
                     const Mat &gx, const Mat &gy, const Mat &im, double scale=8);
    static void plot(Mat &show, const vector<Point> &centers,
                     const Mat &orient, const Mat &im, double factor);
    static void plot(Mat &show, const Mat &orient, const Mat &im, double factor);
};

#endif
