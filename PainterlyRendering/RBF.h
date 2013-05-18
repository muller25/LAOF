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
    static void rbf_center(vector<Point> &centers, const Mat &gm,
                           double thres = 100, int wsize = 5, int range = 10);
    static void rbf_solver(Mat &rbfres, const vector<Point> &centers, const Mat &gm);
    static void rbf(Mat &res, const Mat &src, int radius);
    static void orientation(Mat &field, const Mat &src, int raidus);
    static void orientation(Mat &field, const Mat &rbfres);

    static double phi(double x){return fabs(x);}
    static void plot_gradient(Mat &show, const Mat &gx, const Mat &gy,
                              const Mat &im, double factor);
    static double factor;
    static int wsize;
    static int range;
};

#endif
