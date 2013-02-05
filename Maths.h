#ifndef _MATHS_H
#define _MATHS_H

#include "Constants.h"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
using namespace std;

class Maths
{
public:
    static Mat dx(const Mat &src);
    static Mat dy(const Mat &src);
    static Mat dxx(const Mat &src);
    static Mat dyy(const Mat &src);
    static Mat dxy(const Mat &src);
    static Mat weighted_laplacian(const Mat &flow, const Mat &weight);
    static Mat weighted_laplacian3D(const Mat &pflow, const Mat &flow, const Mat &nflow,
                                    const Mat &weight, const Mat &nweight);
    static int sor_solver(const Mat &A, const Mat &b, Mat &x,
                          int nIters=100, double w=1.8);
};

inline bool match(const Mat& m1, const Mat &m2, bool depth=false);
inline int get_step(const int depth);

#endif
