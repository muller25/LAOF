#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
using namespace cv;

#include <iostream>
using namespace std;

#include "Maths.h"

int main(int argc, char *argv[])
{
    // Mat im1 = imread("car1.jpg", CV_LOAD_IMAGE_COLOR);
    // Mat im2 = imread("car2.jpg", CV_LOAD_IMAGE_COLOR);
    RNG rng(0x7fffffff);
        
    Mat pflow(3, 3, CV_32SC1), flow(3, 3, CV_32SC1), nflow(3, 3, CV_32SC1);
    Mat weight(3, 3, CV_32SC1), nweight(3, 3, CV_32SC1);

    rng.fill(pflow, RNG::UNIFORM, -1, 1);
    rng.fill(flow, RNG::UNIFORM, -10, 10);
    rng.fill(nflow, RNG::UNIFORM, -1, 1);
    rng.fill(weight, RNG::UNIFORM, -10, 10);
    rng.fill(nweight, RNG::UNIFORM, -1, 1);
    
    cout << "***** params *****" << endl;
//    cout << pflow << endl;
    cout << flow << endl;
//    cout << nflow << endl;
    cout << weight << endl;
//    cout << nweight << endl;
    
    Mat lap = Maths::weighted_laplacian<int, int>(flow, weight);
    cout << lap << endl;

    // cout << "***** dx *****" << endl;
    // dx = Maths::dx(src);
    // cout << dx << endl;

    // cout << "***** dy *****" << endl;
    // dy = Maths::dy(src);
    // cout << dy << endl;

    // cout << "***** dxx *****" << endl;
    // dxx = Maths::dxx(src);
    // cout << dxx << endl;

    // cout << "***** dyy *****" << endl;
    // dyy = Maths::dyy(src);
    // cout << dyy << endl;

    // cout << "***** dxy *****" << endl;
    // dxy = Maths::dxy(src);
    // cout << dxy << endl;

    // cout << "***** dxdy *****" << endl;
    // multiply(dx, dy, dxdy, 1);
    // cout << dxdy << endl;
    
    // int b_init[2][1] = {{2}, {5}};
        
    // Mat A = (Mat_<double>(2, 2) << 1, 1, 2, 3);
    // Mat x = Mat(2, 1, CV_64FC1, Scalar(0));
    // Mat b = (Mat_<double>(2, 1) << 2, 5);

    // Maths::sor_solver(A, b, x, 100);
    // cout << "x = " << x.at<double>(0, 0) << ", y = " << x.at<double>(0, 1) << endl;
    
    return 0;
}
