#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>

#include "Maths.h"
using namespace std;
using namespace cv;

int main(int argc, char *argv[])
{
    // Mat im1 = imread("car1.jpg", CV_LOAD_IMAGE_COLOR);
    // Mat im2 = imread("car2.jpg", CV_LOAD_IMAGE_COLOR);

    Mat src(5, 5, CV_8UC1), im0(5, 5, CV_8UC1), im2(5, 5, CV_8UC1);
    RNG rng(12345);
    rng.fill(src, RNG::UNIFORM, -10, 10);
    rng.fill(im0, RNG::UNIFORM, 5, 10);
    rng.fill(im2, RNG::UNIFORM, 0, 5);
    
    cout << "***** src *****" << endl;
    cout << im0 << endl;
    cout << src << endl;
    cout << im2 << endl;

    Mat dz, dxdy, dx, dy, dxx, dyy, dxy, laplace;
    
    laplace = Maths::laplace3D(im0, src, im2);
    cout << laplace << endl;
    
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
