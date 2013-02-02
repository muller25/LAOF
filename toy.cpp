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

    Mat src(3, 3, CV_8UC1);
    RNG rng(0);
    rng.fill(src, RNG::UNIFORM, -10, 10);
    
    cout << src << endl;
    
    Mat dst;
    dst = Maths::dx(src);
    cout << dst << endl;

    dst = Maths::dy(src);
    cout << dst << endl;

    dst = Maths::dxx(src);
    cout << dst << endl;

    dst = Maths::dyy(src);
    cout << dst << endl;

    dst = Maths::dxy(src);
    cout << dst << endl;
    
    // int b_init[2][1] = {{2}, {5}};
        
    // Mat A = (Mat_<double>(2, 2) << 1, 1, 2, 3);
    // Mat x = Mat(2, 1, CV_64FC1, Scalar(0));
    // Mat b = (Mat_<double>(2, 1) << 2, 5);

    // Maths::sor_solver(A, b, x, 100);
    // cout << "x = " << x.at<double>(0, 0) << ", y = " << x.at<double>(0, 1) << endl;
    
    return 0;
}
