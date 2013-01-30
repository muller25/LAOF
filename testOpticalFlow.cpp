#include <opencv2/core/core.hpp>
#include <iostream>

#include "OpticalFlow.h"
using namespace std;
using namespace cv;

int main(int argc, char *argv[])
{
    int b_init[2][1] = {{2}, {5}};
        
    Mat A = (Mat_<double>(2, 2) << 1, 1, 2, 3);
    Mat x = Mat::zeros(2, 1, CV_64FC1);
    Mat b = (Mat_<double>(2, 1) << 2, 5);
    
    OpticalFlow::SOR_solver(A, b, x, 100);
    cout << "x = " << x.at<double>(0, 0) << ", y = " << x.at<double>(0, 1) << endl;
    
    return 0;
}
