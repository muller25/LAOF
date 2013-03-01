#include "GaussianPyramid.h"

#include <cv.h>
#include <highgui.h>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    Mat im = imread("../car1.jpg", CV_LOAD_IMAGE_COLOR);
    const char *name = "test %d";
    char buf[256];
    
    std::cout << "constructing pyramid\n";
    
    GaussianPyramid pyr;
    pyr.ConstructPyramid(im, 0.75, 100);

    std::cout << "show pyramid\n";

    for (int i = 0; i < pyr.nLevels(); i++)
    {
        sprintf(buf, name, i);
        std::cout << "image size " << pyr[i].cols << ", " << pyr[i].rows << std::endl;
        imshow(buf, pyr[i]);
    }
    
    waitKey(0);

    std::cout << "constructing pyramid using opencv\n";
    
    GaussianPyramid pyr2;
    pyr2.ConstructPyramid(im, 50);

    std::cout << "show pyramid\n";
    for (int i = 0; i < pyr2.nLevels(); i++)
    {
        sprintf(buf, name, i);
        std::cout << "image size " << pyr2[i].cols << ", " << pyr2[i].rows << std::endl;
        imshow(buf, pyr2[i]);
    }
    
    waitKey(0);

    std::cout << "constructing pyramid using buildPyramid\n";
    vector<Mat> pyr3;
    buildPyramid(im, pyr3, 4, BORDER_REPLICATE);

    std::cout << "show pyramid\n";
    for (size_t i = 0; i < pyr3.size(); i++)
    {
        sprintf(buf, name, i);
        std::cout << "image size " << pyr3[i].cols << ", " << pyr3[i].rows << std::endl;
        imshow(buf, pyr3[i]);
    }
    
    waitKey(0);
    
    return 0;
}
