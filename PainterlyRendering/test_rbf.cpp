#include "RBF.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <stdio.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2){
        printf("parameters: image\n");
        return -1;
    }

    const double factor = 0.125;
    const int radius = 4;
    
    Mat src, im;
    src = imread(argv[1]);
    resize(src, im, Size(0, 0), factor, factor);
    
    // smooth
    printf("smoothing...\n");
    Mat smooth;
    int ksize = 2 * radius + 1;
    GaussianBlur(im, smooth, Size(ksize, ksize), radius, radius);

    // Mat rbfx, rbfy;
    // vector<Point> centers;
    // RBF::rbf_interpolate(rbfx, rbfy, centers, smooth);

    Mat orient;
    vector<Point> centers;
    RBF::rbf_interpolate(orient, centers, smooth);

    // plot
    Mat rbfres;
    RBF::plot(rbfres, centers, orient, im, 1 / factor);
    imshow("rbf", rbfres);
    waitKey(0);
    
    return 0;
}
