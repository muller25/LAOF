#include "OpticalFlow.h"
#include "Flow2Color.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char *argv[])
{
/*    // load color image
    const char *im1Name = "../car1.jpg";
    const char *im2Name = "../car2.jpg";
    const double factor = 1./255;
    
    Mat im1, im2;
    Mat tmp = imread(im1Name, CV_LOAD_IMAGE_COLOR);
    tmp.convertTo(im1, CV_64FC3, factor); // im2double
    
    tmp = imread(im2Name, CV_LOAD_IMAGE_COLOR);
    tmp.convertTo(im2, CV_64FC3, factor); // im2double
*/
    int nrows = 3, ncols = 3;
    Mat im1(nrows, ncols, CV_64FC3), im2(nrows, ncols, CV_64FC3);
    RNG rng(time(NULL));
    rng.fill(im1, RNG::UNIFORM, 0., 1.);
    rng.fill(im2, RNG::UNIFORM, 0., 1.);
    
    for (int r = 0; r < nrows-1; r++)
        im1.row(r).copyTo(im2.row(r+1));
    
    // init optical flow parameters
    const int nOutIter = 7;
    const int nInIter = 7;
    const int nSORIter = 30;
    const double a_s = 0.012;
    int rows = im1.rows, cols = im1.cols, channels = im1.channels();
    Mat u = Mat::zeros(rows, cols, CV_64F);
    Mat v = Mat::zeros(rows, cols, CV_64F);
    Mat warp(rows, cols, CV_64FC(channels));
    OpticalFlow of;

    of.warpImage(im1, im2, u, v, warp);
    of.compute(im1, im2, warp, u, v, a_s, nOutIter, nInIter, nSORIter);

    cout << "im1 " << im1 << endl;
    cout << "im2 " << im2 << endl;
    cout << u << endl;
    cout << v << endl;
    
/*    warp.convertTo(tmp, CV_8UC(channels));
    imshow("warp image", tmp);

    Mat flowImg, idxImg;
    flow2color(u, v, flowImg, idxImg);
    imshow("flow image", flowImg);
    imshow("unknown flow index image", idxImg);
    waitKey(0);
*/  
    return 0;
}
