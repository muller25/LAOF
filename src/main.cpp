#include "OpticalFlow.h"
#include "Flow2Color.h"
#include "ImageIO.h"

#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char *argv[])
{
    // load color image
    const char *im1Name = "../car1.jpg";
    const char *im2Name = "../car2.jpg";
    
    DImage im1, im2, u, v;
    imread(im1, im1Name);
    imread(im2, im2Name);
    
    // init optical flow parameters
    const double a_s = 0.012;
    const double ratio = 0.75;
    const int minWidth = 5;
    const int nOutIter = 6;
    const int nInIter = 1;
    const int nSORIter = 30;

    cout << "start computing optical flow..." << endl;
    OpticalFlow of;
    of.c2fFlow(u, v, im1, im2, a_s, ratio, minWidth,
               nOutIter, nInIter, nSORIter);

    DImage warp;
    of.warpImage(warp, im1, im2, u, v);
    imwrite("warp.jpg", warp);
    
//    imshow("warp image", warp);

    UCImage flowImg, idxImg;
    flow2color(flowImg, idxImg, u, v);
    imwrite("flow.jpg", flowImg);
    imwrite("idxImg.jpg", idxImg);
    
    imshow("flow image", flowImg);
    imshow("unknown flow index image", idxImg);
    imwait(0);
    
    return 0;
}
