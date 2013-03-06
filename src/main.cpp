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
    
    DImage im1, im2, u, v, ur, vr;
    imread(im1, im1Name);
    imread(im2, im2Name);

    int width = im1.nWidth(), height = im1.nHeight();
    DImage mask1(width, height, 1, 1), mask2(width, height, 1, 1);

    // init optical flow parameters
    const double as = 0.026;
    const double ap = 0.012;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nOutIter = 5;//14;
    const int nInIter = 1;
    const int nSORIter = 10;//30;
    
    cout << "start computing optical flow..." << endl;
    OpticalFlow of;
    // of.c2fFlow(u, v, im1, im2, as, ratio, minWidth,
    //            nOutIter, nInIter, nSORIter);

    of.biC2FFlow(u, v, ur, vr, im1, im2, mask1, mask2,
                 as, ap, ratio, minWidth,
                 nOutIter, nInIter, nSORIter);
    
    imwritef("u.yml", u);
    imwritef("v.yml", v);
    imwritef("ur.yml", ur);
    imwritef("vr.yml", vr);
    
    DImage warp;
    of.warpImage(warp, im1, im2, u, v);
    imwrite("warp.jpg", warp);

    of.warpImage(warp, im2, im1, ur, vr);
    imwrite("warpr.jpg", warp);

    UCImage flowImg, rflowImg, idxImg;
    flow2color(flowImg, idxImg, u, v);
    imwrite("flow.jpg", flowImg);
    flow2color(rflowImg, idxImg, ur, vr);
    imwrite("rflow.jpg", rflowImg);
    imwrite("ridxImg.jpg", idxImg);
    
    imshow("flow image", flowImg);
    imshow("rflow image", rflowImg);

//    imshow("unknown flow index image", idxImg);
    imwait(0);
    
    return 0;
}
