#include "OpticalFlow.h"
#include "Flow2Color.h"
#include "ImageIO.h"

#include <iostream>
#include <ctime>
using namespace std;
/*
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
*/
int main(int argc, char *argv[])
{
    const char *inName = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *outFlow = "/home/iaml/Projects/exp/lena/out/%s%03d.jpg";
    const char *outWarp = "/home/iaml/Projects/exp/lena/out/%swarp%03d.jpg";
    const char *outFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    
    // init optical flow parameters
    const double as = 0.026;
    const double ap = 0.012;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nOutIter = 5;//14;
    const int nInIter = 1;
    const int nSORIter = 10;//30;

    char buf[256];
    DImage im[3], u[3], v[3], ur[3], vr[3], mask[3], warp;
    UCImage flowImg, idxImg;
    OpticalFlow of;
    
    sprintf(buf, inName, 0);
    imread(im[0], buf);
    sprintf(buf, inName, 1);
    imread(im[1], buf);
    sprintf(buf, inName, 2);
    imread(im[2], buf);
    
    int width = im[0].nWidth(), height = im[0].nHeight();
    for (int i = 0; i < 3; ++i)
        mask[i].create(width, height, 1, 1);

    for (int i = 0; i < 2; ++i)
    {
        of.biC2FFlow(u[i], v[i], ur[i], vr[i], im[i], im[i+1], mask[i], mask[i+1],
                     as, ap, ratio, minWidth,
                     nOutIter, nInIter, nSORIter);

        of.warpImage(warp, im[i], im[i+1], u[i], v[i]);
        sprintf(buf, outWarp, "", i);
        imwrite(buf, warp);

        of.warpImage(warp, im[i+1], im[i], ur[i], vr[i]);
        sprintf(buf, outWarp, "r", i);
        imwrite(buf, warp);

        flow2color(flowImg, idxImg, u[i], v[i]);
        sprintf(buf, outFlow, "", i);
        imwrite(buf, flowImg);
        
        flow2color(flowImg, idxImg, ur[i], vr[i]);
        sprintf(buf, outFlow, "r", i);
        imwrite(buf, flowImg);
    }
        
    im[0] = im[1];
    im[1] = im[2];
    mask[0] = mask[1];
    mask[1] = mask[2];
    for (int i = 3; i < 4; ++i)
    {
        sprintf(buf, inName, i);
        imread(im[2], buf);
        
        of.biC2FFlow(u[2], v[2], ur[2], vr[2], im[1], im[2], mask[1], mask[2],
                     as, ap, ratio, minWidth,
                     nOutIter, nInIter, nSORIter);

        of.warpImage(warp, im[1], im[2], u[2], v[2]);
        sprintf(buf, outWarp, "", i);
        imwrite(buf, warp);

        of.warpImage(warp, im[2], im[1], ur[2], vr[2]);
        sprintf(buf, outWarp, "", i);
        imwrite(buf, warp);

        flow2color(flowImg, idxImg, u[2], v[2]);
        sprintf(buf, outFlow, "", i);
        imwrite(buf, flowImg);

        flow2color(flowImg, idxImg, ur[2], vr[2]);
        sprintf(buf, outFlow, "r", i);
        imwrite(buf, flowImg);

        of.stFlow(u[1], v[1], ur[1], vr[1], im[0], im[1], mask[0], mask[1],
                  u[0], v[0], u[2], v[2], ur[0], vr[0], ur[2], vr[2],
                  as, ap, nOutIter, nInIter+2, nSORIter);
        
        of.warpImage(warp, im[0], im[1], u[1], v[1]);
        sprintf(buf, outWarp, "t", i);
        imwrite(buf, warp);

        of.warpImage(warp, im[1], im[0], ur[1], vr[1]);
        sprintf(buf, outWarp, "tr", i);
        imwrite(buf, warp);

        flow2color(flowImg, idxImg, u[1], v[1]);
        sprintf(buf, outFlow, "t", i);
        imwrite(buf, flowImg);

        flow2color(flowImg, idxImg, ur[1], vr[1]);
        sprintf(buf, outFlow, "tr", i);
        imwrite(buf, flowImg);

        // move time window forward
        for (int j = 0; j < 2; j++)
        {
            im[j] = im[j+1];
            u[j] = u[j+1];
            v[j] = v[j+1];
            ur[j] = ur[j+1];
            vr[j] = vr[j+1];
            mask[j] = mask[j+1];
        }
    }
    
    return 0;
}
