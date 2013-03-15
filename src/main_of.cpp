#include "OpticalFlow.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"
#include <vector>

int main(int argc, char *argv[])
{
    const char *inImg = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *outFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    const char *outImg = "/home/iaml/Projects/exp/lena/out/%s%03d.jpg";

    // init optical flow parameters
    const double as = 0.026;
    const double ap = 0.012;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nBiIter = 5;//14;
    const int nIRLSIter = 1;
    const int nSORIter = 10;//30;
  
    char buf[256];
    std::vector<DImage> im(2);
    DImage u1, v1, u2, v2, warp;
    UCImage flowImg, mask1, mask2;
    int cur, next, width, height;
    OpticalFlow of;
    
    sprintf(buf, inImg, 0);
    imread(im[0], buf);
    width = im[0].nWidth(), height = im[0].nHeight();
    mask1.create(width, height, 1, 1);
    mask2.create(width, height, 1, 1);
    cur = 0;
    next = 1;
    for (int i = 1; i < 10; ++i)
    {
        sprintf(buf, inImg, i);
        imread(im[next], buf);

        printf("running biDir optical flow from im%d to im%d...\n", i-1, i);
        of.biC2FFlow(u1, v1, u2, v2, im[cur], im[next], mask1, mask2,
                     as, ap, ratio, minWidth, nBiIter, nIRLSIter, nSORIter);
        
        sprintf(buf, outFile, "u", i-1);
        imwritef(buf, u1);
        sprintf(buf, outFile, "v", i-1);
        imwritef(buf, v1);

        flow2color(flowImg, u1, v1);
        sprintf(buf, outImg, "flow", i-1);
        imwrite(buf, flowImg);

        sprintf(buf, outFile, "ur", i);
        imwritef(buf, u2);
        sprintf(buf, outFile, "vr", i);
        imwritef(buf, v2);

        flow2color(flowImg, u2, v2);
        sprintf(buf, outImg, "rflow", i);
        imwrite(buf, flowImg);

        warpImage(warp, im[cur], im[next], u1, v1);
        sprintf(buf, outImg, "warp", i-1);
        imwrite(buf, warp);

        warpImage(warp, im[next], im[cur], u2, v2);
        sprintf(buf, outImg, "rwarp", i);
        imwrite(buf, warp);
        
        // time window move forward
        cur = 1 - cur;
        next = 1 - next;
    }
    
    return 0;
}
