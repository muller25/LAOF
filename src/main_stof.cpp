#include "OpticalFlow.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"
#include <vector>

char *inPattern, *outDir;
int frameStart, frameEnd, numOfSegs;

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        printf("./em input_pattern output_dir start end #_of_seg\n");
        return 1;
    }
    
    inPattern = argv[1];
    outDir = argv[2];
    frameStart = atoi(argv[3]);
    frameEnd = atoi(argv[4]);
    numOfSegs = atoi(argv[5]);

    // init optical flow parameters
    const double as = 0.026;
    const double ap = 0.012;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nOutIter = 5;//14;
    const int nIRLSIter = 1;
    const int nSORIter = 10;//30;
  
    char buf[256], outFile[128], outImg[128];
    std::vector<DImage> im(2);
    DImage u, v, warp, mask1, mask2;
    UCImage flowImg;
    int cur, next, width, height;
    OpticalFlow of;

    memset(outFile, 0, sizeof(outFile));
    memset(outImg, 0, sizeof(outImg));
    strcat(outFile, outDir);
    strcat(outFile, "%s%03d.yml");
    strcat(outImg, outDir);
    strcat(outImg, "%s%03d.jpg");

    sprintf(buf, inPattern, frameStart);
    imread(im[0], buf);

    // to save time
    DImage tmp;
    imresize(tmp, im[0], 0.5);
    tmp.copyTo(im[0]);
    sprintf(buf, outImg, "im", frameStart);
    imwrite(buf, tmp);
    
    width = im[0].nWidth(), height = im[0].nHeight();
    mask1.create(width, height, 1, 1);
    mask2.create(width, height, 1, 1);

    cur = 0, next = 1;
    for (int i = frameStart+1; i <= frameEnd; ++i)
    {
        sprintf(buf, inPattern, i);
        imread(im[next], buf);

        // to save time
        imresize(tmp, im[next], 0.5);
        tmp.copyTo(im[next]);
        sprintf(buf, outImg, "im", i);
        imwrite(buf, tmp);

        printf("running optical flow from im%d to im%d...\n", i-1, i);
        of.stC2FFlow(u, v, im[cur], im[next], mask1, mask2,
                     as, ap, ratio, minWidth, nOutIter, nIRLSIter, nSORIter);
        
        sprintf(buf, outFile, "u", i-1);
        imwritef(buf, u);
        sprintf(buf, outFile, "v", i-1);
        imwritef(buf, v);

        flow2color(flowImg, u, v);
        sprintf(buf, outImg, "flow", i-1);
        imwrite(buf, flowImg);

        warpImage(warp, im[cur], im[next], u, v);
        sprintf(buf, outImg, "warp", i-1);
        imwrite(buf, warp);
        
        // time window move forward
        cur = 1 - cur;
        next = 1 - next;
    }
    
    return 0;
}
