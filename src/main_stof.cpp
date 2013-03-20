#include "OpticalFlow.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"
#include <vector>

char *inPattern, *outDir;
int frameStart, frameEnd, numOfSegs;

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("./em input_pattern output_dir start end\n");
        return 1;
    }
    
    inPattern = argv[1];
    outDir = argv[2];
    frameStart = atoi(argv[3]);
    frameEnd = atoi(argv[4]);

    // init optical flow parameters
    const double as = 0.026;
    const double at = 0.005;
    const double ratio = 0.75;
    const int minWidth = 20;
    const int nOutIter = 5;//14;
    const int nIRLSIter = 1;
    const int nSORIter = 10;//30;
  
    char buf[256], outFile[128], outImg[128];
    std::vector<DImage> im(3), u(3), v(3), mask(3);
    DImage warp;
    UCImage ucimg;
    int cur, next, prev, width, height, i;
    OpticalFlow of;

    memset(outFile, 0, sizeof(outFile));
    memset(outImg, 0, sizeof(outImg));
    strcat(outFile, outDir);
    strcat(outFile, "%s%03d.yml");
    strcat(outImg, outDir);
    strcat(outImg, "%s%03d.jpg");

    sprintf(buf, inPattern, frameStart);
    imread(im[0], buf);
    sprintf(buf, inPattern, frameStart+1);
    imread(im[1], buf);
    
    // to save time
    DImage tmp;
    imresize(tmp, im[0], 0.5);
    tmp.copyTo(im[0]);
    sprintf(buf, outImg, "im", frameStart);
    imwrite(buf, tmp);

    imresize(tmp, im[1], 0.5);
    tmp.copyTo(im[1]);
    sprintf(buf, outImg, "im", frameStart+1);
    imwrite(buf, tmp);
    
    width = im[0].nWidth(), height = im[0].nHeight();
    mask[0].create(width, height, 1, 1);
    mask[1].create(width, height, 1, 1);
    mask[2].create(width, height, 1, 1);
    
    prev = 0, cur = 1, next = 2;
    for (i = frameStart+2; i <= frameEnd; ++i)
    {
        sprintf(buf, inPattern, i);
        imread(im[next], buf);

        // to save time
        imresize(tmp, im[next], 0.5);
        tmp.copyTo(im[next]);
        sprintf(buf, outImg, "im", i);
        imwrite(buf, tmp);

        of.stC2FFlow(u, v, im, mask, prev, as, at,
                     ratio, minWidth,
                     nOutIter, nIRLSIter, nSORIter);
        
        sprintf(buf, outFile, "u", i-2);
        imwritef(buf, u[prev]);
        sprintf(buf, outFile, "v", i-2);
        imwritef(buf, v[prev]);

        flow2color(ucimg, u[prev], v[prev]);
        sprintf(buf, outImg, "flow", i-2);
        imwrite(buf, ucimg);

        warpImage(warp, im[prev], im[cur], u[prev], v[prev]);
        sprintf(buf, outImg, "warp", i-2);
        imwrite(buf, warp);
        
        // time window move forward
        prev = cur;
        cur = next;
        next = (next + 1) % 3;
    }

    sprintf(buf, outFile, "u", i-2);
    imwritef(buf, u[prev]);
    sprintf(buf, outFile, "v", i-2);
    imwritef(buf, v[prev]);

    flow2color(ucimg, u[prev], v[prev]);
    sprintf(buf, outImg, "flow", i-2);
    imwrite(buf, ucimg);

    warpImage(warp, im[prev], im[cur], u[prev], v[prev]);
    sprintf(buf, outImg, "warp", i-2);
    imwrite(buf, warp);
    
    return 0;
}
