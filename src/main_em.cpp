#include "Flow2Color.h"
#include "ImageIO.h"
#include "LAOF.h"

#include <vector>

char *inPattern, *outDir;
int frameStart, frameEnd;
int numOfSegs;

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
    
    char outFile[128], outImg[128];
    std::vector<DImage> im(2), mask(2), u(2), v(2), ur(2), vr(2);
    char buf[256];
    DImage warp, tmp;
    UCImage ucimg;
    int cur, next, i;

    strcat(outFile, outDir);
    strcat(outFile, "%s%d.yml");
    strcat(outImg, outDir);
    strcat(outImg, "%s%d.jpg");
    strcat(LAOF::out, outDir);
    strcat(LAOF::out, "%d-iter%d-layer%d-%s.jpg");
    
    // load images
    i = frameStart;
    sprintf(buf, inPattern, i++);
    imread(im[0], buf);

    // to save time
    // imresize(tmp, im[0], 0.5);
    // tmp.copyTo(im[0]);
    
    cur = 0, next = 1;
    for (; i <= frameEnd; ++i)
    {
        sprintf(buf, inPattern, i);
        imread(im[next], buf);

        // to save time
        // imresize(tmp, im[next], 0.5);
        // tmp.copyTo(im[next]);
        
        LAOF::EM(u, v, ur, vr, mask, im, cur, 2);

        printf("****** final dump ******\n");
        sprintf(buf, outFile, "u", i-1);
        imwritef(buf, u[cur]);
        sprintf(buf, outFile, "v", i-1);
        imwritef(buf, v[cur]);

        flow2color(ucimg, u[cur], v[cur]);
        sprintf(buf, outImg, "flow", i-1);
        imwrite(buf, ucimg);

        sprintf(buf, outFile, "ur", i);
        imwritef(buf, ur[cur]);
        sprintf(buf, outFile, "vr", i);
        imwritef(buf, vr[cur]);

        flow2color(ucimg, ur[cur], vr[cur]);
        sprintf(buf, outImg, "rflow", i);
        imwrite(buf, ucimg);

        sprintf(buf, outFile, "layers", i-1);
        imwritef(buf, mask[cur]);

        warpImage(warp, im[cur], im[next], u[cur], v[cur]);
        sprintf(buf, outImg, "warp", i-1);
        imwrite(buf, warp);

        warpImage(warp, im[next], im[cur], ur[cur], vr[cur]);
        sprintf(buf, outImg, "rwarp", i);
        imwrite(buf, warp);
        
        // time window move forward
        cur = 1 - cur;
        next = 1 - next;
    }
    
    return 0;
}
