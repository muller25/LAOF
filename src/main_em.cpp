#include "Flow2Color.h"
#include "ImageIO.h"
#include "LAOF.h"

#include <vector>

int main(int argc, char *argv[])
{
    const char *inImg = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *outFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    const char *outImg = "/home/iaml/Projects/exp/lena/out/%s%03d.jpg";
  
    std::vector<DImage> im(3), mask(3), u(3), v(3), ur(3), vr(3);
    char buf[256];
    DImage warp;
    UCImage ucimg;
    int cur, next, i;

    // load images
    sprintf(buf, inImg, 0);
    imread(im[0], buf);
    cur = 0, next = 1;
    for (i = 1; i < 3; ++i)
    {
        sprintf(buf, inImg, i);
        imread(im[next], buf);

        LAOF::EM(u, v, ur, vr, mask, im, cur, 2);
        
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
        cur = next;
        next = (cur + 1) % 3;
    }
    
    return 0;
}
