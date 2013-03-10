#include "MotionLayers.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"

#include <vector>

int main(int argc, char *argv[])
{
    const char *inFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    const char *inIm = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *outFlow = "/home/iaml/Projects/exp/lena/%s%03d.jpg";
    char buf[256];
    DImage u, v, layers, flow, rflow, centers, im;
    UCImage flowImg;
    int idx = 2;
    std::vector<DImage> vec;

    sprintf(buf, inIm, idx);
    imread(im, buf);
    sprintf(buf, inFile, "u", idx);
    imreadf(u, buf);
    sprintf(buf, inFile, "v", idx);
    imreadf(v, buf);
    
    vec.push_back(u);
    vec.push_back(v);
    mergec(flow, vec);
    vec.clear();
    
    sprintf(buf, inFile, "ru", idx-1);
    imreadf(u, buf);
    sprintf(buf, inFile, "rv", idx-1);
    imreadf(v, buf);

    vec.push_back(u);
    vec.push_back(v);
    mergec(rflow, vec);
    vec.clear();
    
    MotionLayers ml;
    int clusters;
    clusters = ml.cluster(centers, layers, im, flow, rflow);

    printf("clusters: %d\n", clusters);
    flow2color(flowImg, layers, layers);

    sprintf(buf, outFlow, "layers", idx);
    imwrite(buf, flowImg);
    
    return 0;
}
