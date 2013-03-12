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
    DImage u, v, layers, flow, rflow, centers, im, gf, gray;
    UCImage flowImg;
    std::vector<DImage> vec;
    MotionLayers ml;
    int clusters;

    for (int idx = 1; idx < 2; ++idx)
    {
        sprintf(buf, inIm, idx-1);
        imread(im, buf);

        sprintf(buf, inFile, "u", idx-1);
        imreadf(u, buf);
        sprintf(buf, inFile, "v", idx-1);
        imreadf(v, buf);

        vec.push_back(u);
        vec.push_back(v);
        mergec(flow, vec);
        vec.clear();

        sprintf(buf, inFile, "ru", idx);
        imreadf(u, buf);
        sprintf(buf, inFile, "rv", idx);
        imreadf(v, buf);

        vec.push_back(u);
        vec.push_back(v);
        mergec(rflow, vec);
        vec.clear();
        
        clusters = ml.cluster(centers, layers, im, flow, rflow);
        printf("frame: %d, clusters: %d\n", idx-1, clusters);

        flow2color(flowImg, layers, layers);
        sprintf(buf, outFlow, "layers", idx-1);
        imwrite(buf, flowImg);

        // filter
        desuarate(gray, im);
        GuidedFilter(gf, layers, gray, 8, 1e-6);
    
        flow2color(flowImg, gf, gf);
        sprintf(buf, outFlow, "layersgf", idx-1);
        imwrite(buf, flowImg);
    }
    
    return 0;
}
