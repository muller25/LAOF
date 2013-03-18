#include "MotionLayers.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"
#include "Maths.h"

#include <vector>

int main(int argc, char *argv[])
{
    const char *inFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    const char *inIm = "/home/iaml/Projects/exp/lena/in/lena.avi.%03d.bmp";
    const char *outFlow = "/home/iaml/Projects/exp/lena/%s%03d.jpg";
    char buf[256];
    DImage u, v, layers, flow, centers, im, fInfo, tmp;
    UCImage flowImg;
    std::vector<DImage> vec;
    MotionLayers ml;
    int clusters, width, height;

    for (int idx = 1; idx < 2; ++idx)
    {
        sprintf(buf, inIm, idx-1);
        imread(im, buf);

        width = im.nWidth(), height = im.nHeight();
        sprintf(buf, inFile, "u", idx-1);
        assert(imreadf(u, buf));
        sprintf(buf, inFile, "v", idx-1);
        assert(imreadf(v, buf));

        vec.push_back(u);
        vec.push_back(v);
        mergec(flow, vec);
        vec.clear();

        ml.flowInfo(fInfo, flow);
        clusters = ml.cluster(centers, layers, fInfo, width, height, 2, 5, 3);
        
        printf("frame: %d, clusters: %d\n", idx-1, clusters);

        flow2color(flowImg, layers, layers);
        sprintf(buf, outFlow, "layers", idx-1);
        imwrite(buf, flowImg);

        coverLabels(tmp, im, flowImg);
        sprintf(buf, outFlow, "merge", idx-1);
        imwrite(buf, tmp);

        // refine
        printf("refine cluster...\n");
        ml.refine(centers, layers, clusters, im, flow, fInfo);

        flow2color(flowImg, layers, layers);
        sprintf(buf, outFlow, "layers-refine", idx-1);
        imwrite(buf, flowImg);

        coverLabels(tmp, im, flowImg);
        sprintf(buf, outFlow, "merge-refine", idx-1);
        imwrite(buf, tmp);
    }
    
    return 0;
}
