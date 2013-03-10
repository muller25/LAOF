#include "MotionLayers.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"

#include <vector>

int main(int argc, char *argv[])
{
    const char *inFile = "/home/iaml/Projects/exp/lena/out/%s%03d.yml";
    const char *outFlow = "/home/iaml/Projects/exp/lena/%s%03d.jpg";
    char buf[256];
    DImage u, v, layers, flow, rflow, centers, zeros;
    UCImage flowImg;
    int idx = 1;
    
    std::vector<DImage> vec, rvec;
    
    sprintf(buf, inFile, "u", idx);
    imreadf(u, buf);
    sprintf(buf, inFile, "v", idx);
    imreadf(v, buf);

    // flow2color(flowImg, u, v);
    // imshow("flow", flowImg);
    // imwait(0);
    
    vec.push_back(u);
    vec.push_back(v);
    printf("merging channels...");
    mergec(flow, vec);
    printf("done\n");
    
    sprintf(buf, inFile, "ru", idx-1);
    imreadf(u, buf);
    sprintf(buf, inFile, "rv", idx-1);
    imreadf(v, buf);

    // flow2color(flowImg, u, v);
    // imshow("rflow", flowImg);
    // imwait(0);

    rvec.push_back(u);
    rvec.push_back(v);

    printf("merging channels...");
    mergec(rflow, rvec);
    printf("done\n");

    MotionLayers ml;
    int clusters;
    clusters = ml.flowCluster(centers, layers, flow, rflow);

    zeros.create(layers.nWidth(), layers.nHeight());
    printf("clusters: %d\n", clusters);
    flow2color(flowImg, layers, zeros);

    sprintf(buf, outFlow, "layers", 0);
    imwrite(buf, flowImg);

    return 0;
}
