#include "MotionLayers.h"
#include "ImageIO.h"
#include "Image.h"
#include "Flow2Color.h"
#include "Maths.h"
#include "ImageProcess.h"
#include "Utils.h"
#include <vector>

/*
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
*/

int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        printf("./layers inPattern inFile outDir start end #ofLabels\n");
        return 1;
    }

    char *inPattern = argv[1];
    char *inFile = argv[2];
    char *outDir = argv[3];
    int frameStart = atoi(argv[4]);
    int frameEnd = atoi(argv[5]);
    int nlabels = atoi(argv[6]);

    char buf[256], outImg[128];
    memset(outImg, 0, sizeof(outImg));
    strcat(outImg, outDir);
    strcat(outImg, "%s%03d.jpg");

    std::vector<DImage> im(2), u(2), v(2);
    DImage centers, layers, tmp;
    UCImage ucimg;
    MotionLayers ml;
    int cur, next, width, height;
    
    sprintf(buf, inPattern, frameStart);
    imread(im[0], buf);
    sprintf(buf, inFile, "u", frameStart);
    imreadf(u[0], buf);
    sprintf(buf, inFile, "v", frameStart);
    imreadf(v[0], buf);

    width = im[0].nWidth(), height = im[0].nHeight();
    cur = 0, next = 1;
    for (int i = frameStart; i < frameEnd; ++i)
    {
        sprintf(buf, inPattern, i+1);
        imread(im[next], buf);
        sprintf(buf, inFile, "u", i+1);
        imreadf(u[next], buf);
        sprintf(buf, inFile, "v", i+1);
        imreadf(v[next], buf);

        printf("running kmeans...\n");
        ml.kcluster(centers, layers, nlabels, u[cur], v[cur]);
        printf("kmeans done\n");
        
        // show trusted points
        ucimg.create(width, height);
        for (int j = 0; j < layers.nSize(); ++j)
            if (layers[j] >= 0) ucimg[j] = 1;

        sprintf(buf, outImg, "trusted_points", i);
        imwrite(buf, ucimg);

        // show kmeans cluster result
        showLayers(ucimg, layers);
        sprintf(buf, outImg, "kcluster", i);
        imwrite(buf, ucimg);
        
        // refine layers
        printf("refining layers...\n");
        ml.refine(centers, layers, nlabels, im[cur], im[next],
                  u[cur], v[cur]);
        printf("done\n");
        
        // show refined layers
        showLayers(ucimg, layers);
        sprintf(buf, outImg, "refine", i);
        imwrite(buf, ucimg);

        // show image under each layer
        for (int l = 0; l < nlabels; ++l)
        {
            genLayerMask(ucimg, layers, l);
            cut(tmp, im[cur], ucimg);
            sprintf(buf, outImg, "seg", i*10+l);
            imwrite(buf, tmp);
        }
        
        cur = 1 - cur;
        next = 1 - next;
    }
    
    return 0;
}
