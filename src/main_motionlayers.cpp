#include "MotionLayers.h"

#include <cstdio>
#include <vector>
using std::vector;

#include <cv.h>
#include <highgui.h>
using namespace cv;

bool imreadf(Mat &m, const char *filename)
{
    cv::FileStorage fs(filename, cv::FileStorage::READ);

    if (!fs.isOpened())
    {
        fs.release();
        return false;
    }

    fs["matrix"] >> m;
    fs.release();
    
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        printf("./layers imPattern inFile outDir start end #ofLabels\n");
        return 1;
    }

    char *imPattern = argv[1];
    char *inFile = argv[2];
    char *outDir = argv[3];
    int frameStart = atoi(argv[4]);
    int frameEnd = atoi(argv[5]);
    int nlabel = atoi(argv[6]);

    char buf[256], outImg[128];
    memset(outImg, 0, sizeof(outImg));
    strcat(outImg, outDir);
    strcat(outImg, "%s%03d.jpg");

    MotionLayers ml;
    Mat im, flow, u, v, lbl;
    for (int i = frameStart; i < frameEnd; ++i)
    {
        sprintf(buf, imPattern, i);
        im = imread(buf);
        sprintf(buf, inFile, "u", i);
        imreadf(u, buf);
        sprintf(buf, inFile, "v", i);
        imreadf(v, buf);
        Mat arr[] = {u, v};
        merge(arr, 2, flow);

        printf("init...\n");
        ml.init(im, flow);
        printf("done\n");

        printf("init segment...\n");
        ml.initSegment(nlabel);
        printf("done\n");

        printf("refine segment...\n");
        ml.refineSegment(nlabel);
        printf("done\n");
    }
    
    return 0;
}
