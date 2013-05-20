#include "PainterlyService.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3){
        cout << "usage: ./pr inImage outDir" << endl;
        return 1;
    }
    
    cout << "loading image " << argv[1] << endl;

    Mat src = imread(argv[1]);
    PainterlyService ps;

    ps.setSourceImage(src);
    Mat canvas = Mat::zeros(src.rows, src.cols, src.type());
    ps.render(canvas);

    char buf[256];
    sprintf(buf, "%srender.jpg", argv[2]);
    imwrite(buf, canvas);

    ps.fixEdges(canvas);
    sprintf(buf, "%srender-fix-edge.jpg", argv[2]);
    imwrite(buf, canvas);
    
/*
    IplImage *edgeFix = RenderingImage::operateEdge(dstImage, 100);
    cvShowImage("edge fix", edgeFix);

    IplImage *lightFix = RenderingImage::operateLight(dstImage, 25.0/100);
    cvShowImage("light fix", lightFix);
*/    
    
	return 0;
}
