#include "PainterlyService.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3){
        cout << "usage: ./pr inImage outImage" << endl;
        return 1;
    }
    
    cout << "loading image " << argv[1] << endl;

    Mat src = imread(argv[1]);
    PainterlyService ps;

    ps.setSourceImage(src);
    ps.render(src);
    
    imwrite(argv[2], src);
/*
    IplImage *edgeFix = RenderingImage::operateEdge(dstImage, 100);
    cvShowImage("edge fix", edgeFix);

    IplImage *lightFix = RenderingImage::operateLight(dstImage, 25.0/100);
    cvShowImage("light fix", lightFix);
*/    
    
	return 0;
}
