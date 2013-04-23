#include "PainterlyService.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2){
        cout << "usage: ./pr imagePath" << endl;
        return 1;
    }
    
    cout << "loading image " << argv[1] << endl;

    Mat src = imread(argv[1]);
    PainterlyService ps;

    ps.setSourceImage(src);
    ps.render();
    
    imwrite("rendered.jpg", ps.getRenderedImage());
/*
    IplImage *edgeFix = RenderingImage::operateEdge(dstImage, 100);
    cvShowImage("edge fix", edgeFix);

    IplImage *lightFix = RenderingImage::operateLight(dstImage, 25.0/100);
    cvShowImage("light fix", lightFix);
*/    
    
	return 0;
}
