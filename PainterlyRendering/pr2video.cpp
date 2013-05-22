#include "PainterlyService.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 5){
        cout << "usage: ./pr2video image outDir start end" << endl;
        return 1;
    }

    char in[256], out[256], buf[256];
    strcpy(in, argv[1]);
    strcpy(out, argv[2]);
    strcat(out, "%s%03d.jpg");
    
    int frameStart = atoi(argv[3]);
    int frameEnd = atoi(argv[4]);
    Mat src, canvas;
    PainterlyService ps;    

    for (int frame = frameStart; frame <= frameEnd; ++frame)
    {
        sprintf(buf, in, frame);
        cout << "loading image " << buf << endl;

        src = imread(buf);
        ps.setSourceImage(src);
        canvas = Mat::zeros(src.rows, src.cols, src.type());
        ps.renderOverRef(canvas);

        sprintf(buf, out, "out", frame);
        imwrite(buf, canvas);

        ps.fixEdges(canvas);
        sprintf(buf, out, "out-edge", frame);
        imwrite(buf, canvas);
    }
    
	return 0;
}
