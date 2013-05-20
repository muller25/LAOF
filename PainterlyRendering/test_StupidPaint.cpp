#include "StupidPaint.h"

#include <iostream>
using namespace std;

#include <cv.h>
#include <highgui.h>
using namespace cv;

int main(int argc, char *argv[])
{
    if (argc < 3) {
        cout << "parameters: input_image output_image" << endl;
        return -1;
    }

    Mat src, dst;

    src = imread(argv[1]);
    StupidPaint::loadTexture();
    StupidPaint::paint(dst, src);

    imshow("render", dst);
    waitKey(0);    

    imwrite(argv[2], dst);

    return 0;
}
