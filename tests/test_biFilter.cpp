#include <cv.h>
#include <highgui.h>
using namespace cv;

int main(int argc, char *argv[])
{
    Mat img, dimg, res;
    img = imread("/home/iaml/Projects/exp/lena/layers000.jpg");
    img.convertTo(dimg, CV_64FC3, 1./255);
    dilate(dimg, res, Mat(), Point(-1, -1), 3);
    erode(res, img, Mat(), Point(-1, -1), 1);

    imshow("res", res);
    waitKey(0);
    
    return 0;
}
