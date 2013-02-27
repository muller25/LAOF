#include "Flow2Color.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

void test_makeColorWheel();
void test_computeColor();
void test_flow2color();
    
int main(int argc, char *argv[])
{
//    test_makeColorWheel();
    
//    test_computeFlowColor();

    test_flow2color();
    
    return 0;
}

void test_makeColorWheel()
{
    cout << "make color wheel... ";

    Mat wheel = makeColorWheel();

    cout << "done!" << endl;

    cout << "show color wheel... " << endl;
    int wrows = wheel.rows;
    double *pw = (double *)wheel.data;
    int wstep = wheel.step / sizeof(double);

    int factor = 50;
    int wirows = factor * (wrows / 10 + 1), wicols = factor * 10, wcha = wheel.cols;
    Mat wImg(wirows, wicols, CV_8UC(wcha));
    uchar *pwi = (uchar *)wImg.data;
    int wistep = wImg.step / sizeof(uchar), wioff;
    
    for (int wr = 0; wr < wrows; wr++)
    {
        int r = wr / 10;
        int c = wr % 10;
        for (int wir = 0; wir < factor; wir++)
        {
            for (int wic = 0; wic < factor; wic++)
            {
                wioff = (r * factor + wir) * wistep + (c * factor + wic) * wcha;
                for (int k = 0; k < wcha; k++)
                    pwi[wioff+k] = (uchar)pw[wr * wstep + k];
            }
        }
    }
    imshow("color wheel", wImg);
    waitKey(0);
    cout << "done!" << endl;
}

void test_computeColor()
{
    cout << "show flow to color image... ";

    Mat u(300, 300, CV_64F), v(300, 300, CV_64F);
    u.setTo(5), v.setTo(5);
    Mat fImg = computeColor(u, v);
    imshow("flow image", fImg);
    waitKey(0);

    cout << "done!" << endl;
}

void test_flow2color()
{
    Mat flowImg, idxImg;
    Mat u(300, 300, CV_64F), v(300, 300, CV_64F);

    cout << "test flow to color..." << endl;

    u.setTo(5), v.setTo(10);
    flow2color(u, v, flowImg, idxImg);
    imshow("flow2color flow image", flowImg);
    waitKey(0);

    imshow("flow2color idx image", idxImg);
    waitKey(0);

    // test idx image
    u.setTo(DBL_MAX), v.setTo(DBL_MAX);
    flow2color(u, v, flowImg, idxImg);
    imshow("flow2color flow image2", flowImg);
    waitKey(0);

    imshow("flow2color idx image2", idxImg);
    waitKey(0);
    
    cout << "done!" << endl;
}

