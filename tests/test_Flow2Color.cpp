#include "Flow2Color.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    cout << "make color wheel... ";

    Mat wheel = makeColorWheel();

    cout << "done!" << endl;

    int wrows = wheel.rows;
    double *pw = (double *)wheel.data;
    int wstep = wheel.step / sizeof(double);

    int wirows = 50, wicols = wrows * wirows, wcha = wheel.cols;
    Mat wImg(wirows, wicols, CV_8UC(wcha));
    uchar *pwi = (uchar *)wImg.data;
    int wistep = wImg.step / sizeof(uchar), wioff;
    
    for (int wr = 0; wr < wrows; wr++)
    {
        for (int wir = 0; wir < wirows; wir++)
        {
            for (int wic = 0; wic < wirows; wic++)
            {
                wioff = wir * wistep + ((wr*wirows) + wic) * wcha;
                for (int k = 0; k < wcha; k++)
                    pwi[wioff+k] = (uchar)pw[wr * wstep + k];
            }
        }
    }
    
    cout << "show color wheel" << endl;
    imshow("color wheel", wImg);
    waitKey(0);
    
    // cout << "show flow to color image" << endl;
    
    return 0;
}
