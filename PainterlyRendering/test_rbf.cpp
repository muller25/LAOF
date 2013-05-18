#include "RBF.h"
#include <stdio.h>

#include <cv.h>
#include <highgui.h>
using namespace cv;

int main(int argc, char *argv[])
{
    if (argc < 2){
        printf("parameters: image\n");
        return -1;
    }

    Mat res, im;
    im = imread(argv[1]);
    RBF::rbf(res, im, 4);
    
    return 0;
}
