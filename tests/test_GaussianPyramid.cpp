#include "GaussianPyramid.h"
#include "Image.h"
#include "ImageIO.h"

#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    DImage im;
    imread(im, "../car1.jpg");
    const char *name = "test %d";
    char buf[256];
    
    std::cout << "constructing pyramid\n";
    
    GaussianPyramid pyr;
    pyr.ConstructPyramid(im, 0.75, 100);

    std::cout << "show pyramid\n";

    for (int i = 0; i < pyr.nLevels(); i++)
    {
        sprintf(buf, name, i);
        cout << "image size " << pyr[i].nWidth() << ", " << pyr[i].nHeight() << endl;
        imshow(buf, pyr[i]);
    }

    imwait(0);
    return 0;
}
