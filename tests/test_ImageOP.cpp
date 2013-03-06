#include "ImageProcess.h"
#include "ImageIO.h"
#include "Image.h"

#include <cstdio>

int main(int argc, char *argv[])
{
    DImage im, blur;

    im.create(3, 3, 3, 0.5);
    imprint(im);
    
    printf("imwritef()...");
    imwritef("imwritef.yml", im);
    printf("done\n");
    
    printf("imread()...");
    imreadf(im, "imwritef.yml");
    printf("done\n");

    imprint(im);
    
    // printf("imread()...");
    // imread(im, "../car1.jpg");
    // printf("done\n");

    // printf("imshow()...");
    // imshow("imshow", im);
    // imwait(0);
    // printf("done\n");

    // printf("imwrite()...");
    // imwrite("imwrite.jpg", im);
    // printf("done\n");

    // printf("gaussian blur...");
    // gSmooth(blur, im, 3, 9);
    // imshow("gaussian blur", blur);
    // imwait(0);
    // printf("done\n");

    // printf("image resize: down sample...");
    // imresize(blur, im, 0.75);
    // imshow("down sample", blur);
    // imwait(0);
    // printf("done\n");

    // printf("image resize: up sample...");
    // imresize(blur, im, 1.5);
    // imshow("up sample", blur);
    // imwait(0);
    // printf("done\n");

    // printf("image resize: down sample2...");
    // imresize(blur, im, 100, 150);
    // imshow("down sample2", blur);
    // imwait(0);
    // printf("done\n");

    // printf("image resize: up sample2...");
    // imresize(blur, im, im.nWidth()*1.5, im.nHeight()*1.4);
    // imshow("up sample2", blur);
    // imwait(0);
    // printf("done\n");

    // printf("desuarate...");
    // desuarate(blur, im);
    // imshow("desuarate", blur);
    // imwait(0);
    // printf("done\n");

    // printf("im2double...");
    // UCImage uc;
    // imread(uc, "../car2.jpg");
    // im2double(im, uc);
    // printf("uc: %d .. %d, dm: %.2f .. %.2f\n", uc.min(), uc.max(), im.min(), im.max());
    // imshow("im2double uc", uc);
    // imshow("im2double dm", im);
    // imwait(0);
    // printf("done\n");
    
    return 0;
}
