#include "RenderingImage.h"
#include <cv.h>
#include <highgui.h>
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 5){
        cout << "usage: ./vr image flow out start end" << endl;
        return 1;
    }

    char im[256]. flow[256], out[256], buf[256];
    int frameStart, frameEnd;
    
    memset(im, 0, sizeof(im));
    memset(flow, 0, sizeof(flow));
    memset(out, 0, sizeof(out));
    strcpy(im, argv[1]);
    strcpy(flow, argv[2]);
    strcpy(out, argv[3]);
    frameStart = atoi(argv[4]);
    frameEnd = atoi(argv[5]);

    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;
    IplImage *srcImage = cvLoadImage(buf);
    CvSize srcSize = cvGetSize(srcImage);

    for (int frame = frameStart; frame <= frameEnd; ++frame){
        
        // 将彩色图像转换为灰度图像
    IplImage *grayImage = cvCreateImage(srcSize, IPL_DEPTH_8U, 1);
    IplImage *smoothImage = cvCreateImage(srcSize, IPL_DEPTH_8U, 1);
    cvCvtColor(srcImage, grayImage, CV_BGR2GRAY);    
    cvSmooth(grayImage, smoothImage, CV_BLUR, 3, 3, 0);
    cvNot(smoothImage, grayImage);
    
    // 对灰度图像进行边缘检测
    float edgeThreshold = 100;
    IplImage *edgeImage = cvCreateImage(srcSize, IPL_DEPTH_8U, 3);
    cvCanny(grayImage, smoothImage, edgeThreshold, edgeThreshold*3, 3);
    cvZero(edgeImage);
    cvCopy(srcImage, edgeImage, smoothImage);

    cvReleaseImage(&grayImage);
    cvReleaseImage(&smoothImage);

    // 设置渲染参数
    RenderingImage::max_stroke_length = 16;
    RenderingImage::min_stroke_length = 2;
    RenderingImage::threshold = 46;
    RenderingImage::grid_size = 2.0;
    memset(RenderingImage::brushMapPath, 0, sizeof(RenderingImage::brushMapPath));
    strcpy(RenderingImage::brushMapPath, argv[2]);
    
    memset(RenderingImage::R, 0, sizeof(int) * RenderingImage::nlayer);
    // 不绘制 R=32
    int radius = 1;
    for (int i = RenderingImage::nlayer-1; i >= 1; i--){
        RenderingImage::R[i] = radius;
        radius <<= 1;
    }
    
    IplImage* dstImage = RenderingImage::Processing(srcImage, edgeImage);
    cvShowImage("rendered", dstImage);
/*
    IplImage *edgeFix = RenderingImage::operateEdge(dstImage, 100);
    cvShowImage("edge fix", edgeFix);

    IplImage *lightFix = RenderingImage::operateLight(dstImage, 25.0/100);
    cvShowImage("light fix", lightFix);
*/    
    cvWaitKey();

    cvReleaseImage(&srcImage);
    cvReleaseImage(&edgeImage);
    cvReleaseImage(&dstImage);
    }
    
	return 0;
}
