#include "Lighting_Effect.h"

#include <iostream>
#include <math.h>
using namespace std;

#define Pii 3.1415926
//#define PAA 0.255
#define mins(x,y) (((x)>(y))?(y):(x))

Lighting_Effect::Lighting_Effect(void){}

Lighting_Effect::~Lighting_Effect(void){}

IplImage* Lighting_Effect::Add_Light(IplImage * src,IplImage * height_map,IplImage * dst,float para){
	int w1 = src->width;
	int h1 = src->height;
	IplImage * dst_i = cvCreateImage(cvGetSize(src),src->depth,3);

	IplImage * gray,*sobel_x,*sobel_y;
	gray = cvCreateImage(cvGetSize(src),src->depth,1);
	sobel_x=cvCreateImage(cvGetSize(src),IPL_DEPTH_16S,1);
	sobel_y=cvCreateImage(cvGetSize(src),IPL_DEPTH_16S,1);
	cvCvtColor(src,gray,CV_BGR2GRAY);

	cvSobel(gray,sobel_x,1,0,3);//x方向求导
	cvSobel(gray,sobel_y,0,1,3);//y方向求导

	IplImage *sobel8_x=cvCreateImage(cvGetSize(sobel_x),IPL_DEPTH_8U,1);
	IplImage *sobel8_y=cvCreateImage(cvGetSize(sobel_y),IPL_DEPTH_8U,1);

    cvConvertScaleAbs(sobel_x,sobel8_x,1,0);//
	cvConvertScaleAbs(sobel_y,sobel8_y,1,0);

	int count =0;
	for(int i=0;i<h1;i++){
		for(int j=0;j<w1;j++){
			CvScalar result;
			CvScalar xx = cvGet2D(sobel8_x,i,j);
			CvScalar yy = cvGet2D(sobel8_y,i,j);
			double d = atan2(yy.val[0],xx.val[0]) - Pii/2.0;
			double ss = sin(d);double cc = cos(d);
			/*if(count<1000){
				//cout<<yy.val[0]<<" "<<xx.val[0]<<" ";
				cout<<d<<" ";
				count++;
			}*/
			CvScalar h_inr1;
			if(ss>cc){
				if(abs(ss -cc)>0.1){
					if(j+1<w1)
						h_inr1 = cvGet2D(height_map,i+0,j+1);
					else
						h_inr1 = cvGet2D(height_map,i+0,j+0);
				}
				else{
					if(i+1<h1&&j+1<w1)
						h_inr1 = cvGet2D(height_map,i+1,j+1);
					else if(i+1<h1&&j+1>=w1){
						h_inr1 = cvGet2D(height_map,i+1,j+0);
					}
					else if(i+1>=h1&&j+1<w1){
						h_inr1 = cvGet2D(height_map,i+0,j+1);
					}
					else{
						h_inr1 = cvGet2D(height_map,i+0,j+0);
					}
				}
			}
			else{
				if(abs(ss -cc)>0.1){
					if(i+1<h1)
						h_inr1 = cvGet2D(height_map,i+1,j+0);
					else
						h_inr1 = cvGet2D(height_map,i+0,j+0);
				}
				else{
					if(i+1<h1&&j+1<w1)
						h_inr1 = cvGet2D(height_map,i+1,j+1);
					else if(i+1<h1&&j+1>=w1){
						h_inr1 = cvGet2D(height_map,i+1,j+0);
					}
					else if(i+1>=h1&&j+1<w1){
						h_inr1 = cvGet2D(height_map,i+0,j+1);
					}
					else{
						h_inr1 = cvGet2D(height_map,i+0,j+0);
					}
				}
		    }
			
			/*if(i+1<h1){
				h_inr1 = cvGet2D(height_map,i+1,j+0);
			}
			else{
				h_inr1 = cvGet2D(height_map,i+0,j+0);
			}*/

			CvScalar h_inr2 = cvGet2D(height_map,i,j);
			double increase = h_inr1.val[0] - h_inr2.val[0];

			result = cvGet2D(dst,i,j);
			//result.val[0] = result.val[0] + PAA*increase*(mins((double)result.val[0],(255.0-result.val[0]))/127.5);
			//result.val[1] = result.val[1] + PAA*increase*(mins((double)result.val[1],(255.0-result.val[1]))/127.5);
			//result.val[2] = result.val[2] + PAA*increase*(mins((double)result.val[2],(255.0-result.val[2]))/127.5);
			result.val[0] = result.val[0] + para*increase*(mins((double)result.val[0],(255.0-result.val[0]))/127.5);
			result.val[1] = result.val[1] + para*increase*(mins((double)result.val[1],(255.0-result.val[1]))/127.5);
			result.val[2] = result.val[2] + para*increase*(mins((double)result.val[2],(255.0-result.val[2]))/127.5);
			cvSet2D(dst_i,i,j,result);
		}
	}
	cvReleaseImage(&gray);
	cvReleaseImage(&sobel_x);
	cvReleaseImage(&sobel_y);
    cvReleaseImage(&sobel8_x);
    cvReleaseImage(&sobel8_y);
	return dst_i;
}
