#include "ImportantMap.h"

#include <cmath>

double * ImportantMap::important_energy;
ImportantMap::ImportantMap(void){}
ImportantMap::~ImportantMap(void){}

void ImportantMap::compute_important_map(double * important_energy,IplImage * src_image){
    //VLDEnable();
	IplImage *gray,*sobel_x,*sobel_y;
	int h = src_image->height;int w = src_image->width;

	gray = cvCreateImage(cvGetSize(src_image),src_image->depth,1);
	sobel_x = cvCreateImage(cvGetSize(src_image),IPL_DEPTH_16S,1);
	sobel_y = cvCreateImage(cvGetSize(src_image),IPL_DEPTH_16S,1);
	cvCvtColor(src_image,gray,CV_BGR2GRAY);

	cvSobel(gray,sobel_x,1,0,3);//x方向求导
	cvSobel(gray,sobel_y,0,1,3);//y方向求导
		
	IplImage *sobel8_x=cvCreateImage(cvGetSize(sobel_x),IPL_DEPTH_8U,1);
	IplImage *sobel8_y=cvCreateImage(cvGetSize(sobel_y),IPL_DEPTH_8U,1);

	cvConvertScaleAbs(sobel_x,sobel8_x,1,0);
	cvConvertScaleAbs(sobel_y,sobel8_y,1,0);

	for(int i = 0;i<h;i++){
		for(int j =0;j<w;j++){
		    CvScalar xx = cvGet2D(sobel8_x,i,j);
			CvScalar yy = cvGet2D(sobel8_y,i,j);
			int gx = (uchar)xx.val[0];
			int gy = (uchar)yy.val[0];
			important_energy[i*w+j] = sqrt(pow((double)gx,2)+pow((double)gy,2));
		}
	}

	/*释放内存*/
	cvReleaseImage(&gray);
	cvReleaseImage(&sobel_x);
	cvReleaseImage(&sobel_y);
	cvReleaseImage(&sobel8_x);
	cvReleaseImage(&sobel8_y);
	//VLDDisable();
}
