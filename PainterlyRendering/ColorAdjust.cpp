#include"ColorAdjust.h"

ColorAdjust::ColorAdjust(){

}

ColorAdjust::~ColorAdjust(){
     
}

IplImage* ColorAdjust::hsv_ad( IplImage* src_image,int h,int s,int v,float nPercent){ 

IplImage* dst_image = cvCloneImage(src_image); 
cvCvtColor( src_image, dst_image, CV_BGR2HSV ); 

int hh=h - 256; 
int ss = s - 256;
int vv = v - 256;
int i; 
float max_value = 0; 

uchar lut[256][3]; 
CvMat* lut_mat; 

lut_mat = cvCreateMatHeader( 1, 256, CV_8UC3 ); 
cvSetData( lut_mat, lut,0); 

for( i = 0; i < 256; i++ ) { 
	int h1 = (i+hh); 
	int v1 = (i+vv);
	int s1 = (i+ss);
	if( h1 < 0 ) 
		h1 = 0; 
	if( h1 > 255 ) 
		h1 = 255; 
	if( v1 < 0 ) 
		v1 = 0; 
	if( v1 > 255 ) 
		v1 = 255; 
	if( s1 < 0 ) 
		s1 = 0; 
	if( s1 > 255 ) 
		s1 = 255; 

	lut[i][1] = (uchar)h1; 
	lut[i][0] = (uchar)s1; 
	lut[i][2] = (uchar)v1; 
} 

cvLUT( dst_image, dst_image, lut_mat ); 
cvCvtColor( dst_image, dst_image, CV_HSV2BGR );
	int x,y;
    float val;
    for (i = 0; i < 3; i++)//彩色图像需要处理3个通道，灰度图像这里可以删掉
    {
        for (y = 0; y < dst_image->height; y++)
        {
            for (x = 0; x < dst_image->width; x++)
            {
 
                val = ((uchar*)(dst_image->imageData + dst_image->widthStep*y))[x*3+i];
                val = 128 + (val - 128) * nPercent;
                //对灰度值的可能溢出进行处理
                if(val>255) val=255;
                if(val<0) val=0;
                ((uchar*)(dst_image->imageData + dst_image->widthStep*y))[x*3+i] = (uchar)val;
            }
        }
    }

 cvReleaseMat(&lut_mat);
 return dst_image;

} 
