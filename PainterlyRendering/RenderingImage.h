#ifndef RENDERINGIMAGE_H
#define RENDERINGIMAGE_H

#include <cv.h>

//#include "RBF.h"
class RenderingImage {

public:
	RenderingImage();
	~RenderingImage();
    static IplImage* Processing(IplImage * src,IplImage* edgeImage);
    static void getStrokeOrientation(IplImage* src_image,double * orientation);
    static IplImage* operateLight(IplImage* res,float para);
    static IplImage* opereateEdge(IplImage* res,int para);

	//static int max_len;
	static int max_stroke_length;
	static int min_stroke_length;
	static int threshold;
	static float grid_size;
    static const int nlayer = 6;
	static int R[nlayer];
    static char brushMapPath[256];
};

#endif
