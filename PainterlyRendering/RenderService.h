#ifndef RENDERSERVICE_H
#define RENDERSERVICE_H

#pragma once
/******************************************************************************/
/*** C Headers                                                              ***/
/******************************************************************************/
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "gaussianblur.h"
#include "queue.h"
#include "splinestroke.h"
#include "cubicbspline.h"
#include "color_models.h"
#include "PainterlyStyle.h"
#include "cv.h"
#include "highgui.h"

#include <map>



using std::map;

typedef struct Brush_color_2{
	unsigned int r,g,b;
}Brush_color_2;

class RenderService
{
public:

	RenderService(void){}
	~RenderService(void){}

	static IplImage * src_image;            // source image
	static IplImage * dst_image;            // destination image
	static IplImage * ref_image;            // reference image
	static int * dif_image;                 // difference image
    static float animate_time ;           // current animation frame

	static IplImage * brush;

	//static IplImage * brush_tex;//纹理笔刷
	static IplImage * brush_tex1;//纹理笔刷
	static IplImage * brush_tex2;//纹理笔刷
	static IplImage * brush_tex3;//纹理笔刷
	static IplImage * brush_tex4;//纹理笔刷
	static IplImage * brush_tex5;//纹理笔刷
    static IplImage * brush_tex6;//纹理笔刷

	static IplImage * edge_image; //边缘图
	static IplImage * height_maps;//新高度图
	static IplImage * new_msk_maps;//新透明图

	static int * count_pass;//通过像素点数次
	static int * sum_pass;//通过像素点像素值之和
	static float * depth_image;// depth image 

	static Queue * strokes_queue;//笔划队列


	/*渲染函数*/
	static int descend(const void * a, const void * b);
	static SplineStroke* RenderService::make_spline_stroke(int x0, int y0, int R, const IplImage * ref_image);
	static void paint_layer(IplImage * dst_image,IplImage * ref_image,int R);
	static void RenderService::difference_image(int * dif_image, const IplImage * ref_image, const IplImage * dst_image, int R);
	static void RenderService::loadStyles(string name, PainterlyStyle style);
    static void getStrokeSize(float aver_gradient,int R,int *Rs);

	static map<string, PainterlyStyle> styles;
	static PainterlyStyle currentStyle;


};

#endif