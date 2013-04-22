#ifndef PAINTERLYSERVICE_H
#define PAINTERLYSERVICE_H

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

/******************************************************************************/
/*** OpenGL & GLUT Headers                                                  ***/
/******************************************************************************/
//#include <GL/gl.h>
/******************************************************************************/
/*** Multi-Threading Headers                                                ***/
/******************************************************************************/
//#include <pthread.h>

/******************************************************************************/
/*** Painterly Rendering Headers                                            ***/
/******************************************************************************/

#include "GaussianBlur.h"
#include "SplineStroke.h"
#include "CubicBSpline.h"
#include "ColorModels.h"
#include "PainterlyStyle.h"

#include <cv.h>
#include <highgui.h>
#include <map>
#include <vector>
using std::map;
using std::vector;

class PainterlyService
{
public:
	PainterlyService(void){}
	~PainterlyService(void){}
	/******************************************************************************/
	/*** Images                                                                 ***/
	/******************************************************************************/
	static IplImage * src_image;            // source image
	static IplImage * dst_image;            // destination image
	static IplImage * ref_image;            // reference image
	//static IplImage * buf_image;            // buffer image
	static int * dif_image;                        // difference image
	static IplImage * tex_image; //笔刷纹理
	static IplImage * msk_image; //遮罩图
	static IplImage * map_image;//高度图
	static IplImage * brush_tex;//纹理笔刷
	static IplImage * brush_test;//测试笔刷是怎样的
	static IplImage * edge_image; //边缘图
	static IplImage * height_maps;//新高度图
	static IplImage * new_msk_maps;//新透明图
	static double * grad_orient;
	static IplImage * gradx_maps; //x方向梯度
	static IplImage * grady_maps; //y方向梯度
	static int * count_pass;//通过像素点数次
	static int * sum_pass;//通过像素点像素值之和
	static bool useTexture;      

	static IplImage* color_map;//颜色值
	// static unsigned char * depth_image;            // depth image
	static float * depth_image;                    // depth image

	/******************************************************************************/
	/*** Drawing States                                                         ***/
	/******************************************************************************/
	static int NUM_STROKES_PER_FRAME ;     // number of strokes per frame
	static int ANIMATE_STROKE ;             // animation flag
	static float ANIMATE_SPEED;          // animation speed
	static int layer_done ;                 // current layer done
	static int paint_done;                 // painting done
	static float animate_time ;           // current animation frame
	static int is_locked_for_animation ;    // queue locked status
	static int first_run;                  // first time entering draw_scene
	static time_t current_time;            // random number generator seed
	static int output_written;             // output file written       

	/******************************************************************************/
	/*** Strokes Queue                                                          ***/
	/******************************************************************************/
	// static Queue * strokes_queue;


	/******************************************************************************/
	/*** Posix Threads and Mutex Locks                                          ***/
	/******************************************************************************/
	//static pthread_t paint_thread;                 // handle for paint thread
	//static pthread_mutex_t mutex_lock_queue;	    // mutex lock for queue
	//static pthread_mutex_t mutex_lock_next_layer;  // mutex lock for next layer

	/******************************************************************************/
	/*** Descend Function For QSORT                                             ***/
	/******************************************************************************/
	static int descend(const void * a, const void * b);

	static SplineStroke *make_spline_stroke(int x0, int y0, int R, const IplImage *ref_image);

	static void difference_image(int *dif_image, const IplImage *ref_image, 
		const IplImage *dst_image, int R);

	static void generate_strokes(IplImage *dst_image, IplImage *ref_image, int R,
                                 vector<SplineStroke> &strokes_queue);
	static void paint_layer(IplImage *dst_image, IplImage *ref_image, int R,
                            const vector<SplineStroke> &strokes_queue);

	static void loadStyles(string, PainterlyStyle);
	static CvPoint bernstein(float u, CvPoint *p);
    static CvPoint pointAdd(CvPoint p, CvPoint q);
    static CvPoint pointTimes(float c, CvPoint p);

	static map<string, PainterlyStyle> styles;
	static PainterlyStyle currentStyle;
};

#endif
