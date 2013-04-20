#ifndef _GAUSSIAN_BLUR_H_
#define _GAUSSIAN_BLUR_H_
#pragma once

#include<cv.h>

class GaussianBlur
{
public:
	GaussianBlur(void){}
	~GaussianBlur(void){}

	/******************************************************************************/
	/*** gaussian blur                                                          ***/
	/*** input:  AUX_RGBImageRec * src_image                                    ***/
	/***         AUX_RGBImageRec * dst_image                                    ***/
	/*** return: nothing                                                        ***/
	/*** modify: AUX_RGBImageRec * dst_image                                    ***/
	/******************************************************************************/
	static void gaussian_blur(IplImage * dst_image, const IplImage * src_image, float sigma);
};

#endif
