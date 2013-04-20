#ifndef COLORMODELS_H
#define COLORMODELS_H
#pragma once

class ColorModels
{
public:
	ColorModels(){}
	~ColorModels(){}

	/******************************************************************************/
	/*** Convert Color Model HSV to RGB                                         ***/
	/***                                                                        ***/
	/*** Input: h: hue                                                          ***/
	/***        s: saturation                                                   ***/
	/***        v: value                                                        ***/
	/*** Modifies: *r: red                                                      ***/
	/***           *g: green                                                    ***/
	/***           *b: blue                                                     ***/
	/******************************************************************************/
	static void hsv_to_rgb(double * r, double * g, double * b, double h, double s, double v);


	/******************************************************************************/
	/*** Convert Color Model RGB to HSV                                         ***/
	/***                                                                        ***/
	/*** Input: r: red                                                          ***/
	/***        g: green                                                        ***/
	/***        b: blue                                                         ***/
	/*** Modifies: *h: hue                                                      ***/
	/***           *s: saturation                                               ***/
	/***           *v: value                                                    ***/
	/******************************************************************************/
	static void rgb_to_hsv(double r, double g, double b, double * h, double * s, double * v);

	/******************************************************************************/
	/*** Returns the smallest number amongst the three                          ***/
	/******************************************************************************/
	static double min_rgb(double r, double g, double b);

	/******************************************************************************/
	/*** Returns the greatest number amongst the three                          ***/
	/******************************************************************************/
	static double max_rgb(double r, double g, double b);

};

#endif
