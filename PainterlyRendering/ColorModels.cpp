#include "ColorModels.h"

/******************************************************************************/
/*** C Headers                                                              ***/
/******************************************************************************/
#include <stdlib.h>

/******************************************************************************/
/*** Returns the greatest number amongst the three                          ***/
/******************************************************************************/
double ColorModels::max_rgb(double r, double g, double b)
{
	if (r >= g && r >= b)
		return r;

	if (g >= r && g >= b)
		return g;

	return b;
}

/******************************************************************************/
/*** Returns the smallest number amongst the three                          ***/
/******************************************************************************/
double ColorModels::min_rgb(double r, double g, double b)
{
	if (r <= g && r <= b)
		return r;

	if (g <= r && g <= b)
		return g;

	return b;
}


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
void ColorModels::hsv_to_rgb(double *r, double *g, double * b, double h, double s, double v)
{

	if (s <= 0.0)
		s  = 0.0;

	if (s >= 1.0)
		s = 1.0;

	if (v <= 0.0)
		v  = 0.0;

	if (v >= 1.0)
		v = 1.0;

	if (s <= 0.0001)
	{
		if (h < 0.0 || h > 360.0)
		{
			*r = v;
			*g = v;
			*b = v;

			return;
		}
		else
		{
			*r = v;
			*g = v;
			*b = v;

			return;
		}
	}

	else
	{
		double f, p, q, t;
		int i;

		if (h >= 359.999)
			h = 0.0;

		h /= 60.0;
		i = (int) h;
		f = h - i;
		p = v * (1.0 - s);
		q = v * (1.0 - (s * f));
		t = v * (1.0 - (s * (1.0 - f)));

		switch(i)
		{
		case 0:
			*r = v;
			*g = t;
			*b = p;
			break;

		case 1:
			*r = q;
			*g = v;
			*b = p;
			break;

		case 2:
			*r = p;
			*g = v;
			*b = t;
			break;

		case 3:
			*r = p;
			*g = q;
			*b = v;
			break;

		case 4:
			*r = t;
			*g = p;
			*b = v;
			break;

		case 5:
			*r = v;
			*g = p;
			*b = q;
			break;
		}
	}
}


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
void ColorModels::rgb_to_hsv(double r, double g, double b, double * h, double * s, double * v)
{
	double max = max_rgb(r, g, b);
	double min = min_rgb(r, g, b);

	*v = max;
	*s = (max >= 0.0001)? ((max - min) / max) : 0.0;

	if (*s <= 0.0001)
	{
		*h = -1.0;
		return;
	}

	else
	{
		double delta = max - min;

		if (r + 0.0001 >= max)
			*h = (g - b) / delta;

		else if (g + 0.0001 >= max)
			*h = 2.0 + (b - r) / delta;


		else if (b + 0.0001 >= max)
			*h = 4.0 + (r - g) / delta;

		*h *= 60.0;

		if (*h < 0.0)
			*h += 360.0;
	}
}
