#include "PainterlyStyle.h"

void PainterlyStyle::initParameter(
				   int num_layers,
				   int * brush_radius,
				   float blur_factor,
				   int max_stroke_length,
				   int min_stroke_length,
				   float curvature_filter,
				   int threshold,
				   float alpha,
				   float grid_size,
				   float jitter_r,
				   float jitter_g,
				   float jitter_b,
				   float jitter_hue,
				   float jitter_sat,
				   float jitter_val) 
{
	this->num_layers = num_layers;
	this->brush_radius = brush_radius;
	this->blur_factor = blur_factor;
	this->max_stroke_length = max_stroke_length;
	this->min_stroke_length = min_stroke_length;
	this->curvature_filter = curvature_filter;
	this->threshold = threshold;
	this->alpha = alpha;
	this->grid_size = grid_size;
	this->jitter_r = jitter_r;
	this->jitter_g = jitter_g;
	this->jitter_b = jitter_b;
	this->jitter_hue = jitter_hue;
	this->jitter_sat = jitter_sat;
	this->jitter_val = jitter_val;
}
