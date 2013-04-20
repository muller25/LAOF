#ifndef PAINTERLYSTYLE_H
#define PAINTERLYSTYLE_H

#pragma once
#include <string>

using std::string;

class PainterlyStyle
{
public:
	PainterlyStyle(void){}
	~PainterlyStyle(void){
		/*if (brush_radius != NULL)
			delete brush_radius;
			brush_radius = NULL;*/
	}

	virtual void initParameter(
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
		float jitter_val) ;

	virtual void setName(string name) {this->name = name;}
	virtual string toString() const {return name;}

public:
	int num_layers;            // number of brushes
	//int r[] = {8, 4, 2};            // brush sizes
	int * brush_radius;             // brush sizes
	float blur_factor;        // blur factor
	int max_stroke_length;     // maximum stroke length
	int min_stroke_length;      // minimum stroke length
	float curvature_filter;   // curvature filter
	int threshold ;            // error threshold
	float alpha ;              // opacity
	float grid_size ;          // grid size
	float jitter_r ;           // jitter red
	float jitter_g ;           // jitter green
	float jitter_b ;           // jitter blue
	float jitter_hue;         // jitter hue
	float jitter_sat;         // jitter saturation
	float jitter_val;         // jitter value

	// [9/2/2009 pluo]
	string name;
};

#endif
