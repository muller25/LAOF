#ifndef _PainterlyService_H
#define _PainterlyService_H

#include <cv.h>
using namespace cv;

#include <vector>
using std::vector;

#include "PainterlyStyle.h"
#include "SplineStroke.h"

class PainterlyService
{
public:
	PainterlyService();
	virtual ~PainterlyService(){clear();};
    virtual void clear();

    void setSourceImage(const Mat &src);
    void setTexture(const char *path = NULL);
    void getEdgeMap();
    void getStrokeOrientation();

	SplineStroke *make_spline_stroke(int x0, int y0, int R, const Mat &ref);
    void difference_image(int *diff, const Mat &ref, int R);
    void generate_strokes(vector<SplineStroke> &strokes_queue, const Mat &ref, int R);
    void paint_layer(const Mat &ref, int R, const vector<SplineStroke> &strokes_queue);
    void render();

    Mat &getRenderedImage() const{return m_dst;}
    Mat &getRenderedImage() {return m_dst;}
    
private:
	Mat m_src, m_dst, m_edge_map, m_texture;
    int m_width, m_height;
    int *m_count_pass, *m_sum_pass;
    double *m_grad_orient;
    bool m_init;
    
    // 渲染参数
	PainterlyStyle currentStyle;
    const int m_nlayer = 6;
    const int m_brush_radius[m_nlayer];
};

#endif
