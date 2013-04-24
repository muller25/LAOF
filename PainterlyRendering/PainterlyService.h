#ifndef _PainterlyService_H
#define _PainterlyService_H

#include <cv.h>
using namespace cv;

#include <list>
using namespace std;

#include "PainterlyStyle.h"
#include "SplineStroke.h"

#ifndef _PI_
#define _PI_
const double PI = atan(1.0) * 4;
#endif

class PainterlyService
{
public:
    PainterlyService();
	virtual ~PainterlyService(){
        clear();
        if (m_brush_radius != NULL) delete []m_brush_radius;
        m_brush_radius = NULL;
    };
    virtual void clear();

    void setSourceImage(const Mat &src);
    void setTexture(const char *path = NULL);
    
    void getEdgeMap();
    void getStrokeOrientation();
    void make_spline_stroke(SplineStroke &spline_stroke, int x0, int y0, int R,
                            const Mat &dst, const Mat &ref);

    void generate_strokes(list<SplineStroke> &strokes_queue, int R,
                          Mat canvas, const Mat &ref);
    void paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int R, const Mat &ref);
    void paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int layerId);
    void paint_layer(Mat &canvas, const list<SplineStroke> *strokes_queue, int nlayer);
    
    void render(Mat &canvas);
    void render(Mat &canvas, list<SplineStroke> *strokes_queue, int nlayer);
    void render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId);
    void render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId, Mat &ref);
    
    inline int nlayer() const{return m_nlayer;}
    
private:
	Mat m_src, m_edge_map, m_texture, m_height_maps;
    int m_width, m_height;
    int *m_count_pass, *m_sum_pass;
    double *m_grad_orient;
    bool m_init;
    
    // 渲染参数
	PainterlyStyle currentStyle;
    int *m_brush_radius;
    int m_nlayer;
};

#endif
