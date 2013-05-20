#ifndef _PainterlyService_H
#define _PainterlyService_H

#include <cv.h>
using cv::Mat;

#include <list>
#include <vector>
using std::list;
using std::vector;

#include "PainterlyStyle.h"
#include "SplineStroke.h"

class PainterlyService
{
public:
    PainterlyService();
	virtual ~PainterlyService(){m_init = false;};

    void setSourceImage(const Mat &src);
    void setTexture(const char *path = NULL);
    void getStrokeOrientation();
    void getEdgeMap();

    void make_spline_stroke(SplineStroke &spline_stroke, int x0, int y0, int R,
                            const Mat &dst, const Mat &ref);

    void generate_strokes(list<SplineStroke> &strokes_queue, int R, Mat canvas, const Mat &ref);
    void paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int R, const Mat &ref);
    void paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int layerId);
    void paint_layer(Mat &canvas, const list<SplineStroke> *strokes_queue, int nlayer);
    
    void render(Mat &canvas);
    void render(Mat &canvas, list<SplineStroke> *strokes_queue, int nlayer);
    void render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId);
    void render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId, Mat &ref);
    
    inline int nlayer() const{return m_nlayer;}
    
private:
	Mat m_src, m_edge_map, m_texture, m_orient, m_reference;
    int m_width, m_height;
    bool m_init;
    
    // 渲染参数
	PainterlyStyle currentStyle;
    vector<int> m_brush_radius;
};

#endif
