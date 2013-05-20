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
    void make_spline_stroke(SplineStroke &spline_stroke, int x0, int y0, int R,
                            const Mat &canvas, const Mat &coverage);

    void generate_strokes(list<SplineStroke> &strokes_queue, int R, const Mat &canvas);
    void paint_layer(Mat &canvas, const list<SplineStroke> &strokes_queue, int layerId);
    void paint_layer(Mat &canvas, const list<SplineStroke> *strokes_queue, int nlayer);
    
    void render(Mat &canvas);
    void render(Mat &canvas, list<SplineStroke> *strokes_queue, int nlayer);
    void render(Mat &canvas, list<SplineStroke> &strokes_queue, int layerId);
    void fixEdges(Mat &canvas);

    inline int nLayers() const{return m_brush_radius.size();}
    
private:
	Mat m_src, m_texture, m_orient, m_reference;
    int m_width, m_height;
    bool m_init;
    
    // 渲染参数
	PainterlyStyle currentStyle;
    vector<int> m_brush_radius;
};

#endif
