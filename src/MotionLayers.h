#ifndef _MotionLayers_H
#define _MotionLayers_H

#include <cv.h>
using cv::Mat;
using cv::Point;
using cv::Vec;

#include <vector>
using std::vector;

#include "types.h"

#define DEBUG

typedef vector<PERCISION> SPPATCH;

class MotionLayers
{
public:
    void init(const Mat &im, const Mat &flow);
    void covariance(vector<PERCISION> &cov, vector< vector<Point> > &pos);
    void extract(vector< vector<Point> > &pos, int nlabel, const Mat &label);
    void extractAsHist(vector<SPPATCH> &patch, const vector< vector<Point> > &pos);
    void extractAsPixel(vector<SPPATCH> &patch, const vector< vector<Point> > &pos);
    void initSegment(int ncluster = 2);
    void refineSegment(int ncluster);
    void getLabel(Mat &lbl) const{m_label.copyTo(lbl);}
    
    static double kmdist(const PERCISION *p1, const PERCISION *p2, int start, int end);
    static double kldist(const PERCISION *p1, const PERCISION *p2, int start, int end);
    static double skldist(const PERCISION *p1, const PERCISION *p2, int start, int end);
    
private:
    Mat m_im, m_flow, m_spLbl, m_label;
    int m_width, m_height, m_nSPLbl;
    vector<SPPATCH> m_spHist, m_spPixel;
    vector< vector<Point> > m_sp2pos;
    vector<PERCISION> m_cov;
};

#endif
