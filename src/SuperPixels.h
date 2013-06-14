#ifndef _SuperPixels_H
#define _SuperPixels_H

#include <vector>
using std::vector;

#include <cv.h>
using cv::Mat;
using cv::Point;

typedef int Value;
typedef int TotalValue;
typedef int Var;

#define CV_VALUE CV_32S

class SuperPixels
{
public:
    SuperPixels(int patch_size = 30, int numIter = 2, int TYPE = 0, Value lambda = 10)
    {
        m_patch_size = patch_size;
        m_numIter = numIter;
        m_lambda = lambda;
        m_TYPE = TYPE;
        m_init = false;
    }
    
    virtual ~SuperPixels(){clear();}
    void setSourceImage(Mat &im);
    void clear(){m_init = false;}
    void PlaceSeeds();
    Value computeEnergy();
    void generate();
    void getBounds(int &seedX, int &seedY, int &startX, int &startY,
                   int &endX, int &endY, int label);
    void expandOnLabel(int label,
                       vector<int> &changeMask, vector<int> &changeMaskNew);
    void initializeLabeling();
    float computeImageVariance();
    void purturbSeeds();
    void computeWeights();
    void computeWeights(Mat &weights, int incrX, int incrY);
    void loadHEdge(const char *name);
    void loadVEdge(const char *name);
    void loadD1Edge(const char *name);
    void loadD2Edge(const char *name);
    void loadEdge(Mat &weights, const char *name);
    int countLabels();
    int saveSegmentationColor(const char *name);
    int saveSegmentationEdges(const char *name);
    int getLabels(Mat &label) {
        int nlabel = countLabels();
        m_label.copyTo(label);
        
        return nlabel;
    }
    
private:
    bool m_init;
    int m_patch_size, m_width, m_height, m_numIter, m_TYPE;
    float m_variance;
    Value m_lambda;
    Mat m_im, m_gray, m_label;
    Mat m_horizWeights, m_vertWeights, m_diag1Weights, m_diag2Weights;
    vector<Point> m_seeds;
};

#endif
