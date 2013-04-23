#include "PainterlyService.h"
#include "SplineStroke.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
#include <list>
using namespace std;

bool imreadf(Mat &m, const char *filename);
void propagate(list<SplineStroke> *strokes_queue, int nlayer,
               const Mat &src, const Mat &u, const Mat &v);
bool removeStrokes(const SplineStroke &stroke);

int main(int argc, char *argv[])
{
    if (argc != 6){
        cout << "usage: ./vr image flow out start end" << endl;
        return 1;
    }

    char im[256], flow[256], out[256], buf[256];
    int frameStart, frameEnd;

    memset(im, 0, sizeof(im));
    memset(flow, 0, sizeof(flow));
    memset(out, 0, sizeof(out));
    strcpy(im, argv[1]);
    strcpy(flow, argv[2]);
    strcpy(out, argv[3]);
    frameStart = atoi(argv[4]);
    frameEnd = atoi(argv[5]);

    Mat canvas, src, render, u, v;
    PainterlyService ps;
    const int nlayer = ps.nlayer();
    const int preRenderInterval = 6;
    int preload, frame;
    list<SplineStroke> strokes_queue[nlayer], pre_strokes_queue[nlayer];
    
    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.render(strokes_queue, nlayer);
    ps.getRenderedImage(render, src);
    sprintf(buf, out, frameStart);
    imwrite(buf, render);
    
    for (int count = 1; count <= (frameEnd - frameStart); ++count)
    {
        frame = count + frameStart;
        if (count % preRenderInterval == 1){
            preload = frameStart + count + preRenderInterval;
            if (preload > frameEnd) preload = frameEnd;
            
            sprintf(buf, im, preload);
            cout << "pre-rendering image " << buf << endl;

            src = imread(buf);
            ps.setSourceImage(src);
            ps.render(pre_strokes_queue, nlayer);
            ps.getRenderedImage(canvas, src);
        }

        sprintf(buf, im, frame);
        src = imread(buf);
        sprintf(buf, flow, "u", frame-1);
        imreadf(u, buf);
        sprintf(buf, flow, "v", frame-1);
        imreadf(v, buf);
        propagate(strokes_queue, nlayer, src, u, v);
        
        ps.setSourceImage(src);
        ps.render(strokes_queue, nlayer);
        ps.getRenderedImage(render, canvas);
        sprintf(buf, out, frame);
        imwrite(buf, render);
    }
    
	return 0;
}

bool imreadf(Mat &m, const char *filename)
{
    printf("reading file %s... ", filename);
    
    FileStorage fs(filename, cv::FileStorage::READ);

    if (!fs.isOpened())
    {
        fs.release();
        printf("failed\n");
        return false;
    }
    
    fs["matrix"] >> m;
    fs.release();
    printf("done\n");
    return true;
}

void propagate(list<SplineStroke> *strokes_queue, int nlayer,
               const Mat &src, const Mat &u, const Mat &v)
{
    assert(src.channels() == 3 && u.type() == CV_64F);

    cout << "propagating strokes..." << endl;
    
    const int threshold = 30;
    const double fadeStep = 0.2;
    
    int sstep = src.step1(), fstep = u.step1();
    uchar *ps = src.data;
    double *pu = (double*)u.data, *pv = (double*)v.data;
    
    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        for (list<SplineStroke>::iterator iter = strokes_queue[i].begin(); iter != strokes_queue[i].end(); ++iter)
        {
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);
                int foffset = p.y * fstep + p.x;
                int ioffset = p.y * sstep + p.x * 3;
                p.x += cvRound(pu[foffset]);
                p.y += cvRound(pv[foffset]);
                
                iter->setControlPoint(k, p);

                int diff = iter->ColorR() + iter->ColorG() + iter->ColorB();
                diff = abs(diff - ps[ioffset] - ps[ioffset+1] - ps[ioffset+2]);
                if (diff >= threshold) iter->fadeOut(fadeStep);
            }
        }

        cout << "remove strokes..." << endl;
        strokes_queue[i].remove_if(removeStrokes);
    }
}

bool removeStrokes(const SplineStroke &stroke)
{
    return stroke.isTransparent();
}
