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
               const Mat &u, const Mat &v);
void removeStrokes(list<SplineStroke> *strokes_queue, int nlayer, const Mat &src);
void addStrokes(list<SplineStroke> &strokes_queue, const list<SplineStroke> &pre_strokes_queue,
                const Mat &taken);

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

    Mat canvas, src, render, u, v, taken;
    PainterlyService ps;
    const int nlayer = ps.nlayer();
    const int preRenderInterval = 4;
    int preload, frame;
    list<SplineStroke> strokes_queue[nlayer], pre_strokes_queue[nlayer];
    
    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    src.copyTo(render);
    ps.render(render, strokes_queue, nlayer);
    sprintf(buf, out, "out", frameStart);
    imwrite(buf, render);

    canvas = Mat::zeros(src.rows, src.cols, CV_8UC3);
    for (int count = 1; count <= (frameEnd - frameStart); ++count)
    {
        frame = count + frameStart;
        if (count % preRenderInterval == 1)
        {
            cout << "adding strokes..." << endl;
            for (int i = 0; i < nlayer; ++i)
            {
                if (strokes_queue[i].empty() || pre_strokes_queue[i].empty()) continue;
                
                taken = Mat::zeros(src.rows, src.cols, src.type());
                ps.paint_layer(taken, strokes_queue[i], i);
                addStrokes(strokes_queue[i], pre_strokes_queue[i], taken);
            }
            
            preload = frameStart + count + preRenderInterval;
            if (preload > frameEnd) preload = frameEnd;
            
            sprintf(buf, im, preload);
            cout << "pre-rendering image " << buf << endl;

            src = imread(buf);
            ps.setSourceImage(src);
            src.copyTo(canvas);
            ps.render(canvas, pre_strokes_queue, nlayer);

            sprintf(buf, out, "pre", preload);
            imwrite(buf, canvas);
        }

        sprintf(buf, im, frame);
        src = imread(buf);
        sprintf(buf, flow, "u", frame-1);
        imreadf(u, buf);
        sprintf(buf, flow, "v", frame-1);
        imreadf(v, buf);

        propagate(strokes_queue, nlayer, u, v);
        removeStrokes(strokes_queue, nlayer, src);

        ps.setSourceImage(src);
        canvas.copyTo(render);
        ps.paint_layer(render, strokes_queue, nlayer);
        sprintf(buf, out, "out", frame);
        imwrite(buf, render);

        render.setTo(0);
        ps.paint_layer(render, strokes_queue, nlayer);
        sprintf(buf, out, "nobg", frame);
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
               const Mat &u, const Mat &v)
{
    assert(u.type() == CV_64F);

    cout << "propagating strokes...";

    int offset, step = u.step1(), width = u.cols, height = u.rows;
    double *pu = (double*)u.data, *pv = (double*)v.data;
    
    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        for (list<SplineStroke>::iterator iter = strokes_queue[i].begin(); iter != strokes_queue[i].end(); ++iter)
        {
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);
                if (p.y >= height || p.y < 0 || p.x >= width || p.x < 0) continue;

                offset = p.y * step + p.x;
                p.x += cvRound(pu[offset]);
                p.y += cvRound(pv[offset]);
                iter->setControlPoint(k, p);
            }
        }
    }

    cout << "done" << endl;
}

void removeStrokes(list<SplineStroke> *strokes_queue, int nlayer, const Mat &src)
{
    assert(src.channels() == 3);
    cout << "removing strokes...";
    
    const int threshold = 60;
    const double fadeStep = 0.2;
    int offset, diff, color, step = src.step1(), width = src.cols, height = src.rows;
    uchar *ptr = src.data;

    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        for (list<SplineStroke>::iterator iter = strokes_queue[i].begin(); iter != strokes_queue[i].end();)
        {
            color = iter->ColorR() + iter->ColorG() + iter->ColorB();
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);

                // out of boundary
                if (p.y >= height || p.y < 0 || p.x >= width || p.x < 0) continue;
                
                offset = p.y * step + p.x * 3;
                diff = abs(ptr[offset] + ptr[offset+1] + ptr[offset+2] - color);
                if (diff > threshold) {
                    iter->fadeOut(fadeStep);
                    break;
                }
            }

            if (iter->isTransparent()) iter = strokes_queue[i].erase(iter);
            else iter++;
        }
    }

    cout << "done" << endl;
}

void addStrokes(list<SplineStroke> &strokes_queue, const list<SplineStroke> &pre_strokes_queue,
                const Mat &taken)
{
    assert(taken.channels() == 3);
    
    int offset, width = taken.cols, height = taken.rows, step = taken.step1();
    uchar *p = taken.data;

    for (list<SplineStroke>::const_iterator iter = pre_strokes_queue.begin(); iter != pre_strokes_queue.end(); ++iter)
    {
        for (int k = 0; k < iter->nPoints(); ++k)
        {
            Point point = iter->get(k);
            if (point.x >= width || point.y >= height || point.x < 0 || point.y < 0) continue;

            offset = point.y * step + point.x * 3;
            if (p[offset] != 0 || p[offset+1] != 0 || p[offset+2] != 0) continue;

            strokes_queue.push_front(*iter);
            break;
        }
    }
}
