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
void addStrokes(list<SplineStroke> *strokes_queue, const list<SplineStroke> *pre_strokes_queue,
                int nlayer, const Mat &taken);

int main(int argc, char *argv[])
{
    if (argc != 6){
        cout << "usage: ./vr image flow outDir start end" << endl;
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
    strcat(out, "%s%03d.jpg");
    frameStart = atoi(argv[4]);
    frameEnd = atoi(argv[5]);

    Mat canvas, src, render, u, v;
    PainterlyService ps;
    const int nlayer = ps.nLayers();
    const int preRenderInterval = 10;
    int preload, frame;
    list<SplineStroke> strokes_queue[nlayer], pre_strokes_queue[nlayer];
    
    preload = frameStart + preRenderInterval;
    if (preload > frameEnd) preload = frameEnd;
            
    sprintf(buf, im, preload);
    cout << "pre-rendering image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.render(canvas, pre_strokes_queue, nlayer);
    src.copyTo(canvas);
    ps.paint_layer(canvas, strokes_queue, nlayer);
    ps.fixEdges(canvas);
    sprintf(buf, out, "pre", preload);
    imwrite(buf, canvas);

    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.render(render, strokes_queue, nlayer);
    canvas.copyTo(render);
    ps.paint_layer(render, strokes_queue, nlayer);
    ps.fixEdges(render);
    
    // render = Mat::zeros(src.rows, src.cols, src.type());
    // src.setTo(0);
    // ps.paint_layer(render, strokes_queue[ll], ll, src);
    
    sprintf(buf, out, "out", frameStart);
    imwrite(buf, render);

    // canvas = Mat::zeros(src.rows, src.cols, CV_8UC3);
    for (int count = 1; count <= (frameEnd - frameStart); ++count)
    {
        frame = count + frameStart;
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
        ps.fixEdges(render);
        sprintf(buf, out, "out", frame);
        imwrite(buf, render);

        render.setTo(0);
        ps.paint_layer(render, strokes_queue, nlayer);
        ps.fixEdges(render);
//        ps.paint_layer(render, strokes_queue[ll], ll, canvas);
        sprintf(buf, out, "nobg", frame);
        imwrite(buf, render);

        if (count % preRenderInterval == 0)
        {
            addStrokes(strokes_queue, pre_strokes_queue, nlayer);
            
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
    assert(src.type() == CV_8UC3);
    cout << "removing strokes...";
    
    const int threshold = 20 * 3;
    const double fadeStep = -0.1;
    int width = src.cols, height = src.rows;

    for (int i = 0; i < nlayer; ++i)
    {
        if (strokes_queue[i].empty()) continue;

        list<SplineStroke>::iterator iter = strokes_queue[i].begin();
        int R = iter->getRadius();
        Vec3i bcolor = iter->getColor();
        while (iter != strokes_queue[i].end())
        {
            int pointDiff = 0;
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                // color diff of brush covered areas
                Point p = iter->get(k);
                int count = 0, diff = 0;
                for (int hh = -R; hh <= R; ++hh)
                    for (int ww = -R; ww <= R; ++ww)
                    {
                        int hhh = p.y + hh;
                        int www = p.x + ww;
                        if (hhh >= height || www >= width || hhh < 0 || www < 0) continue;

                        Vec3i scolor = src.at<Vec3b>(hhh, www);
                        Vec3i diffcolor = scolor - bcolor;
                        diff += norm(diffcolor, NORM_L1);
                        count++;
                    }

                if (count <= 0) count = 1;
                diff /= count;
                if (diff >= threshold) ++pointDiff;
            }

            // more than half of control points are diff to underneath color
            if (pointDiff >= 0.5 * iter->nPoints()) iter->changeOpacity(fadeStep);
            if (iter->isTransparent()) iter = strokes_queue[i].erase(iter);
            else iter++;
        }
    }

    cout << "done" << endl;
}

void addStrokes(list<SplineStroke> *strokes_queue, const list<SplineStroke> *pre_strokes_queue,
                int nlayer, Mat &src)
{
    assert(taken.channels() == 3);
    cout << "adding strokes..." << endl;

    int width = src.cols, height = src.rows;
    Mat coverage(height, width, CV_32S);
    for (int i = 0; i < nlayer; ++i)
    {
        if (pre_strokes_queue[i].empty()) continue;

        coverage.setTo(0);
        getCoverage(coverage, strokes_queue[i]);
        for (list<SplineStroke>::const_iterator iter = pre_strokes_queue[i].begin(); iter != pre_strokes_queue[i].end(); ++iter)
        {
            int R = iter->getRadius(), covered = 0;
            int area = (2 * R + 1) * (2 * R + iter->nPoints());
            Point prePoint = iter->get(0);
            for (int h = -R; h <= R; ++h)
                for (int w = -R; w <= R; ++w)
                {
                    
                }
            
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point point = iter->get(k);
                if (point.x >= width || point.y >= height || point.x < 0 || point.y < 0) continue;

                

                strokes_queue[i].push_front(*iter);
                break;
            }
        }
    }
}

void getCoverage(Mat &coverage, const list<SplineStroke> &strokes_queue)
{
    assert(coverage.data != NULL && coverage.type() == CV_32S);

    coverage.setTo(0);
    int width = coverage.cols, height = coverage.rows;
    for (list<SplineStroke>::const_iterator iter = strokes_queue.begin(); iter != strokes_queue.end(); iter++)
    {
        int R = iter->getRadius();
        for (int i = 0; i < iter->nPoints(); ++i)
        {
            Point point = iter->get(i);
            for(int hh = -R; hh <= R; hh++)
                for(int ww = -R; ww <= R; ww++){
                    int bhh = point.y + hh, bww = point.x + ww;
                    if (bhh < 0 || bhh >= height || bww >= width || bww < 0) continue;

                    coverage.at<int>(bhh, bww) += 1;
                }
        }
    }
}
