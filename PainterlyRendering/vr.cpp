#include "PainterlyService.h"
#include "SplineStroke.h"

#include <cv.h>
#include <highgui.h>
using namespace cv;

#include <iostream>
#include <list>
using namespace std;

#define DEBUG

bool imreadf(Mat &m, const char *filename);
void propagate(list<SplineStroke> *next_queue, const list<SplineStroke> *cur_queue,
               size_t nlayer, const Mat &u, const Mat &v);
void removeStrokes(list<SplineStroke> *strokes_queue, size_t nlayer, const Mat &src);
void removeStrokes2(list<SplineStroke> *strokes_queue, size_t nlayer, const Mat &src);
void addStrokes(list<SplineStroke> *strokes_queue, const list<SplineStroke> *pre_strokes_queue,
                const size_t nlayer, Mat &src);
void updateStrokes(list<SplineStroke> *queue, size_t nlayer, const Mat &src);

char im[256], flow[256], out[256];
int frameStart, frameEnd, frame;

int main(int argc, char *argv[])
{
    if (argc != 6){
        cout << "usage: ./vr image flowDir outDir start end" << endl;
        return 1;
    }

    memset(im, 0, sizeof(im));
    memset(flow, 0, sizeof(flow));
    memset(out, 0, sizeof(out));
    strcpy(im, argv[1]);
    strcpy(flow, argv[2]);
    strcat(flow, "%s%03d.yml");
    strcpy(out, argv[3]);
    strcat(out, "%s%03d.jpg");
    frameStart = atoi(argv[4]);
    frameEnd = atoi(argv[5]);

    char buf[256];
    Mat canvas, src, u, v, smooth;
    PainterlyService ps;
    const size_t nlayer = ps.nLayers();
    const int interval = 10;
    list<SplineStroke> strokes_queue[interval][nlayer], new_queue[nlayer], add_queue[nlayer];

    sprintf(buf, im, frameStart);
    cout << "loading image " << buf << endl;

    src = imread(buf);
    ps.setSourceImage(src);
    ps.generate_strokes(strokes_queue[0], nlayer);
    canvas = Mat::zeros(src.rows, src.cols, CV_8UC3);

    int pid, qstart = 0, qend = 1, dumpedFrame = -1;
    for (int count = 1; count <= (frameEnd - frameStart); ++count, qend = (qend+1) % interval)
    {
        frame = count + frameStart;
        sprintf(buf, im, frame);
        cout << "loading image " << buf << endl;

        src = imread(buf);
        sprintf(buf, flow, "u", frame-1);
        imreadf(u, buf);
        sprintf(buf, flow, "v", frame-1);
        imreadf(v, buf);
        
        // stroke propagation
        pid = (qend-1+interval) % interval;
//        cout << "******** qstart = " << qstart << ", qend = " << qend << ", pid = " << pid << " ********" << endl;
        propagate(new_queue, strokes_queue[pid], nlayer, u, v);
#ifdef DEBUG
        canvas.setTo(0);
        ps.paint_layer(canvas, new_queue, nlayer);
        sprintf(buf, out, "propagate", frame);
        imwrite(buf, canvas);
#endif

        // update strokes
        updateStrokes(new_queue, nlayer, src);
#ifdef DEBUG
        canvas.setTo(0);
        ps.paint_layer(canvas, new_queue, nlayer);
        sprintf(buf, out, "update", frame);
        imwrite(buf, canvas);
#endif
        
        // remove strokes
        removeStrokes(new_queue, nlayer, src);
#ifdef DEBUG
        canvas.setTo(0);
        ps.paint_layer(canvas, new_queue, nlayer);
        sprintf(buf, out, "remove", frame);
        imwrite(buf, canvas);
#endif

        // add strokes
        addStrokes(add_queue, new_queue, nlayer, src);
#ifdef DEBUG
        canvas.setTo(0);
        ps.paint_layer(canvas, add_queue, nlayer);
        sprintf(buf, out, "add", frame);
        imwrite(buf, canvas);
#endif
        for (size_t l = 0; l < nlayer; ++l)
            new_queue[l].insert(new_queue[l].begin(), add_queue[l].begin(), add_queue[l].end());
        
        int qlen = (qend > qstart ? qend-qstart+1 : qend+interval-qstart);
        double opacity = 0, fadeStep = 1. / (double)qlen;
        for (int c = 0; c < qlen && qlen > interval / 2; ++c)
        {
            int qid = (c + qstart) % interval;
            opacity += fadeStep;
//            cout << "opacity: " << opacity << endl;
            for (size_t l = 0; l < nlayer; ++l)
                for (list<SplineStroke>::iterator iter = add_queue[l].begin(); iter != add_queue[l].end(); iter++)
                {
                    iter->setAlpha(opacity);
                    strokes_queue[qid][l].push_front(*iter);
                }
        }
        
        // dump
        if (qstart == qend)
        {
            cout << "dumping qstart = " << qstart << endl;
            dumpedFrame = frame - interval;
            cout << "write to frame " << dumpedFrame << endl;
            sprintf(buf, im, dumpedFrame);
            src = imread(buf);
            ps.setSourceImage(src);
            ps.paintOverRef(canvas, strokes_queue[qstart], nlayer);
            blur(canvas, smooth, Size(3, 3));

            sprintf(buf, out, "out", dumpedFrame);
            imwrite(buf, canvas);
            ps.fixEdges(canvas);
            sprintf(buf, out, "out-edge", dumpedFrame);
            imwrite(buf, canvas);
            sprintf(buf, out, "out-smooth", dumpedFrame);
            imwrite(buf, smooth);
            ps.fixEdges(smooth);
            sprintf(buf, out, "out-smooth-edge", dumpedFrame);
            imwrite(buf, smooth);
            qstart = (qstart + 1) % interval;
        }

        for (size_t l = 0; l < nlayer; ++l)
            strokes_queue[qend][l] = new_queue[l];
    }

    for (frame = dumpedFrame+1; frame <= frameEnd; ++frame)
    {
        cout << "dumping cache " << frame << endl;
        sprintf(buf, im, frame);
        src = imread(buf);
        ps.setSourceImage(src);
        ps.paintOverRef(canvas, strokes_queue[qstart], nlayer);
        blur(canvas, smooth, Size(3, 3));

        sprintf(buf, out, "out", frame);
        imwrite(buf, canvas);
        ps.fixEdges(canvas);
        sprintf(buf, out, "out-edge", frame);
        imwrite(buf, canvas);
        sprintf(buf, out, "out-smooth", frame);
        imwrite(buf, smooth);
        ps.fixEdges(smooth);
        sprintf(buf, out, "out-smooth-edge", frame);
        imwrite(buf, smooth);

        qstart = (qstart + 1) % interval;
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

void propagate(list<SplineStroke> *next_queue, const list<SplineStroke> *cur_queue,
               size_t nlayer, const Mat &u, const Mat &v)
{
    assert(u.type() == CV_64F);

    cout << "propagating strokes..." << endl;

    int offset, step = u.step1(), width = u.cols, height = u.rows;
    double *pu = (double*)u.data, *pv = (double*)v.data;
    for (size_t i = 0; i < nlayer; ++i)
    {
        next_queue[i].clear();
        for (list<SplineStroke>::const_iterator iter = cur_queue[i].begin(); iter != cur_queue[i].end(); ++iter)
        {
            SplineStroke ss;
            ss.set(iter->getRadius(), iter->getAngle(), iter->getColor());
            ss.setAlpha(iter->getAlpha());
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);
                if (p.y >= height || p.y < 0 || p.x >= width || p.x < 0) continue;

                offset = p.y * step + p.x;
                p.x += cvRound(pu[offset]);
                p.y += cvRound(pv[offset]);
                if (p.y >= height || p.y < 0 || p.x >= width || p.x < 0) continue;

                ss.add(p);
            }

            if (ss.nPoints() > 0) next_queue[i].push_back(ss);
        }
//        cout << "queue size: " << next_queue[i].size() << endl;
    }

    cout << "done" << endl;
}

void updateStrokes(list<SplineStroke> *queue, size_t nlayer, const Mat &src)
{
    const double factor = 0.1;
    for (size_t l = 0; l < nlayer; ++l)
        for (list<SplineStroke>::iterator iter = queue[l].begin(); iter != queue[l].end(); ++iter)
        {
            Vec3i color(0, 0, 0);
            Vec3i bcolor = iter->getColor();
            for (int k = 0; k < iter->nPoints(); ++k)
            {
                Point p = iter->get(k);
                color += src.at<Vec3b>(p.y, p.x);
            }

            color /= iter->nPoints();
            iter->setColor(factor * color + (1-factor) * bcolor);
        }
}

void removeStrokes(list<SplineStroke> *strokes_queue, size_t nlayer, const Mat &src)
{
    assert(src.type() == CV_8UC3);
    cout << "removing strokes..." << endl;
    
    const int threshold = 90 * 3;
    const double fadeStep = 0.1;
    int width = src.cols, height = src.rows;
    for (size_t i = 0; i < nlayer; ++i)
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
            if (pointDiff >= 0.5 * iter->nPoints()) iter->changeOpacity(-fadeStep);
            if (pointDiff < 0.1 * iter->nPoints())  iter->changeOpacity(fadeStep);
            
            if (iter->isTransparent() || iter->nPoints() == 0) iter = strokes_queue[i].erase(iter);
            else iter++;
        }
    }

    cout << "done" << endl;
}

// remove strokes by distance
void removeStrokes2(list<SplineStroke> *strokes_queue, size_t nlayer, const Mat &src)
{
    assert(src.type() == CV_8UC3);
    cout << "removing strokes..." << endl;

    const double fadeStep = 0.1;
    int i, j, hit;
    vector<bool> taken;
    list<SplineStroke>::iterator iter, iter2;
    for (size_t l = 0; l < nlayer; ++l)
    {
        // find out cllopase strokes
        taken.resize(strokes_queue[l].size(), false);
        for (iter = strokes_queue[l].begin(), i = 0; iter != strokes_queue[l].end(); iter++, i++)
        {
            if (iter->nPoints() == 0) continue;
            
            int threshold = iter->getRadius() * iter->getRadius() / 2;
            Point p1 = iter->get(0);
            iter2 = iter;
            for (iter2++, j = i+1, hit = 0; iter2 != strokes_queue[l].end(); iter2++, j++)
            {
                if (iter2->nPoints() == 0) continue;

                Point p2 = iter2->get(0);
                int dist = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
                if (dist < threshold && !taken[j]){
                    taken[j] = true;
                    iter2->changeOpacity(fadeStep);
                    hit++;
                }
            }

            if (hit > 0) iter->changeOpacity(fadeStep);
        }

        // remove
        iter = strokes_queue[l].begin();
        while (iter != strokes_queue[l].end()){
            if (iter->isTransparent() || iter->nPoints() == 0) iter = strokes_queue[l].erase(iter);
            else iter++;
        }
    }

    cout << "done" << endl;
}

void addStrokes(list<SplineStroke> *new_queue, const list<SplineStroke> *strokes_queue,
                const size_t nlayer, Mat &src)
{
    cout << "adding strokes..." << endl;

    int width = src.cols, height = src.rows;
    Mat canvas = Mat::zeros(height, width, CV_8UC3);
    PainterlyService ps;
    ps.setSourceImage(src);
    ps.paint_layer(canvas, strokes_queue, nlayer);

    Mat gray, smooth, mask;
    cvtColor(canvas, gray, CV_BGR2GRAY);
    smooth = (gray > 0);
    morphologyEx(smooth, mask, MORPH_CLOSE, Mat(5, 5, CV_8U), Point(-1, -1), 3);

#ifdef DEBUG
    char buf[256];
    sprintf(buf, out, "add-mask", frame);
    imwrite(buf, mask);
#endif

    list<SplineStroke> tmp_queue[nlayer];
    list<SplineStroke>::iterator iter;
    canvas.setTo(0);
    src.copyTo(canvas, mask);
    ps.render(canvas, tmp_queue, nlayer);
    for (size_t l = 0; l < nlayer; ++l)
    {
        new_queue[l].clear();
        for (iter = tmp_queue[l].begin(); iter != tmp_queue[l].end(); iter++)
        {
            Mat taken = Mat::zeros(height, width, CV_8U);
            int area = 0, covered = 0, R = iter->getRadius();
            for (int i = 0; i < iter->nPoints(); ++i)
            {
                Point p = iter->get(i);
                for (int h = -R; h <= R; ++h)
                    for (int w = -R; w <= R; ++w)
                    {
                        int hh = p.y + h, ww = p.x + w;
                        if (ww < 0 || hh >= height || ww >= width || hh < 0) continue;

                        if (taken.at<uchar>(hh, ww) != 0) continue;

                        area++;
                        if (mask.at<uchar>(hh, ww) == 0) covered++;
                        taken.at<uchar>(hh, ww) = 1;
                    }
            }

            if (covered >= area * 0.5) new_queue[l].push_back(*iter);
        }
    }
    
    cout << "done" << endl;
}
